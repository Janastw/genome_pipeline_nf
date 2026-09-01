#!/usr/bin/env nextflow

include { GUNC_DOWNLOADDB }          from './modules/nf-core/gunc/downloaddb/main.nf'
include { GUNC_RUN }                 from './modules/nf-core/gunc/run/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from './modules/nf-core/checkm2/databasedownload/main.nf'
include { CHECKM2_PREDICT }          from './modules/nf-core/checkm2/predict/main.nf'
include { FASTANI }                  from './modules/nf-core/fastani/main.nf'
include { BAKTA_BAKTA }              from './modules/nf-core/bakta/bakta/main.nf'
include { KOFAMSCAN }                from './modules/nf-core/kofamscan/main.nf'
include { EGGNOGMAPPER }             from './modules/nf-core/eggnogmapper/main.nf'

include { QC_TRIAGE }          from './local_modules/qc/qc_triage/main.nf'
include { DUPLICATE_HANDLING } from './local_modules/qc/duplicate_handling/main.nf'
include { KO_COUNT_MATRIX }    from './local_modules/faa_processing/create_ko_matrix/ko_count_matrix.nf'
include { COUNT_AMINO_ACIDS }  from './local_modules/faa_processing/count_amino_acids/count_aa.nf'
include { EXTRACT_16S_RRNA }   from './local_modules/faa_processing/extract_16s_rrna/extract_16s_rrna.nf'
include { MERGE_FAA }          from './local_modules/faa_processing/merge_faa/merge_faa.nf'


workflow {
    genome_manifest = params.genome_manifest

    samples_ch = Channel
        .fromPath(genome_manifest)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = [ id: row.sample_id ]
            def fasta = file(row.assembly_path)
            if ( !fasta.exists() ) {
                throw new IllegalArgumentException("Fasta file not found: ${fasta}")
            }
            tuple(meta, fasta)
        }

    // --- QC ---

    if (params.gunc_db) {
        gunc_db_ch = Channel.fromPath(params.gunc_db)
    } else {
        gunc_db_ch = GUNC_DOWNLOADDB("progenomes_2.1").db
    }

    gunc_run_ch = GUNC_RUN(
        samples_ch,
        gunc_db_ch
    ).maxcss_level_tsv

    if (params.checkm2_db) {
        checkm2_db_ch = Channel
            .fromPath(params.checkm2_db)
            .map { path -> tuple([id: "checkm2_db"], path) }.first()
    } else {
        checkm2_db_ch = CHECKM2_DATABASEDOWNLOAD("").database
    }

    checkm2_predict_ch = CHECKM2_PREDICT(
        samples_ch,
        checkm2_db_ch
    ).checkm2_tsv

    all_gunc_ch = gunc_run_ch.map { it[1] }
        .collectFile(name: "all_gunc_summary.tsv", keepHeader: true)

    all_checkm2_ch = checkm2_predict_ch.map { it[1] }
        .collectFile(name: "all_checkm2_summary.tsv", keepHeader: true)

    qc_summary = QC_TRIAGE(
        genome_manifest,
        all_gunc_ch,
        all_checkm2_ch,
    )

    duplicate_genome_manifest = DUPLICATE_HANDLING(
        qc_summary.filtered_genome_manifest,
        params.priority_tsv
    ).filtered_deduplicated_genome_manifest

    qc_samples_ch = duplicate_genome_manifest
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = [ id: row.sample_id ]
            def fasta = file(row.assembly_path)
            if ( !fasta.exists() ) {
                throw new IllegalArgumentException("Fasta file not found: ${fasta}")
            }
            tuple(meta, fasta)
        }

    // --- ANI ---

    query_list_ch = qc_samples_ch
        .map { meta, fasta -> fasta.toString() }
        .collectFile(name: 'query_pathways.txt', newLine: true)

    ref_list_ch = qc_samples_ch
        .map { meta, fasta -> fasta.toString() }
        .collectFile(name: 'ref_pathways.txt', newLine: true)

    FASTANI(
        [[], []],
        [[], []],
        query_list_ch,
        ref_list_ch
    )

    // --- Annotation ---

    db_ch = Channel.value(file(params.bakta_db))

    BAKTA_BAKTA(
        qc_samples_ch,
        db_ch,
        [],
        [],
        [],
        [],
    )

    ch_faa_dir = MERGE_FAA(BAKTA_BAKTA.out.faa.map { it[1] }.collect())

    COUNT_AMINO_ACIDS(
        file("${projectDir}/local_modules/faa_processing/count_amino_acids/count_aa.py"),
        BAKTA_BAKTA.out.faa.map { it[1] }.collect()
    )

    // --- Functional profiling ---

    kofamscan_profiles_ch = Channel.value(file(params.kofamscan_profiles))
    kofamscan_ko_list_ch  = Channel.value(file(params.kofamscan_ko_list))

    KOFAMSCAN(
        BAKTA_BAKTA.out.faa,
        kofamscan_profiles_ch,
        kofamscan_ko_list_ch
    )

    KO_COUNT_MATRIX(
        file("${projectDir}/local_modules/faa_processing/create_ko_matrix/create_ko_count_matrix.py"),
        KOFAMSCAN.out.txt.map { it[1] }.collect()
    )

    eggnogmapper_db_ch = Channel.value(file(params.eggnogmapper_db))
        .map { path -> tuple("diamond", path.resolve("eggnog_proteins.dmnd")) }
        .first()

    EGGNOGMAPPER(
        BAKTA_BAKTA.out.faa,
        eggnogmapper_db_ch,
        file(params.eggnogmapper_db)
    )

    // --- 16S extraction ---

    EXTRACT_16S_RRNA(
        BAKTA_BAKTA.out.fna.join(BAKTA_BAKTA.out.gff)
    )

}
