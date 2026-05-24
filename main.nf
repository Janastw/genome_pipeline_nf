#!/usr/bin/env nextflow

include { GUNC_DOWNLOADDB } from './modules/nf-core/gunc/downloaddb/main.nf'
include { GUNC_RUN } from './modules/nf-core/gunc/run/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from './modules/nf-core/checkm2/databasedownload/main.nf'
include { CHECKM2_PREDICT } from './modules/nf-core/checkm2/predict/main.nf'
include { FASTANI } from './modules/nf-core/fastani/main.nf'
include { BAKTA_BAKTADBDOWNLOAD } from './modules/nf-core/bakta/baktadbdownload/main.nf'
include { BAKTA_BAKTA } from './modules/nf-core/bakta/bakta/main.nf'
include { KOFAMSCAN } from './modules/nf-core/kofamscan/main.nf'

include { QC_TRIAGE } from './modules/qc/qc_triage/main.nf'
include { DUPLICATE_HANDLING } from './modules/qc/duplicate_handling/main.nf'
include { KO_COUNT_MATRIX } from './modules/faa_processing/create_ko_matrix/ko_count_matrix.nf'
include { COUNT_AMINO_ACIDS } from './modules/faa_processing/count_amino_acids/count_aa.nf'


workflow {
    genome_manifest = params.genome_manifest

    samples_ch = Channel
        .fromPath(genome_manifest)
        .splitCsv(header: true, sep: "\t")
        .map { row -> 
            def meta = [ id: row.sample_id ]
            def fasta = file(row.assembly_path)
            // Validates that the assembly file exists before proceeding, needed for downstream processes that rely on the assembly file being present
            if ( !fasta.exists() ) {
                throw new IllegalArgumentException("Fasta file not found: ${fasta}")
            }
            tuple(meta, fasta)
        }

    if (params.gunc_db) {
        gunc_db_ch = Channel
            .fromPath(params.gunc_db)
            .map { path -> tuple([id: "gunc_db"], path) }.first()
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
    
    all_checkm2_ch = checkm2_predict_ch.map { it[1]}
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
            // Validates that the assembly file exists before proceeding, needed for downstream processes that rely on the assembly file being present
            if ( !fasta.exists() ) {
                throw new IllegalArgumentException("Fasta file not found: ${fasta}")
            }
            tuple(meta, fasta)
        }
    
    // qc_samples_list_ch = qc_samples_ch
    //     .map { meta, fasta -> fasta.toString() }
    //     .collectFile(name: 'qc_samples_pathways.txt', newLine: true)

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

    // FASTANI(
    //     [[], []],
    //     [[], []],
    //     qc_samples_list_ch,
    //     qc_samples_list_ch
    // )
    
    // fastani_ch = samples_ch.combine(qc_samples_ch).map { tuple1, tuple2 -> 
    //     def meta1 = tuple1[0]
    //     def meta2 = tuple2[0]
    //     def fasta1 = tuple1[1]
    //     def fasta2 = tuple2[1]
    //     return meta1, fasta1, meta2, fasta2
    // }

    // FASTANI(fastani_ch, [], [])
    // 1. Combine the channel with itself to get every possible pair
    // fastani_ch = samples_ch
    //     .combine(samples_ch)
    //     .map { meta1, fasta1, meta2, fasta2 ->
    //         // Create a new combined meta map so the output file has a unique name
    //         def new_meta = [ id: "${meta1.id}_vs_${meta2.id}" ]
            
    //         // Return: meta, query, reference
    //         tuple(new_meta, fasta1, fasta2)
    //     }

    // // 2. Feed it into your process
    // FASTANI(fastani_ch)

    db = file("/project/cdonnat/shared/databases/bakta_5.1/db")



    // if (params.bakta_db) {
    //     bakta_db_ch = Channel
    //         .fromPath(params.bakta_db)
    //         .map { path -> tuple([id: "bakta_db"], path) }
    //         .first()
    // } else {
    //     bakta_db_ch = BAKTA_BAKTADBDOWNLOAD().db
    // }

    // bakta_db_ch.view()

    db_ch = Channel.value(file(params.bakta_db))

    BAKTA_BAKTA(
        qc_samples_ch,
        db_ch,
        [],
        [],
        [],
        [],
    )

    COUNT_AMINO_ACIDS(
        file("${projectDir}/modules/faa_processing/count_amino_acids/count_aa.py"),
        BAKTA_BAKTA.out.faa.map { it[1] }.collect()
    )


    kofamscan_profiles_ch = Channel.value(file(params.kofamscan_profiles))
    kofamscan_ko_list_ch = Channel.value(file(params.kofamscan_ko_list))

    KOFAMSCAN(
        BAKTA_BAKTA.out.faa,
        kofamscan_profiles_ch,
        kofamscan_ko_list_ch
    )

    KO_COUNT_MATRIX(
         file("${projectDir}/modules/faa_processing/create_ko_matrix/create_ko_count_matrix.py"),
         KOFAMSCAN.out.txt.map{ it[1] }.collect()
    )

}