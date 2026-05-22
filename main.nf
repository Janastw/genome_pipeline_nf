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

    // mergecheckm_ch = gunc_run_ch.map { it[0] }.combine(checkm2_db_ch).map { tuple(it[0], it[1]) }
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

    // FASTANI(
    //     qc_samples_ch.map { it[1] },
    //     qc_samples_ch.map { it[1] },
    //     "",
    //     ""
    // )

    // checkm2_passing = checkm2_predict_ch[1]
    //     .splitCsv(header: true, sep: "\t")
    //     .filter( row -> row.contamination >= 5.0 || row.completeness <= 90.0 )

    // checkm2_summary = checkm2_predict_ch[1]
    //     .splitCsv(header: true, sep: "\t")
    //     .filter( row -> row.contamination >= 5.0 || row.completeness <= 90.0 )

    // GUNC_MERGECHECKM(
    //     gunc_run_ch
    // )

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

    db_ch          = Channel.value(file(params.bakta_db))
    // proteins_ch    = Channel.value(file("${params.bakta_db}/expert-protein-sequences.dmnd"))
    // regions_ch     = Channel.value(file("${params.bakta_db}/ncRNA-regions"))
    // hmms_ch        = Channel.value(file("${params.bakta_db}/pfam"))
    // prodigal_tf_ch = Channel.value([]) 

    BAKTA_BAKTA(
        qc_samples_ch,
        db_ch,
        [],
        [],
        [],
        [],
    )


    kofamscan_profiles_ch = Channel.value(file(params.kofamscan_profiles))
    kofamscan_ko_list_ch = Channel.value(file(params.kofamscan_ko_list))

    KOFAMSCAN(
        BAKTA_BAKTA.out.faa,
        kofamscan_profiles_ch,
        kofamscan_ko_list_ch
    )

}