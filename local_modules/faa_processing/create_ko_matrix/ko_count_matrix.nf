#!/usr/bin/bash Nextflow




process KO_COUNT_MATRIX {
    // publishDir "${params.output_dir}/ko_count_matrix", mode: "copy"
    conda "${moduleDir}/environment.yml"

    input:
    path(ko_count_matrix_script)
    path(kofamscan_outputs)

    output:
    path "ko_count_matrix.tsv", emit: ko_count_matrix

    script:
    """
    python3 ${ko_count_matrix_script} \\
        -f ${kofamscan_outputs} \\
        -o ko_count_matrix.tsv
    """

    stub:
    """
    touch ko_count_matrix.tsv
    """
}

workflow run_KO_COUNT_MATRIX {
    take:
        ko_count_matrix_script_ch
        kofamscan_outputs_ch
    main:
        KO_COUNT_MATRIX(ko_count_matrix_script_ch, kofamscan_outputs_ch)
    emit:
        ko_count_matrix = KO_COUNT_MATRIX.out.ko_count_matrix
}

workflow {
    ko_count_matrix_script_ch = Channel.fromPath("${moduleDir}/ko_count_matrix.py")
    kofamscan_outputs_ch = Channel.fromPath("${params.kofamscan_output_dir}/*.tsv").collect()

    KO_COUNT_MATRIX(ko_count_matrix_script_ch, kofamscan_outputs_ch)
}
