#!/usr/bin/bash Nextflow




process MERGE_FAA {
    executor 'local'

    publishDir "${params.output_dir}", mode: "copy"

    input:
    path(faa_files)

    output:
    path "faa", emit: faa_dir

    script:
    """
    mkdir -p faa
    mv *.faa faa
    """

    stub:
    """
    mkdir -p faa
    touch temp.faa
    """
}

workflow {
    ch_faa_list = BAKTA_BAKTA.out.faa.map { meta, faa -> faa }.collect()
    
    ch_faa_dir = MERGE_FAA(ch_faa_list) 
}
