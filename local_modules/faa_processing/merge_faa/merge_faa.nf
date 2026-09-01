#!/usr/bin/bash Nextflow




process MERGE_FAA {
    executor 'local'

    publishDir "${params.output_dir}/faa", mode: "copy"

    input:
    path(faa_files)

    output:
    path "faa", emit: faa_dir

    script:
    """
    mv *.faa .
    """

    stub:
    """
    touch temp.faa
    """
}

workflow {
    ch_faa_list = BAKTA_BAKTA.out.faa.map { meta, faa -> faa }.collect()
    
    ch_faa_dir = MERGE_FAA(ch_faa_list) 
}
