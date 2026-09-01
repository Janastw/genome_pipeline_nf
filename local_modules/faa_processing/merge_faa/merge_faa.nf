#!/usr/bin/env nextflow

process MERGE_FAA {
    executor 'local'

    input:
    path("faa/*") 

    output:
    path "faa", emit: faa_dir

    script:
    """
    echo "Successfully bundled .faa files"
    """

    stub:
    """
    mkdir -p faa
    touch faa/stub_file.faa
    """
}

workflow {
    ch_faa_list = BAKTA_BAKTA.out.faa.map { meta, faa -> faa }.collect()
    
    ch_faa_dir = MERGE_FAA(ch_faa_list) 
}
