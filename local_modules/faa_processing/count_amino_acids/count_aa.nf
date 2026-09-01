// Taken from BaGPipe
// Annotate fasta genome assemblies using Prokka
process COUNT_AMINO_ACIDS {
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/qiime2/core:2023.7'
        : 'quay.io/qiime2/core:2023.7'}"


    // publishDir "${params.output_dir}/amino_acid_counts", mode:'copy', overwrite: true

    input:
    path(count_aa_script)
    path(faa_files)

    output:
    path "amino_acid_counts.csv", emit: aa_counts_csv
    path "amino_acid_proportions.csv", emit: aa_proportions_csv
    path "amino_acid_clr.csv", emit: aa_clr_csv

    script:
    """
    python3 ${count_aa_script} ${faa_files} -o amino_acid_counts.csv

    python3 -c "import pandas as pd; \
        df = pd.read_csv('amino_acid_counts.csv', index_col=0); \
        df = df.div(df.sum(axis=1), axis=0); \
        df.to_csv('amino_acid_proportions.csv'); \
        from skbio.stats.composition import clr; \
        df_clr = pd.DataFrame(clr(df), index=df.index, columns=df.columns); \
        df_clr.to_csv('amino_acid_clr.csv')"
    """

    stub:
    """
    touch amino_acid_counts.csv amino_acid_proportions.csv amino_acid_clr.csv
    """
}

workflow RUN_COUNT_AMINO_ACIDS {

    take:
        count_aa_script_ch
        faa_files_ch
    main:
        COUNT_AMINO_ACIDS(count_aa_script_ch, faa_files_ch)
    emit:
        aa_counts_csv = COUNT_AMINO_ACIDS.out.aa_counts_csv
        aa_proportions_csv = COUNT_AMINO_ACIDS.out.aa_proportions_csv
        aa_clr_csv = COUNT_AMINO_ACIDS.out.aa_clr_csv
}

workflow {
    count_aa_script_ch = Channel.fromPath("${moduleDir}/count_aa.py")
    faa_files_ch = Channel.fromPath("${params.faa_input_dir}/*.faa").collect()

    COUNT_AMINO_ACIDS(count_aa_script_ch, faa_files_ch)
}



