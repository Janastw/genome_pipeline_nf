#!/usr/bin/env nextflow

process QC_TRIAGE {
    container 'quay.io/hdc-workflows/python-pandas:v1.2.1_latest'

    input:
    path genome_manifest
    path all_gunc_summary
    path all_checkm2_summary

    output:
    path "1_complete_summary.tsv"
    path "2_filtered_summary.tsv"
    path "3_filtered_genome_manifest.tsv" , emit: filtered_genome_manifest

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    genome_manifest = pd.read_csv("${genome_manifest}", sep=chr(9))
    checkm2_df = pd.read_csv("${all_checkm2_summary}", sep=chr(9))
    gunc_df = pd.read_csv("${all_gunc_summary}", sep=chr(9))
    gunc_df["genome"] = gunc_df["genome"].str.rsplit(".", n=1).str[-2]
    genome_manifest["filename"] = genome_manifest["assembly_path"].str.rsplit("/", n=1).str[-1].str.rsplit(".", n=1).str[0]


    checkm2_df.columns = checkm2_df.columns.str.replace(r'^([A-Z])', r'checkm2_\\1', regex=True)
    checkm2_df.columns = checkm2_df.columns.str.lower()

    gunc_cols = ['n_genes_called', 'n_genes_mapped', 'n_contigs', 'taxonomic_level',
                'proportion_genes_retained_in_major_clades', 'genes_retained_index',
                'clade_separation_score', 'contamination_portion', 'n_effective_surplus_clades',
                'mean_hit_identity', 'reference_representation_score', 'pass.GUNC']

    gunc_df = gunc_df.rename(columns={col: f'gunc_{col}' for col in gunc_cols})

    summary_df = checkm2_df.merge(gunc_df, left_on="checkm2_name", right_on="genome", how="inner").merge(genome_manifest, left_on="checkm2_name", right_on="filename", how="inner")

    summary_df = summary_df[['sample_id', 'assembly_path', 'checkm2_completeness', 'checkm2_contamination', 
                'checkm2_completeness_model_used', 'checkm2_translation_table_used', 'checkm2_coding_density',
                'checkm2_contig_n50', 'checkm2_average_gene_length', 'checkm2_genome_size', 'checkm2_gc_content',
                'checkm2_total_coding_sequences', 'checkm2_total_contigs', 'checkm2_max_contig_length',
                'checkm2_additional_notes', 'gunc_n_genes_called', 'gunc_n_genes_mapped', 'gunc_n_contigs',
                'gunc_taxonomic_level', 'gunc_proportion_genes_retained_in_major_clades', 'gunc_genes_retained_index',
                'gunc_clade_separation_score', 'gunc_contamination_portion', 'gunc_n_effective_surplus_clades',
                'gunc_mean_hit_identity', 'gunc_reference_representation_score', 'gunc_pass.GUNC']]

    summary_df.to_csv("1_complete_summary.tsv", sep='\t', index=False)
    filtered_summary_df = summary_df[(summary_df["checkm2_completeness"] >= 90) & (summary_df["checkm2_contamination"] <= 5) & (summary_df["gunc_pass.GUNC"] == True)]
    filtered_summary_df.to_csv("2_filtered_summary.tsv", sep='\t', index=False)
    filtered_genome_manifest_df = filtered_summary_df[["sample_id", "assembly_path"]]
    filtered_genome_manifest_df.to_csv("3_filtered_genome_manifest.tsv", sep='\t', index=False)
    """
}