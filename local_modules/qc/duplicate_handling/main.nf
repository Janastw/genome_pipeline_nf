#!/usr/bin/env nextflow

process DUPLICATE_HANDLING {
    container 'quay.io/hdc-workflows/python-pandas:v1.2.1_latest'

    input:
    path filtered_genome_manifest
    path priority_tsv

    output:
    path "1_duplicate_genomes.tsv",                     emit: duplicate_genomes
    path "2_removed_genome_manifest.tsv",               emit: removed_genome_manifest
    path "3_filtered_deduplicated_genome_manifest.tsv", emit: filtered_deduplicated_genome_manifest

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    # Load inputs
    manifest_df = pd.read_csv("${filtered_genome_manifest}", sep='\\t')
    priority_df = pd.read_csv("${priority_tsv}", sep='\\t')

    # Left join the priority data onto the filtered genome manifest
    # Assumes a shared column like 'sample_id' or 'assembly_path' exists to join on. 
    merged_df = pd.merge(manifest_df, priority_df, left_on="assembly_path", right_on="assembly_path", how="left")

    # Fill missing values in the priority column with 0 and ensure it is an integer
    # Replace 'priority' with the actual column name in your priority TSV if different
    priority_col = 'priority' 
    merged_df[priority_col] = merged_df[priority_col].fillna(0).astype(int)

    # Identify duplicates (based on sample_id or whichever column defines a duplicate genome)
    duplicate_mask = merged_df.duplicated(subset=["sample_id"], keep=False)
    all_duplicates_df = merged_df[duplicate_mask]
    all_duplicates_df.sort_values(by="sample_id").to_csv("1_duplicate_genomes.tsv", sep='\\t', index=False)

    # Sort by priority descending so the highest number is at the top of each duplicate group
    merged_df = merged_df.sort_values(by=priority_col, ascending=False)

    # Keep the first occurrence (highest priority), drop the rest to deduplicate
    deduplicated_df = merged_df.drop_duplicates(subset=["sample_id"], keep="first")
    deduplicated_df[["sample_id", "assembly_path"]].to_csv("3_filtered_deduplicated_genome_manifest.tsv", sep='\\t', index=False)

    # Track which rows were removed during the deduplication process
    removed_df = merged_df[~merged_df.index.isin(deduplicated_df.index)]
    removed_df.to_csv("2_removed_genome_manifest.tsv", sep='\\t', index=False)
    """
}