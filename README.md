# genome_pipeline_nf

A Nextflow DSL2 pipeline for prokaryotic genome QC, annotation, and functional profiling. Designed for HPC execution on Midway3 (University of Chicago RCC) via Slurm + Conda/Singularity.

## Pipeline overview

```
Input TSV manifest (sample_id, assembly_path)
  |
  +-- GUNC_RUN              Contamination screening
  +-- CHECKM2_PREDICT       Completeness / contamination
  |
  +-- QC_TRIAGE             Filter: completeness >= 90%, contamination <= 5%, pass GUNC
  +-- DUPLICATE_HANDLING    Dereplicate by priority TSV
  |
  +-- FASTANI               All-vs-all ANI matrix
  |
  +-- BAKTA_BAKTA           Genome annotation
       |
       +-- EXTRACT_16S_RRNA    16S rRNA extraction (barrnap + bedtools)
       +-- EGGNOGMAPPER         Functional annotation (COG / eggNOG)
       +-- KOFAMSCAN            KEGG KO assignment
            |
            +-- KO_COUNT_MATRIX  Per-genome KO count matrix
  |
  +-- COUNT_AMINO_ACIDS     Amino acid composition (counts, proportions, CLR)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.10.1
- Conda (for local module environments) or Singularity/Apptainer (for nf-core modules on HPC)
- Slurm (for HPC execution)

## Quick start

### Run directly from GitHub

```bash
nextflow run Janastw/genome_pipeline_nf \
  -c /path/to/your/midway3.config \
  --genome_manifest samples.tsv \
  --bakta_db /path/to/bakta_db \
  --kofamscan_profiles /path/to/profiles.tar.gz \
  --kofamscan_ko_list /path/to/ko_list.gz \
  --eggnogmapper_db /path/to/eggnog_db \
  --priority_tsv priority.tsv
```

### Run from a local clone

```bash
git clone https://github.com/Janastw/genome_pipeline_nf.git
cd genome_pipeline_nf

# Copy and fill in your HPC config
cp configs/midway3.config.template configs/midway3.config
# Edit configs/midway3.config with your SLURM account, scratch path, etc.

nextflow run main.nf \
  -c configs/midway3.config \
  --genome_manifest samples.tsv \
  --bakta_db /path/to/bakta_db \
  --kofamscan_profiles /path/to/profiles.tar.gz \
  --kofamscan_ko_list /path/to/ko_list.gz \
  --eggnogmapper_db /path/to/eggnog_db \
  --priority_tsv priority.tsv
```

## Input format

### Genome manifest (`--genome_manifest`)

Tab-separated file with a header row:

```
sample_id	assembly_path
genome_001	/path/to/genome_001.fna
genome_002	/path/to/genome_002.fna
```

All `assembly_path` values must be absolute paths accessible from compute nodes.

### Priority TSV (`--priority_tsv`)

Tab-separated file used to break ties when deduplicating genomes with the same `sample_id`. Must have at minimum an `assembly_path` column and a `priority` column (integer; higher = preferred).

```
assembly_path	priority
/path/to/genome_001.fna	2
/path/to/genome_002.fna	1
```

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `--genome_manifest` | Yes | Path to input TSV manifest |
| `--bakta_db` | Yes | Path to Bakta database directory |
| `--kofamscan_profiles` | Yes | Path to KOfam HMM profiles archive |
| `--kofamscan_ko_list` | Yes | Path to KOfam KO metadata file |
| `--eggnogmapper_db` | Yes | Path to eggNOG-mapper database directory (must contain `eggnog_proteins.dmnd`) |
| `--priority_tsv` | Yes | Path to priority TSV for deduplication |
| `--gunc_db` | No | Path to existing GUNC ProGenomes DB (auto-downloaded if omitted) |
| `--checkm2_db` | No | Path to existing CheckM2 DB (auto-downloaded if omitted) |
| `--outdir` | No | Output directory (default: `results/`) |

## HPC setup (Midway3)

1. Copy the config template and fill in your values:
   ```bash
   cp configs/midway3.config.template configs/midway3.config
   ```
2. Edit `configs/midway3.config`:
   - Replace `YOUR_SLURM_ACCOUNT` with your SLURM account name
   - Replace `YOUR_PARTITION` with your Slurm partition name
   - Replace `YOUR_CNETID` with your CNetID
   - Set `cacheDir` to your Singularity image cache path
3. Pass it at runtime with `-c configs/midway3.config`. Both `nextflow.config` and `configs/midway3.config` are gitignored - never commit them.

## Outputs

```
results/
  gunc/                  GUNC contamination TSVs per genome
  checkm2/               CheckM2 quality TSVs per genome
  qc_triage/             QC summary and filtered manifest
  duplicate_handling/    Deduplication report and cleaned manifest
  fastani/               All-vs-all ANI matrix
  bakta/<sample_id>/     Annotation files (.gff, .fna, .faa, .tsv, ...)
  eggnogmapper/          eggNOG functional annotations
  kofamscan/             KOfam assignment files per genome
  ko_count_matrix/       ko_count_matrix.tsv
  amino_acid_counts/     amino_acid_counts.csv, proportions, CLR
  extract_16s_rrna/      16S rRNA sequences per genome
```

## Development

### Adding nf-core modules

```bash
nf-core modules install <module/name>
```

### Local modules

Custom modules live under `local_modules/`. Each subdirectory should contain a `.nf` process file and an `environment.yml` for Conda.

### CI

GitHub Actions runs on every push and pull request to `main`. It checks:
- All `include` paths in `main.nf` resolve to existing files
- All modules listed in `modules.json` have a corresponding `main.nf`
- Local module files contain valid process or workflow blocks
