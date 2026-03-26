# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

snoDB3 is a Snakemake-based bioinformatics pipeline that processes small nucleolar RNA (snoRNA) data and generates PostgreSQL database tables for the snoDB3 web database. The pipeline handles snoRNA identification, expression data, host gene relationships, rRNA modifications, and ENCODE eCLIP data.

## Commands

### Running the Pipeline

```bash
# Run complete pipeline (create conda environments automatically)
snakemake --use-conda

# Run with specific number of cores
snakemake --use-conda -j <n_cores>

# Dry run to check what will be executed
snakemake --dry-run

# Run a specific rule
snakemake -j1 <rule_name> --use-conda
```

### Database Setup

The pipeline uses Docker with the `psql13` image for PostgreSQL operations. Ensure Docker is running before executing database-related rules.

## Architecture

### Workflow Structure

```
Snakefile              # Main entry point, defines target outputs and includes rule files
config.yaml            # All paths and file names configuration
rules/
  preprocessing.smk    # Rules for adding new snoRNAs and updating expression data
  psql.smk             # Rules for generating PostgreSQL database tables
scripts/               # Python scripts executed by Snakemake rules
envs/
  python.yaml          # Conda environment specification
```

### Key Data Flow

1. **Preprocessing** (`rules/preprocessing.smk`):
   - `add_alphonse_snoRNA`: Adds new snoRNA entries from BED file
   - `update_snorna_expression`: Updates expression matrix with new data

2. **PostgreSQL Table Generation** (`rules/psql.smk`):
   - Each rule generates a `data_table.tsv` and `data_script.sql` in `data/psql/<table_name>/`
   - Scripts generate SQL with: schema selection, table drop, table creation, data import
   - Docker containers mount the data directory and execute SQL scripts

### PostgreSQL Tables Generated

- `external_ids`: Cross-references to external databases (Ensembl, RefSeq, HGNC, etc.)
- `genomic_location`: Chromosome coordinates
- `basic_features`: Conservation scores, sequences, expression levels
- `specie`: Species information
- `host_features`: Host gene relationships
- `snoRNA_expression` / `host_expression`: Expression matrices
- `targets` / `target_grouped_by`: rRNA modification targets
- `encode_eclip`: ENCODE eCLIP data
- `snoRNA_boxes`: C/D and H/ACA box annotations (cd_data_table, haca_data_table, sca_data_table)
- `rRNAs`, `rRNA_modifications`, `conversion_18S`, `conversion_28S`: rRNA modification data

### Script Pattern

Python scripts in `scripts/` follow a consistent pattern:
1. Read input files via `snakemake.input`
2. Process data using pandas
3. Generate `data_table.tsv` (tab-separated, no header, no index)
4. Generate `data_script.sql` with DROP TABLE, CREATE TABLE, and `\copy` commands
5. Copy `psql_container.sh` to output directory
6. Execute Docker container to load data into PostgreSQL

### Conda Environment

The environment (`envs/python.yaml`) includes: matplotlib, numpy, pandas, seaborn, requests, beautifulsoup4, lxml, pybedtools. Snakemake automatically creates environments in `.snakemake/conda/`.

## Configuration

All file paths are centralized in `config.yaml`. When adding new data sources:
1. Add path in `path:` section
2. Add filename in appropriate section
3. Reference in rules as `config['path']['<path_key>']` and `config['<section>']['<file_key>']`

## Database Schema

The database schema is documented in `database_organisation/snoDB3_schema.svg` (or `.png`).
