# Plan: Extract phastCons Scores for snoRNA Regions via bigWigAverageOverBed

## Context

Need average phastCons100way conservation scores for 56 snoRNA regions. The bigWig file (~50GB if converted to bedGraph) can be queried **directly** using `bigWigAverageOverBed` from UCSC tools — no conversion needed. This will be integrated into an existing Snakemake workflow in a separate project folder.

## Tool

**`bigWigAverageOverBed`** — UCSC utility that takes a bigWig + BED file and outputs mean (and other stats) per region. Reads only the needed portions of the bigWig, so it's fast.

## Inputs

- **bigWig file**: `hg38.phastCons100way.bw` (phastCons scores)
- **BED file**: regular BED file with 56 snoRNA genomic coordinates

## Output

- Tab-separated file with one row per snoRNA region containing: `chrom, chromStart, chromEnd, name, size, covered, sum, mean0, mean`

## Snakemake Integration Plan

### 1. Conda environment file (`envs/phastcons.yaml`)

```yaml
channels:
  - bioconda
  - conda-forge
dependencies:
  - ucsc-bigwigtobedgraph  # provides bigWigAverageOverBed
```

### 2. Snakemake rule

```python
rule phastcons_average:
    input:
        bigwig="path/to/hg38.phastCons100way.bw",
        bed="path/to/snornas.bed"
    output:
        "results/phastcons/snornas_phastcons_avg.tsv"
    conda:
        "envs/phastcons.yaml"
    shell:
        "bigWigAverageOverBed {input.bigwig} {input.bed} {output}"
```

### 3. Add to `rule all`

Add the output path to the `rule all` input list so the pipeline knows to produce it.

## Verification

1. Activate the conda env manually: `conda env create -f envs/phathcons.yaml`
2. Run the command by hand: `bigWigAverageOverBed hg38.phastCons100way.bw snornas.bed test_output.tsv`
3. Check output has 56 rows (one per snoRNA) with reasonable mean values (0-1 range for phastCons)
4. Run the full Snakemake pipeline and confirm the rule executes correctly
