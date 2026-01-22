# BOLD Gene Extractor

Filters BOLD (Barcode of Life Data System) TSV data files to extract rows matching specified marker genes, producing a TSV output suitable for downstream bioinformatics analysis.

## Overview

This script parses BOLD TSV files (both web downloads and data package exports) and extracts:
- All rows where the `marker_code` column matches one or more specified gene names
- Header row is preserved in output
- Case-insensitive gene matching
- Support for multiple genes in a single run

## Requirements

- Python 3.6+
- No external dependencies (uses standard library only)

## Installation

```bash
# Make executable on Unix systems
chmod +x bold_gene_extract.py
```

## Usage

### Basic Usage

```bash
# Extract rbcL sequences
python bold_gene_extract.py -g rbcL input.tsv output.tsv

# Extract multiple related genes (comma-separated)
python bold_gene_extract.py -g rbcL,rbcLa input.tsv output.tsv

# Extract multiple genes (multiple -g flags)
python bold_gene_extract.py -g rbcL -g matK -g ITS2 input.tsv output.tsv
```

### With Verbose Output

```bash
python bold_gene_extract.py -v -g rbcL,rbcLa large_dataset.tsv rbcL_only.tsv
```

### Gene Name Matching

Gene matching is case-insensitive:

```bash
# These are all equivalent
python bold_gene_extract.py -g rbcL input.tsv output.tsv
python bold_gene_extract.py -g RBCL input.tsv output.tsv
python bold_gene_extract.py -g RbcL input.tsv output.tsv
```

### Pipeline Integration (stdin/stdout)

```bash
# From gzipped file
zcat input.tsv.gz | python bold_gene_extract.py -g rbcL - - > output.tsv

# Pipe to another tool
python bold_gene_extract.py -g COI-5P input.tsv - | cut -f1,5,70 > subset.tsv
```

## Command Line Arguments

| Argument | Description |
|----------|-------------|
| `input` | Input BOLD TSV file (positional, use "-" for stdin) |
| `output` | Output TSV file (positional, use "-" for stdout) |
| `-g, --gene` | **Required.** Target gene name(s) to extract. Can be specified multiple times or as comma-separated values. Case-insensitive. |
| `-v, --verbose` | Print progress and statistics to stderr |
| `-d, --delimiter` | Field delimiter (default: tab) |

## Input Format

The script supports BOLD TSV files from:
- **Web downloads**: Direct downloads from BOLD website with additional metadata columns
- **Data packages**: Extracted TSV files from BOLD public data packages

Both formats share the core `marker_code` column used for filtering. The web download format includes additional columns (e.g., `geopol_denorm.country_iso3`, `marker_count`, `pre_md5hash`) that are not present in the data package format.

### Key Columns

The script requires the `marker_code` column to be present. Common marker codes include:

| Marker Code | Description |
|-------------|-------------|
| COI-5P | Cytochrome c oxidase subunit I (animal barcode) |
| rbcL | Ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit |
| rbcLa | rbcL variant/amplicon |
| matK | Maturase K |
| ITS | Internal transcribed spacer |
| ITS2 | Internal transcribed spacer 2 |
| trnH-psbA | trnH-psbA intergenic spacer |

## Output Format

Tab-separated file with the same columns as the input file. Only rows matching the specified gene(s) are included. The header row is always preserved.

## Behaviour

### Filtering Logic

- Rows are included if `marker_code` matches any of the target genes (case-insensitive)
- Empty rows are skipped
- Rows with insufficient columns generate a warning (with `-v`) and are skipped
- Header row is always copied to output

### Performance

- Memory-efficient streaming processing (no loading entire file into memory)
- Progress output every 100,000 rows when using `-v`
- Suitable for files with millions of rows

## HPC Cluster Usage

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=bold_extract
#SBATCH --output=bold_extract_%j.out
#SBATCH --error=bold_extract_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load python/3.9

python bold_gene_extract.py \
    -v \
    -g rbcL,rbcLa \
    /data/bold_sequences.tsv \
    /results/rbcL_extracted.tsv
```

### PBS/Torque Example

```bash
#!/bin/bash
#PBS -N bold_extract
#PBS -l walltime=01:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
module load python/3.9

python bold_gene_extract.py \
    -v \
    -g COI-5P \
    bold_animals.tsv \
    COI_data.tsv
```

### Snakemake Integration

```python
rule extract_gene:
    input:
        "data/raw/{dataset}.tsv"
    output:
        "data/processed/{dataset}_{gene}.tsv"
    params:
        gene="{gene}"
    log:
        "logs/extract_{dataset}_{gene}.log"
    shell:
        """
        python scripts/bold_gene_extract.py \
            -v \
            -g {params.gene} \
            {input} \
            {output} \
            2> {log}
        """

# Example usage in workflow
rule all:
    input:
        expand("data/processed/{dataset}_{gene}.tsv",
               dataset=["bold_plants", "bold_animals"],
               gene=["rbcL", "matK", "COI-5P", "ITS2"])
```

### Nextflow Integration

```nextflow
process extractGene {
    input:
    path tsv
    val gene

    output:
    path "${tsv.baseName}_${gene}.tsv"

    script:
    """
    python bold_gene_extract.py \
        -v \
        -g ${gene} \
        ${tsv} \
        ${tsv.baseName}_${gene}.tsv
    """
}
```

## Performance Notes

- Streaming processing: memory usage is minimal regardless of input file size
- Typical throughput: 1-5 million rows per minute depending on I/O
- For very large files, consider splitting input first if parallel processing is needed
- Progress output every 100,000 rows when using `-v`

## Error Handling

- Missing `marker_code` column generates a clear error with available column names
- Empty input files are detected and reported
- Rows with insufficient columns are skipped with optional warning

## Examples

### Extract plant barcodes

```bash
# Extract rbcL and matK for plants
python bold_gene_extract.py -v -g rbcL,rbcLa,matK \
    bold_plants_complete.tsv \
    bold_plants_barcodes.tsv
```

### Extract animal COI barcodes

```bash
# Extract COI-5P sequences
python bold_gene_extract.py -v -g COI-5P \
    bold_animals.tsv \
    bold_COI.tsv
```

### Process multiple markers separately

```bash
for gene in rbcL matK ITS2 trnH-psbA; do
    python bold_gene_extract.py -v -g $gene \
        sequences.tsv \
        ${gene}_extracted.tsv
done
```

### Filter and count unique species

```bash
# Extract rbcL and count unique species
python bold_gene_extract.py -g rbcL input.tsv - | \
    cut -f22 | sort | uniq -c | sort -rn > species_counts.txt
```

## Comparison with NCBI GenBank Extractor

This tool complements the `ncbi_gb_extract.py` script for GenBank files:

| Feature | bold_gene_extract.py | ncbi_gb_extract.py |
|---------|---------------------|-------------------|
| Input format | BOLD TSV | GenBank flat file |
| Gene matching | marker_code column | CDS gene qualifier |
| Processing | Row filtering | Record parsing |
| Memory usage | Streaming (low) | All records in memory |
| Dependencies | None | BioPython |

## License

MIT

## Author

Ben Price, Natural History Museum London
