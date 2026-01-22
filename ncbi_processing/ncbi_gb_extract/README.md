# NCBI GenBank to TSV Extractor

Extracts structured data from NCBI GenBank flat files for a specified target gene, producing a TSV format suitable for downstream bioinformatics analysis.

## Overview

This script parses GenBank-formatted files and extracts:
- Header metadata (LOCUS, DEFINITION, ACCESSION, VERSION, KEYWORDS)
- Organism and taxonomic lineage
- Reference information (authors, title, journal)
- Source feature qualifiers (organism, strain, collection_date, geo_loc_name, etc.)
- CDS feature qualifiers for the target gene (product, protein_id, EC_number, etc.)
- Nucleotide sequence extracted from ORIGIN based on CDS coordinates

## Requirements

- Python 3.6+
- BioPython (`pip install biopython`)

## Installation

```bash
# Install BioPython if not already available
pip install biopython

# Make executable on Unix systems
chmod +x ncbi_gb_extract.py
```

## Usage

### Basic Usage

```bash
# Extract rbcL gene data from a single file
python ncbi_gb_extract.py -g rbcL input.gb output.tsv

# Process all GenBank files in a folder
python ncbi_gb_extract.py -g COI -i /path/to/gb_files/ -o results.tsv
```

### With Verbose Output

```bash
python ncbi_gb_extract.py -v -g matK sequences.gb matK_data.tsv
```

### Gene Name Matching

Gene matching is case-insensitive:

```bash
# These are all equivalent
python ncbi_gb_extract.py -g rbcL input.gb output.tsv
python ncbi_gb_extract.py -g RBCL input.gb output.tsv
python ncbi_gb_extract.py -g RbcL input.gb output.tsv
```

### Pipeline Integration (stdin/stdout)

```bash
# From gzipped file
zcat input.gb.gz | python ncbi_gb_extract.py -g rbcL - - > output.tsv

# Pipe to another tool
python ncbi_gb_extract.py -g COI input.gb - | cut -f1,5,10 > subset.tsv
```

## Command Line Arguments

| Argument | Description |
|----------|-------------|
| `input` | Input GenBank file (positional, use "-" for stdin) |
| `output` | Output TSV file (positional, use "-" for stdout) |
| `-i, --input-dir` | Input directory containing GenBank files (.gb, .gbk, .genbank) |
| `-o, --output-file` | Output TSV file (alternative to positional) |
| `-g, --gene` | **Required.** Target gene name to extract (case-insensitive) |
| `-v, --verbose` | Print progress and skipped records to stderr |

## Output Format

Tab-separated file with dynamically determined columns based on input data.

### Core Columns

| Column | Description | Example |
|--------|-------------|---------|
| locus_name | GenBank locus identifier | PX577658 |
| locus_length | Sequence length in bp | 1168 |
| locus_mol_type | Molecule type | DNA |
| locus_topology | Linear or circular | linear |
| locus_division | GenBank division | PLN |
| locus_date | Submission date | 19-JAN-2026 |
| definition | Sequence description | Klebsormidium sp. strain VKM Al-481... |
| accession | Primary accession number | PX577658 |
| version | Accession with version | PX577658.1 |
| keywords | Associated keywords | . |
| organism | Source organism name | Klebsormidium sp. |
| taxonomy | Full taxonomic lineage | Eukaryota; Viridiplantae; ... |
| ref_authors | Reference authors | Temraleeva,A., Krivina,E. ... |
| ref_title | Reference title | Direct Submission |
| ref_journal | Reference journal | Submitted (13-NOV-2025) ... |

### Source Feature Columns

Prefixed with `source_`:

| Column | Description | Example |
|--------|-------------|---------|
| source_organism | Organism name | Klebsormidium sp. |
| source_organelle | Organelle source | plastid:chloroplast |
| source_mol_type | Molecule type | genomic DNA |
| source_strain | Strain identifier | VKM Al-481 |
| source_db_xref | Database cross-reference | taxon:13781 |
| source_geo_loc_name | Geographic location | Russia: Lyakhovsky Islands... |
| source_collection_date | Collection date | 2024 |

### CDS Feature Columns

Prefixed with `cds_`:

| Column | Description | Example |
|--------|-------------|---------|
| cds_location | Feature coordinates | 1..1168 or complement(100..500) |
| cds_gene | Gene name | rbcL |
| cds_product | Protein product | ribulose-1,5-bisphosphate carboxylase... |
| cds_protein_id | Protein accession | YCM66143.1 |
| cds_codon_start | Reading frame start | 1 |
| cds_transl_table | Translation table | 11 |
| cds_EC_number | Enzyme Commission number | 4.1.1.39 |
| nucleotide_sequence | Extracted DNA sequence | CTAGCTGCATTTCGGATGAC... |

### Multiple CDS Features

If a record contains multiple CDS features for the same gene (e.g., in complete genomes), additional columns are created with numeric suffixes:

- `cds_location`, `cds_location_2`, `cds_location_3`...
- `cds_protein_id`, `cds_protein_id_2`...
- `nucleotide_sequence`, `nucleotide_sequence_2`...

A warning is logged to stderr when this occurs.

## Behaviour

### Record Filtering

- Only records containing a CDS feature matching the target gene are included
- Records without the target gene are skipped and logged to stderr (with `-v`)
- Gene matching is case-insensitive

### Sequence Extraction

- Nucleotide sequences are extracted from ORIGIN based on CDS coordinates
- Handles `complement()` locations (reverse complement)
- Handles `join()` locations (spliced features)
- Handles partial indicators (`<`, `>`)

### Dynamic Columns

- Columns are discovered dynamically from the input data
- New qualifiers encountered in later records are added as new columns
- All records are held in memory until processing is complete to ensure consistent column ordering

## HPC Cluster Usage

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=gb_extract
#SBATCH --output=gb_extract_%j.out
#SBATCH --error=gb_extract_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load python/3.9
module load biopython  # or: pip install --user biopython

python ncbi_gb_extract.py \
    -v \
    -g rbcL \
    -i /data/genbank_files/ \
    -o /results/rbcL_extracted.tsv
```

### PBS/Torque Example

```bash
#!/bin/bash
#PBS -N gb_extract
#PBS -l walltime=04:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
module load python/3.9

python ncbi_gb_extract.py \
    -v \
    -g COI \
    ncbi_sequences.gb \
    COI_data.tsv
```

### Snakemake Integration

```python
rule extract_gene:
    input:
        "data/raw/{dataset}.gb"
    output:
        "data/processed/{dataset}_{gene}.tsv"
    params:
        gene="{gene}"
    log:
        "logs/extract_{dataset}_{gene}.log"
    shell:
        """
        python scripts/ncbi_gb_extract.py \
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
               dataset=["ncbi_plants", "ncbi_animals"],
               gene=["rbcL", "matK", "COI"])
```

## Performance Notes

- Memory usage scales with the number of processed records (all held in memory for consistent column ordering)
- For very large datasets (millions of records), consider splitting input files
- Progress output every 1000 records when using `-v`
- Typical throughput: ~10,000-50,000 records/minute depending on record complexity

## Error Handling

- Invalid GenBank files are logged and skipped
- Missing qualifiers result in empty cells (not errors)
- Malformed locations fall back to string representation

## Examples

### Extract rbcL from NCBI download

```bash
# Download sequences from NCBI
esearch -db nucleotide -query "rbcL[gene] AND plants[organism]" | \
    efetch -format gb > plant_rbcl.gb

# Extract to TSV
python ncbi_gb_extract.py -v -g rbcL plant_rbcl.gb plant_rbcl.tsv
```

### Process multiple markers

```bash
for gene in rbcL matK trnL ITS; do
    python ncbi_gb_extract.py -v -g $gene \
        sequences.gb \
        ${gene}_extracted.tsv
done
```

### Filter output columns

```bash
# Extract only accession, organism, and sequence
python ncbi_gb_extract.py -g rbcL input.gb - | \
    csvtk cut -t -f accession,organism,nucleotide_sequence > minimal.tsv
```

## License

MIT

## Author

Ben Price, Natural History Museum London
