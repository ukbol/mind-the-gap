# MIDORI FASTA to TSV Processor

Converts MIDORI reference library FASTA files into structured TSV format for downstream bioinformatics analysis.

## Overview

This script parses MIDORI-formatted FASTA files containing taxonomic annotations and extracts:
- Accession numbers
- Taxonomic hierarchy (kingdom through species)
- NCBI taxonomy IDs (taxid)
- DNA sequences

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Installation

No installation required. Simply copy `process_midori.py` to your working directory or add to your PATH.

```bash
chmod +x process_midori.py  # Make executable on Unix systems
```

## Usage

### Basic Usage

```bash
python process_midori.py input.fasta output.tsv
```

### With Progress Output

```bash
python process_midori.py -v MIDORI2_LONGEST_NUC_GB264_lrRNA_QIIME.fasta lrRNA_processed.tsv
```

### Pipeline Integration (stdin/stdout)

```bash
# From gzipped file
zcat input.fasta.gz | python process_midori.py - - > output.tsv

# Pipe to another tool
python process_midori.py input.fasta - | cut -f1,7,8 > species_list.tsv
```

## Output Format

Tab-separated file with the following columns:

| Column    | Description                                           | Example                    |
|-----------|-------------------------------------------------------|----------------------------|
| accession | GenBank/EMBL accession number                         | KF807066                   |
| kingdom   | Taxonomic kingdom                                     | Metazoa                    |
| phylum    | Taxonomic phylum                                      | Chordata                   |
| class     | Taxonomic class                                       | Amphibia                   |
| order     | Taxonomic order                                       | Anura                      |
| family    | Taxonomic family                                      | Dendrobatidae              |
| genus     | Taxonomic genus                                       | Silverstoneia              |
| species   | Species epithet (binomial)                            | Silverstoneia nubicola     |
| taxid     | NCBI taxonomy ID of lowest identified rank            | 384849                     |
| sequence  | DNA sequence                                          | CGCCTCTTGTT...             |

## Input Format

Expects MIDORI FASTA format with headers structured as:

```
>ACCESSION.version.start.end rank_name_taxid;rank_name_taxid;...;species_Genus species_taxid
SEQUENCE
```

Example:
```
>KF807066.1.<1.>569 root_1;...;kingdom_Metazoa_33208;...;species_Silverstoneia nubicola_384849
CGCCTCTTGTTTATCAATAAGAGGTCCAGCCTGCCC...
```

## HPC Cluster Usage

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=midori_process
#SBATCH --output=midori_%j.out
#SBATCH --error=midori_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load python/3.9  # Adjust for your cluster

python process_midori.py \
    -v \
    MIDORI2_LONGEST_NUC_GB264_lrRNA_QIIME.fasta \
    lrRNA_processed.tsv
```

### PBS/Torque Example

```bash
#!/bin/bash
#PBS -N midori_process
#PBS -l walltime=02:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
module load python/3.9

python process_midori.py \
    -v \
    MIDORI2_LONGEST_NUC_GB264_lrRNA_QIIME.fasta \
    lrRNA_processed.tsv
```

### Snakemake Integration

```python
rule process_midori:
    input:
        "data/raw/{marker}_QIIME.fasta"
    output:
        "data/processed/{marker}_taxonomy.tsv"
    log:
        "logs/process_midori_{marker}.log"
    shell:
        """
        python scripts/process_midori.py -v {input} {output} 2> {log}
        """
```

## Performance

- Memory efficient: processes sequences one at a time (streaming)
- Typical throughput: ~50,000 sequences/minute on standard hardware
- Suitable for files with millions of sequences

## Notes

- Empty taxonomic ranks are represented as empty strings in the output
- The taxid corresponds to the lowest taxonomic level present in the header
- Sequences with missing or malformed headers will have empty taxonomy fields
- Works with both Unix (LF) and Windows (CRLF) line endings

## License

MIT

## Author

Ben Price, Natural History Museum London
