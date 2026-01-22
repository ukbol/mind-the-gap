# UNITE FASTA to TSV Processor

Converts UNITE reference library FASTA files into structured TSV format for downstream bioinformatics analysis.

## Overview

This script parses UNITE-formatted FASTA files containing fungal/eukaryotic ITS taxonomic annotations and extracts:
- Accession numbers (UDB identifiers)
- Species Hypothesis (SH) cluster codes
- Taxonomic hierarchy (kingdom through species)
- DNA sequences

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Installation

No installation required. Simply copy `process_unite.py` to your working directory or add to your PATH.

```bash
chmod +x process_unite.py  # Make executable on Unix systems
```

## Usage

### Basic Usage

```bash
python process_unite.py input.fasta output.tsv
```

### With Progress Output

```bash
python process_unite.py -v sh_general_release_all_10.0.fasta unite_processed.tsv
```

### Pipeline Integration (stdin/stdout)

```bash
# From gzipped file
zcat input.fasta.gz | python process_unite.py - - > output.tsv

# Pipe to another tool
python process_unite.py input.fasta - | cut -f1,2,8,9 > species_clusters.tsv
```

## Output Format

Tab-separated file with the following columns:

| Column    | Description                                           | Example                              |
|-----------|-------------------------------------------------------|--------------------------------------|
| accession | UNITE database accession (UDB ID)                     | UDB02393310                          |
| cluster   | Species Hypothesis (SH) cluster code                  | SH0000001.10FU                       |
| kingdom   | Taxonomic kingdom                                     | Alveolata                            |
| phylum    | Taxonomic phylum                                      | Ciliophora                           |
| class     | Taxonomic class                                       | Spirotrichea                         |
| order     | Taxonomic order                                       | Sporadotrichida                      |
| family    | Taxonomic family                                      | Oxytrichidae                         |
| genus     | Taxonomic genus                                       | Oxytrichidae_gen_Incertae_sedis      |
| species   | Species epithet                                       | Oxytrichidae_sp                      |
| sequence  | DNA sequence (ITS region)                             | ACACTAATCCAACTCAACT...               |

## Input Format

Expects UNITE FASTA format with headers structured as:

```
>Name|ACCESSION|CLUSTER|type|k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
SEQUENCE
```

Example:
```
>Ciliophora_sp|UDB02393310|SH0000001.10FU|reps_singleton|k__Alveolata;p__Ciliophora;c__Ciliophora_cls_Incertae_sedis;o__Ciliophora_ord_Incertae_sedis;f__Ciliophora_fam_Incertae_sedis;g__Ciliophora_gen_Incertae_sedis;s__Ciliophora_sp
ACACTAATCCAACTCAACTCAACGAAGCCTTCAGTTGCAGTAGCAGTACTAGCGAAAGCCTGTGCCGCCGCTGCAGCACCAAAAAACTAATTCAAACAAAGGAGCCTAACT...
```

## HPC Cluster Usage

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=unite_process
#SBATCH --output=unite_%j.out
#SBATCH --error=unite_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load python/3.9  # Adjust for your cluster

python process_unite.py \
    -v \
    sh_general_release_all_10.0_29.04.2024.fasta \
    unite_processed.tsv
```

### PBS/Torque Example

```bash
#!/bin/bash
#PBS -N unite_process
#PBS -l walltime=02:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
module load python/3.9

python process_unite.py \
    -v \
    sh_general_release_all_10.0_29.04.2024.fasta \
    unite_processed.tsv
```

### Snakemake Integration

```python
rule process_unite:
    input:
        "data/raw/sh_general_release_{version}.fasta"
    output:
        "data/processed/unite_{version}_taxonomy.tsv"
    log:
        "logs/process_unite_{version}.log"
    shell:
        """
        python scripts/process_unite.py -v {input} {output} 2> {log}
        """
```

## Performance

- Memory efficient: processes sequences one at a time (streaming)
- Typical throughput: ~50,000 sequences/minute on standard hardware
- Suitable for files with millions of sequences

## Notes

- Empty taxonomic ranks are represented as empty strings in the output
- Incertae sedis annotations are preserved as-is (e.g., `Ciliophora_cls_Incertae_sedis`)
- Works with both Unix (LF) and Windows (CRLF) line endings
- Compatible with all UNITE release formats (general, dynamic, etc.)

## UNITE Database

UNITE is a database for molecular identification of fungi and other eukaryotes. The Species Hypothesis (SH) system provides stable, named clusters for environmental sequence identification.

For more information: https://unite.ut.ee/

## License

MIT

## Author

Ben Price, Natural History Museum London
