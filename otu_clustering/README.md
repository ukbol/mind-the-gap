# OTU Clustering

Clusters DNA sequences into Operational Taxonomic Units (OTUs) using VSEARCH. Accepts TSV input with sequence data and outputs an annotated version of the input with an additional OTU_ID column.

## Overview

This script performs sequence clustering for taxonomic analysis and biodiversity assessment. It:
- Reads sequences from TSV files (common output from database exports)
- Handles sequences in either orientation (forward or reverse complement)
- Uses VSEARCH for fast, accurate clustering
- Outputs the original TSV with an OTU_ID column appended
- Preserves all original data - rows without valid sequences have an empty OTU_ID

## Requirements

- Python 3.6+
- VSEARCH (tested with v2.21+)

### Installing VSEARCH

```bash
# Conda (recommended)
conda install -c bioconda vsearch

# Ubuntu/Debian
sudo apt-get install vsearch

# macOS with Homebrew
brew install vsearch

# From source
# See https://github.com/torognes/vsearch
```

## Installation

```bash
# No additional Python dependencies required (uses standard library only)

# Make executable on Unix systems
chmod +x otu_clustering.py
```

## Usage

### Basic Usage

```bash
# Cluster at 99% identity
python otu_clustering.py -t 0.99 input.tsv output.tsv

# Cluster at 97% identity (common for species-level OTUs)
python otu_clustering.py -t 0.97 sequences.tsv otus.tsv
```

### With Custom Column Names

```bash
# Specify accession and sequence column names
python otu_clustering.py -t 0.99 \
    --accession-col recordid \
    --sequence-col nuc \
    input.tsv output.tsv
```

### With Verbose Output

```bash
python otu_clustering.py -v -t 0.99 input.tsv output.tsv
```

### Pipeline Integration (stdin/stdout)

```bash
# From gzipped file
zcat sequences.tsv.gz | python otu_clustering.py -t 0.99 - - > otus.tsv

# Chain with other tools
cat data.tsv | python otu_clustering.py -t 0.99 - - | grep "OTU_"
```

## Command Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `input` | Yes | Input TSV file (use "-" for stdin) |
| `output` | Yes | Output TSV file (use "-" for stdout) |
| `-t, --threshold` | Yes | Similarity threshold (0.0-1.0, e.g., 0.99 for 99%) |
| `--threads` | No | Number of threads for VSEARCH (default: 8) |
| `--strand` | No | Strand orientation: "plus" or "both" (default: both) |
| `--min-length` | No | Minimum sequence length to include (default: 100) |
| `--accession-col` | No | Name of accession column (default: accession) |
| `--sequence-col` | No | Name of sequence column (can specify multiple) |
| `--temp-dir` | No | Temporary directory for intermediate files |
| `-v, --verbose` | No | Print progress and statistics to stderr |

## Input Format

Tab-separated file with at minimum:
- An accession/identifier column
- A sequence column

### Default Column Detection

The script tries these column names in order:
- Accession: `accession`
- Sequence: `sequence`, `nucleotide_sequence`, `nuc`

### Example Input

```tsv
accession	organism	sequence
ABC123	Homo sapiens	ATCGATCGATCGATCG...
DEF456	Mus musculus	GCTAGCTAGCTAGCTA...
GHI789	Rattus norvegicus	ATCGATCGATCGATCG...
JKL012	Canis familiaris	
```

## Output Format

The output is the input TSV with an additional `OTU_ID` column appended at the end. All original columns and data are preserved.

- Rows with valid sequences receive an OTU identifier (e.g., `OTU_000001`)
- Rows without valid sequences (empty, too short, invalid characters) have an empty OTU_ID field

### Example Output

```tsv
accession	organism	sequence	OTU_ID
ABC123	Homo sapiens	ATCGATCGATCGATCG...	OTU_000001
DEF456	Mus musculus	GCTAGCTAGCTAGCTA...	OTU_000002
GHI789	Rattus norvegicus	ATCGATCGATCGATCG...	OTU_000001
JKL012	Canis familiaris		
```

Note: ABC123 and GHI789 share the same OTU because their sequences are identical (or similar above threshold).

## Clustering Parameters

### Threshold Selection

Common threshold values:
- **0.99 (99%)**: Very stringent, often used for intraspecific variation
- **0.97 (97%)**: Traditional bacterial species-level OTUs
- **0.95 (95%)**: Genus-level clustering
- **0.90 (90%)**: Family-level clustering

### Strand Handling

- `--strand plus`: Only compare sequences in forward orientation
- `--strand both`: Compare both forward and reverse complement (default)

Use `both` when sequences may be in different orientations (common with database extracts).

## HPC Cluster Usage

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=otu_cluster
#SBATCH --output=otu_cluster_%j.out
#SBATCH --error=otu_cluster_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

module load vsearch
module load python/3.9

python otu_clustering.py \
    -v \
    -t 0.99 \
    --threads 16 \
    /data/sequences.tsv \
    /results/sequences_with_otus.tsv
```

### PBS/Torque Example

```bash
#!/bin/bash
#PBS -N otu_cluster
#PBS -l walltime=04:00:00
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=16

cd $PBS_O_WORKDIR
module load vsearch
module load python/3.9

python otu_clustering.py \
    -v \
    -t 0.97 \
    --threads 16 \
    sequences.tsv \
    sequences_with_otus.tsv
```

### Snakemake Integration

```python
rule cluster_otus:
    input:
        "data/processed/{dataset}.tsv"
    output:
        "results/otus/{dataset}_t{threshold}.tsv"
    params:
        threshold=lambda w: f"0.{w.threshold}"
    threads: 16
    log:
        "logs/otu_clustering/{dataset}_t{threshold}.log"
    shell:
        """
        python scripts/otu_clustering.py \
            -v \
            -t {params.threshold} \
            --threads {threads} \
            {input} \
            {output} \
            2> {log}
        """

# Example usage in workflow
rule all:
    input:
        expand("results/otus/{dataset}_t{threshold}.tsv",
               dataset=["bold_coi", "ncbi_rbcl"],
               threshold=["99", "97", "95"])
```

### Snakemake with Conda Environment

```yaml
# envs/vsearch.yaml
name: vsearch
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - vsearch>=2.21
  - python>=3.6
```

```python
rule cluster_otus:
    input:
        "data/processed/{dataset}.tsv"
    output:
        "results/otus/{dataset}_otus.tsv"
    params:
        threshold=0.99
    threads: 16
    conda:
        "envs/vsearch.yaml"
    log:
        "logs/otu_clustering/{dataset}.log"
    shell:
        """
        python scripts/otu_clustering.py \
            -v \
            -t {params.threshold} \
            --threads {threads} \
            {input} \
            {output} \
            2> {log}
        """
```

## Performance Notes

- VSEARCH is highly efficient and can cluster millions of sequences
- Memory usage scales with the number and length of sequences
- Threading provides near-linear speedup for large datasets
- Typical throughput: 100,000+ sequences/minute depending on length and threshold

### Optimizing Performance

```bash
# For very large datasets, increase threads
python otu_clustering.py -t 0.99 --threads 32 large_dataset.tsv output.tsv

# Use local SSD for temp files on HPC
python otu_clustering.py -t 0.99 --temp-dir /scratch/$USER/temp input.tsv output.tsv
```

## Algorithm Details

The script uses VSEARCH's `cluster_fast` algorithm:
1. Sequences are sorted by abundance (or length)
2. The first sequence becomes a centroid (OTU seed)
3. Subsequent sequences are compared to existing centroids
4. If similarity >= threshold, sequence joins that OTU
5. Otherwise, sequence becomes a new centroid

This greedy algorithm is fast and produces reproducible results.

## Verbose Output

With `-v`, statistics are printed to stderr:

```
[2025-01-22 10:30:00] INFO: Starting OTU clustering
[2025-01-22 10:30:00] INFO: Input: sequences.tsv
[2025-01-22 10:30:00] INFO: Threshold: 0.99
[2025-01-22 10:30:01] INFO: Total rows: 150500
[2025-01-22 10:30:01] INFO: Valid sequences for clustering: 150000
[2025-01-22 10:30:01] INFO: Rows without valid sequence: 500
[2025-01-22 10:30:15] INFO: VSEARCH clustering completed
[2025-01-22 10:30:15] INFO: Summary statistics:
[2025-01-22 10:30:15] INFO:   Total rows in input: 150500
[2025-01-22 10:30:15] INFO:   Sequences clustered: 150000
[2025-01-22 10:30:15] INFO:   Rows without OTU assignment: 500
[2025-01-22 10:30:15] INFO:   Total OTUs created: 12500
[2025-01-22 10:30:15] INFO:   Singleton OTUs: 3200
[2025-01-22 10:30:15] INFO:   Largest OTU: OTU_000001 (450 sequences)
[2025-01-22 10:30:15] INFO:   Mean OTU size: 12.00
[2025-01-22 10:30:15] INFO:   Median OTU size: 4.0
```

## Error Handling

- Missing VSEARCH binary: Clear error message with installation instructions
- Invalid sequences: Preserved in output with empty OTU_ID
- Missing columns: Lists available columns for debugging
- Empty input: Produces valid output file with header only

## Related Scripts

This script is designed to work with other tools in the mind-the-gap pipeline:
- `ncbi_gb_extract.py`: Extract sequences from GenBank files
- `bold_gene_extract.py`: Extract sequences from BOLD exports

## License

MIT

## Author

Ben Price, Natural History Museum London
