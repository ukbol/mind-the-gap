# GenBank Batch Downloader

Batch downloads GenBank flat files from NCBI using EDirect, with automatic batching and rate-limit compliance. Configured by default to fetch rbcL (ribulose-1,5-bisphosphate carboxylase) records, but can be adapted for any gene or query.

## Overview

This shell script uses NCBI's EDirect command-line tools to:
- Query NCBI Nucleotide for records matching a specified gene
- Count total matching records
- Download in configurable batches (default 5,000 records per file)
- Skip already-downloaded batches for resumability
- Retry failed downloads after a 30-second delay
- Respect NCBI rate limits (3-second delay between batches)

## Requirements

- Bash
- NCBI EDirect tools (`esearch`, `efetch`, `xtract`)

### Installing EDirect

```bash
# Conda (recommended)
conda install -c bioconda entrez-direct

# Manual installation
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

## Usage

The script is designed to run as a SLURM batch job:

```bash
sbatch edirect_genbank_fetch.sh
```

Or run directly:

```bash
bash edirect_genbank_fetch.sh
```

## Configuration

Edit the following variables in the script to customise behaviour:

| Variable | Default | Description |
|----------|---------|-------------|
| `OUTDIR` | `/mnt/shared/scratch/museomix/rbcl_genbank` | Output directory for downloaded files |
| `QUERY` | `(rbcl[Gene Name]) OR (rbc-l[Gene Name]) OR (rubisco[Gene Name])` | NCBI search query |
| `BATCH` | `5000` | Records per download batch |

### Customising the Query

To download a different gene, modify the `QUERY` variable:

```bash
# COI / COX1 (animal barcode)
QUERY='(COI[Gene Name]) OR (COX1[Gene Name]) OR (cytochrome oxidase subunit I[Gene Name])'

# matK (plant barcode)
QUERY='(matK[Gene Name]) OR (maturase K[Gene Name])'

# ITS (fungal barcode)
QUERY='(ITS[Gene Name]) AND fungi[Organism]'
```

## Output

Downloaded files are named sequentially:
```
rbcl_batch_000000.gb
rbcl_batch_005000.gb
rbcl_batch_010000.gb
...
```

Each file contains up to `BATCH` GenBank flat-file records, suitable for parsing with `ncbi_gb_extract.py`.

## SLURM Configuration

The script includes SLURM directives:

```bash
#SBATCH --job-name=16s_gb_download
#SBATCH --output=16s_download_%j.log
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --partition=himem
```

Adjust these for your cluster. Memory and CPU requirements are modest for downloading; the high-memory partition is used for reliability on long-running jobs.

## Rate Limiting

- A 3-second delay between batches respects NCBI's default rate limit of 3 requests/second
- Setting the `NCBI_API_KEY` environment variable increases the limit to 10 requests/second:

```bash
export NCBI_API_KEY="your_api_key_here"
```

Register for an API key at: https://www.ncbi.nlm.nih.gov/account/settings/

## Resumability

The script checks for existing non-empty output files and skips them, making it safe to restart after interruptions without re-downloading completed batches.

## Pipeline Integration

Downloaded GenBank files are parsed with the companion script `ncbi_gb_extract.py`:

```bash
# Download
sbatch edirect_genbank_fetch.sh

# Extract to TSV
python ../ncbi_gb_extract/ncbi_gb_extract.py \
    -v -g rbcL \
    -i /path/to/rbcl_genbank/ \
    -o rbcl_extracted.tsv
```

## License

MIT

## Author

Ben Price, Natural History Museum London
