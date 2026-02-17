# ENA Mitogenome Status Processing

Matches ENA mitogenome accession data against a species list to produce a per-species count of available mitogenome accessions. Designed to complement the gap analysis pipeline by identifying which taxa already have complete mitochondrial genomes in ENA.

## Overview

This script:

1. **Loads a target species list** — Taxa defined by valid names and their synonyms (the same format used by `gap_analysis.py`)
2. **Loads ENA mitogenome metadata** — A TSV export of complete mitochondrial genomes from ENA
3. **Matches taxa to mitogenomes** — Using both valid names and synonyms to bridge nomenclatural differences
4. **Outputs a status file** — Preserving all input columns and appending mitogenome-specific columns

## Requirements

- Python 3.7+
- No external dependencies (uses only standard library)

## Data Acquisition

### Step 1: Fetch mitogenomes from ENA

Download complete mitochondrial genome records from ENA using their search API. Note the limit of 200,000 records:

```bash
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=organelle%3D%22mitochondrion%22%20AND%20description%3D%22complete%20genome%22&fields=scientific_name,tax_id,accession&limit=200000&format=tsv" > mito_species_ena.tsv
```

This produces a TSV file with three columns: `accession`, `scientific_name`, and `tax_id`.

### Step 2: Run the mitogenome status analysis

```bash
python mito_status.py \
    --species-list uksi_species_export.tsv \
    --mito-metadata mito_species_ena.tsv \
    --output mito_gap_analysis.tsv
```

## Usage

### Command Line Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--species-list` | Yes | — | Path to species list TSV (with `taxon_name`/`species` + `synonyms` columns) |
| `--mito-metadata` | Yes | — | Path to ENA mitogenome TSV (with `accession`, `scientific_name`, `tax_id` columns) |
| `--output` | Yes | — | Output path for mitogenome status TSV |
| `--log-level` | No | INFO | Logging level: DEBUG, INFO, WARNING, ERROR |

## Input Files

### Species List Format

Tab-separated file with headers. Uses the same format as `gap_analysis.py`:

| Column | Required | Description |
|--------|----------|-------------|
| `taxon_name` | Yes* | Valid species name (preferred over `species`) |
| `species` | Yes* | Valid species name (used if `taxon_name` not present) |
| `synonyms` | No | Semicolon-separated synonyms (can be empty) |
| *(other columns)* | No | Preserved in output |

\* One of `taxon_name` or `species` must be present.

**Example:**
```
taxon_name	synonyms	kingdom	phylum	family
Gammarus pulex	Gammarus fossarum;Gammarus caparti	Animalia	Arthropoda	Gammaridae
Baetis rhodani		Animalia	Arthropoda	Baetidae
```

### ENA Mitogenome Metadata Format

Tab-separated file as downloaded from the ENA Portal API (see Data Acquisition above):

| Column | Required | Description |
|--------|----------|-------------|
| `accession` | Yes | ENA sequence accession |
| `scientific_name` | Yes | Species name as recorded in ENA |
| `tax_id` | No | NCBI taxonomy ID |

## Output

### Mitogenome Status TSV

The main output file contains all input species-list columns plus:

| Column | Description |
|--------|-------------|
| `mitogenome_count` | Number of mitogenome accessions matching the taxon |
| `species_status` | GREEN (mitogenome available) or BLACK (no mitogenome found) |
| `mitogenome_accessions` | Semicolon-separated list of matching ENA accessions |

### Unmatched Records TSV

An additional file (`<output_stem>_unmatched.tsv`) is automatically written alongside the main output. It contains mitogenome records from ENA that did not match any taxon in the species list, useful for manual review and identifying potential gaps in the species list itself.

| Column | Description |
|--------|-------------|
| `accession` | ENA sequence accession |
| `scientific_name` | Species name as recorded in ENA |
| `tax_id` | NCBI taxonomy ID |

## Status Colour Mapping

| Status | Meaning |
|--------|---------|
| **GREEN** | Mitogenome(s) available in ENA |
| **BLACK** | No mitogenome found in ENA |

## Matching Logic

The script matches ENA records to species-list taxa by normalizing names (lowercasing, replacing underscores with spaces) and checking against both valid names and synonyms. This ensures that nomenclatural differences between ENA and the species list are bridged automatically.

For each taxon, all matching accessions across all names (valid + synonyms) are aggregated into the result.

## Example

```bash
# Download mitogenome data from ENA
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=organelle%3D%22mitochondrion%22%20AND%20description%3D%22complete%20genome%22&fields=scientific_name,tax_id,accession&limit=200000&format=tsv" > mito_species_ena.tsv

# Run mitogenome status analysis
python mito_status.py \
    --species-list uksi_species_export.tsv \
    --mito-metadata mito_species_ena.tsv \
    --output 2026-02-17_ena_mitogenome_gap_analysis.tsv
```

### Example Output

```
taxon_name        synonyms          kingdom   mitogenome_count  species_status  mitogenome_accessions
Gammarus pulex    Gammarus foss...  Animalia  3                 GREEN           AB123456;CD789012;EF345678
Baetis rhodani                      Animalia  0                 BLACK
```

## Summary Statistics

The script logs a summary after processing, including:

- Total taxa and how many have/lack mitogenomes
- Distribution of accession counts (min, median, max)
- Bracket breakdown (1, 2-5, 6-10, 11-50, 51+ accessions)

## Troubleshooting

**"Species list must have a 'taxon_name' or 'species' column"**
- Ensure your input file is tab-separated
- Check that the header row contains `taxon_name` or `species`

**"Mitogenome file missing required columns"**
- The ENA download must contain `accession` and `scientific_name` columns
- Re-download using the curl command in Data Acquisition above

**Encoding errors**
- The script automatically tries UTF-8 then Latin-1
- If issues persist, convert files to UTF-8

**ENA API returns fewer records than expected**
- The API has a limit parameter (set to 200,000 in the example)
- If the total number of matching records exceeds this limit, increase the `limit` parameter in the URL

## Author

Ben Price / Claude
Date: 2025-02-17
