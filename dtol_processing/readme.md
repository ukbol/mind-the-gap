# DToL Status Processor

Matches Darwin Tree of Life (DToL) genome sequencing metadata against a UKSI species list to determine how far each UK species has progressed through the DToL sequencing pipeline.

## Overview

This script cross-references the [Darwin Tree of Life portal](https://portal.darwintreeoflife.org/) metadata export with a UKSI species list (the same format used by `gap_analysis.py`) to produce a per-species status panel. Matching uses both valid names and synonyms so that nomenclatural differences between DToL and UKSI are bridged automatically.

The output is a TSV file that preserves all input species-list columns and appends DToL-specific columns, following the same pattern as the gap analysis output.

## Requirements

- Python 3.7+
- No external dependencies (uses only standard library)

## Usage

### Basic Usage

```bash
python dtol_status.py \
    --species-list uksi_species_export.tsv \
    --dtol-metadata dtol_download.csv \
    --output dtol_gap_analysis.tsv
```

### Command Line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--species-list` | Yes | -- | Path to UKSI species list TSV (with `taxon_name`/`species` + `synonyms` columns) |
| `--dtol-metadata` | Yes | -- | Path to DToL portal CSV export |
| `--output` | Yes | -- | Output path for DToL status TSV |
| `--log-level` | No | INFO | Logging level: DEBUG, INFO, WARNING, ERROR |

## Input Files

### Species List

Tab-separated file in the same format used by `gap_analysis.py`:

| Column | Required | Description |
|--------|----------|-------------|
| `taxon_name` or `species` | Yes | Valid species name |
| `synonyms` | No | Semicolon-separated synonyms |
| *(other columns)* | No | Preserved in output |

### DToL Metadata

CSV export from the [DToL Data Portal](https://portal.darwintreeoflife.org/):

| Column | Required | Description |
|--------|----------|-------------|
| `Organism` | Yes | Scientific name |
| `Current Status` | Yes | Pipeline stage |
| `Common Name` | No | Vernacular name |
| `INSDC ID` | No | INSDC accession |
| `ToL ID` | No | Tree of Life identifier(s), comma-separated |

## Output Files

### Main Output (`--output`)

The input species list with three columns appended:

| Column | Description |
|--------|-------------|
| `dtol_status` | DToL pipeline stage (e.g., `Annotation Complete`) or `Not in DToL` |
| `species_status` | Traffic-light colour: GREEN, BLUE, AMBER, RED, or BLACK |
| `dtol_insdc_ids` | Semicolon-separated INSDC accessions from all matching records |

### Unmatched File (`{output_stem}_unmatched.tsv`)

DToL organisms that did not match any UKSI taxon, written alongside the main output for manual review.

| Column | Description |
|--------|-------------|
| `organism` | DToL organism name |
| `common_name` | Common name |
| `current_status` | Pipeline stage |
| `insdc_id` | INSDC accession |
| `tol_ids` | Semicolon-separated ToL IDs |

## DToL Pipeline Stages

Stages are ordered from furthest along to earliest entry:

| Stage | Status Colour | Meaning |
|-------|---------------|---------|
| Annotation Complete | GREEN | Genome annotated and published |
| Assemblies - Submitted | BLUE | Genome assembly submitted |
| Raw Data - Submitted | AMBER | Raw sequencing data submitted |
| Submitted to BioSamples | RED | Sample registered, sequencing pending |
| Not in DToL | BLACK | Species not in the DToL pipeline |

When multiple DToL records match a single UKSI taxon (e.g., via valid name and a synonym), the record at the furthest pipeline stage is used for the status, while INSDC IDs and ToL IDs are aggregated from all matches.

## Matching Logic

1. Normalize all names (lowercase, underscores to spaces)
2. Check the valid species name against DToL organism names
3. Check each synonym against DToL organism names
4. If multiple DToL records match, select the one at the furthest pipeline stage
5. INSDC IDs and ToL IDs are collected from all matching records

## Example

```bash
# Download metadata from DToL portal as CSV
# https://portal.darwintreeoflife.org/data

# Run status matching
python dtol_status.py \
    --species-list ../uksi_processing/uksi_db/uksi_species_export.tsv \
    --dtol-metadata dtol_portal_download.csv \
    --output dtol_gap_analysis.tsv

# Check the summary statistics in the log output
```

Example summary output:
```
DToL STATUS SUMMARY
============================================================
Total UKSI taxa: 76,640
By pipeline stage:
  Annotation Complete                1,234  (  1.6%)
  Assemblies - Submitted               567  (  0.7%)
  Raw Data - Submitted                 890  (  1.2%)
  Submitted to BioSamples              456  (  0.6%)
  Not in DToL                       73,493  ( 95.9%)
```

## Data Source

DToL metadata is downloaded from the Darwin Tree of Life Data Portal:
https://portal.darwintreeoflife.org/

## License

MIT

## Author

Ben Price, Natural History Museum London
