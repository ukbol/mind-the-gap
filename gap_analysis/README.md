# HPC Gap Analysis for DNA Barcode Library Curation

A high-performance Python script for taxon-centric gap analysis of DNA barcode reference libraries. Designed for HPC environments with support for parallel processing of large datasets (millions of records, hundreds of thousands of taxa).

## Overview

This tool analyzes the quality of DNA barcode reference data by:

1. **Loading a target species list** — Taxa defined by valid names and their synonyms
2. **Scanning records** — Finding all records matching each taxon's names
3. **Analyzing BIN/OTU sharing** — Detecting taxonomic conflicts where different taxa share clusters
4. **Assigning grades and status** — BAGS grades (A-F) and traffic light status (GREEN/AMBER/BLUE/ORANGE/RED/BLACK)
5. **Writing filtered records** — Exporting a subset of the records file containing only records in clusters relevant to the target taxa

## Key Concepts

### Taxon Definition

Each row in the input species list defines a **taxon** comprising:
- **Valid name**: The accepted species name
- **Synonyms**: Alternative names (junior synonyms, misspellings, etc.)

The taxon's "name set" = valid name + all synonyms. Analysis is performed at the taxon level, aggregating all records matching any name in the set.

### Cluster Identity

The script uses either `bin_uri` (preferred) or `otu_id` as the cluster identifier:
- **File-level decision**: If `bin_uri` column exists, it's used for all records; otherwise `otu_id`
- **Multiple clusters**: Records with pipe-separated values (e.g., `BOLD:AAA1234|BOLD:AAA5678`) belong to multiple clusters

### Name Normalisation

Species names are normalised before matching: lowercased and underscores replaced with spaces. This handles differences between naming conventions (e.g. BOLD/UKSI use spaces, UNITE uses underscores).

### Linnaean vs Placeholder Names

When checking for BIN/OTU conflicts, the script distinguishes between:
- **Linnaean binomials** — formally described species (e.g. *Genus species*)
- **Placeholder/provisional names** — names containing `cf.`, `aff.`, `sp.`, digits, or other non-Linnaean markers

This distinction drives the RED vs ORANGE status assignment (see Grading System).

## Installation

### Requirements

- Python 3.7+
- No external dependencies (uses only standard library)

```bash
# Clone the repository
git clone https://github.com/your-repo/mind-the-gap.git
cd mind-the-gap/gap_analysis
```

## Usage

### Basic Usage

```bash
python gap_analysis.py \
    --species-list species.tsv \
    --records result_output.tsv \
    --output gap_analysis.tsv
```

### HPC Usage with Parallel Processing

```bash
python gap_analysis.py \
    --species-list species.tsv \
    --records result_output.tsv \
    --output gap_analysis.tsv \
    --workers 32 \
    --batch-size 2000
```

### Raw BOLD Data (BOLD Mode)

For raw BOLD TSV downloads, use `--bold` to enable BOLD-specific parsing and default filters (Animalia kingdom, COI-5P marker, valid species names only):

```bash
python gap_analysis.py \
    --species-list species.tsv \
    --records BOLD_Public.tsv \
    --output gap_analysis.tsv \
    --bold
```

With custom BOLD filters:

```bash
python gap_analysis.py \
    --species-list species.tsv \
    --records BOLD_Public.tsv \
    --output gap_analysis.tsv \
    --bold --marker matK --kingdom-list Plantae Fungi --no-filter-species
```

### Command Line Options

**Core options:**

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--species-list` | Yes | — | Path to species list TSV file |
| `--records` | Yes | — | Path to records TSV file (e.g., result_output.tsv) |
| `--output` | Yes | — | Output path for gap analysis TSV |
| `--filtered-records` | No | `<output_stem>_filtered_records.tsv` | Output path for filtered records TSV (auto-derived if omitted) |
| `--workers` | No | CPU count | Number of parallel workers |
| `--batch-size` | No | 1000 | Taxa per batch for parallel processing |
| `--log-level` | No | INFO | Logging level: DEBUG, INFO, WARNING, ERROR |

**BOLD mode options** (only meaningful when `--bold` is set):

| Option | Default in BOLD mode | Description |
|--------|----------------------|-------------|
| `--bold` | — | Enable BOLD mode: QUOTE_NONE parsing, field sanitisation, and default filters |
| `--filter-species` / `--no-filter-species` | ON | Include only records with valid species names (skips empty/`None` values) |
| `--filter-kingdom` / `--no-filter-kingdom` | ON | Filter records by kingdom |
| `--kingdom-list <K1> <K2> ...` | `Animalia` | Kingdom(s) to keep when kingdom filter is active |
| `--marker <CODE>` | `COI-5P` | Filter by `marker_code` column (set to empty string to disable) |

## Input Files

### Species List Format

Tab-separated file with headers:

| Column | Required | Description |
|--------|----------|-------------|
| `taxon_name` | Yes (preferred) | Valid species name — takes priority over `species` if both columns are present |
| `species` | Yes (fallback) | Valid species name — used if `taxon_name` column is absent |
| `synonyms` | No | Semicolon-separated synonyms (can be empty) |
| *(other columns)* | No | Preserved verbatim in output |

**Example:**
```
species	synonyms	kingdom	phylum	family
Gammarus pulex	Gammarus fossarum;Gammarus caparti	Animalia	Arthropoda	Gammaridae
Baetis rhodani		Animalia	Arthropoda	Baetidae
Glossiphonia paludosa	Bactracobdella paludosa;Batracobdella paludosa	Animalia	Annelida	Glossiphoniidae
```

### Records File Format

Tab-separated file (e.g., `result_output.tsv`) with:

| Column | Required | Description |
|--------|----------|-------------|
| `species` | Yes (preferred) | Species name — used for matching against species list |
| `organism` | Yes (fallback) | Species name in NCBI-format files — used if `species` column is absent |
| `subspecies` | No | Subspecies name — also matched against the species list if present |
| `bin_uri` or `otu_id` | Yes | Cluster identifier(s), pipe-separated if multiple |
| `country_iso` | No | ISO country code — used to count UK (GB) records |
| `Country` | No | Country name (fallback for UK detection when `country_iso` is absent or empty) |
| `kingdom` | No | Kingdom — used by `--filter-kingdom` in BOLD mode |
| `marker_code` | No | Marker/gene code — used by `--marker` filter in BOLD mode |

**Notes:**
- The cluster column used for analysis is determined at the file level: `bin_uri` is preferred over `otu_id`; both are tracked separately in output columns regardless.
- When `subspecies` is present, records are indexed under both the `species` and `subspecies` values, so a taxon matching a subspecies name will accumulate those records.
- UK record counts use `country_iso == 'GB'` when available, falling back to matching the `Country` column against "United Kingdom"/"UK".

## Output

### Gap Analysis TSV (`--output`)

The output file contains all input columns from the species list, plus:

| Column | Description |
|--------|-------------|
| `number_records` | Total records matching any name in the taxon |
| `gb_records` | Records from the UK (where `country_iso=GB` or `Country` matches "United Kingdom") |
| `bags_grade` | Quality grade: A, B, C, D, E, or F |
| `species_status` | Traffic light status: GREEN, AMBER, BLUE, ORANGE, RED, or BLACK |
| `bin_uris` | Semicolon-separated list of distinct `bin_uri` values found for this taxon |
| `otu_ids` | Semicolon-separated list of distinct `otu_id` values found for this taxon |
| `other_names` | Semicolon-separated names sharing BIN/OTU but not in this taxon's name set |

### Filtered Records TSV (`--filtered-records`)

A second output file is automatically written alongside the gap analysis TSV. It contains a filtered copy of the input records file, keeping only records whose `bin_uri` or `otu_id` appears in a cluster associated with at least one target taxon.

This file is useful for downstream visualisation and manual inspection — it includes both the target species records and any other-name records sharing the same clusters, while omitting unrelated records.

The output path defaults to `<output_stem>_filtered_records.tsv` but can be set explicitly with `--filtered-records`.

## Grading System

### BAGS Grades

| Grade | Meaning | Criteria |
|-------|---------|----------|
| **A** | Excellent | 1 BIN/OTU, ≥11 records, no conflicts |
| **B** | Good | 1 BIN/OTU, 3–10 records, no conflicts |
| **C** | Split | Multiple BIN/OTUs, no conflicts |
| **D** | Low coverage | 1 BIN/OTU, <3 records, no conflicts |
| **E** | Conflict | Other taxa share the BIN/OTU |
| **F** | No data | No records found (or records present but no cluster assignment) |

### Traffic Light Status

| Status | Colour | Meaning | Criteria |
|--------|--------|---------|----------|
| **GREEN** | 🟢 | Clean | Only valid name recorded, no conflicts |
| **AMBER** | 🟡 | Nomenclatural mess | Both valid name AND synonym(s) recorded (needs cleanup) |
| **BLUE** | 🔵 | Valid name absent | Only synonym(s) recorded, valid name missing from database |
| **ORANGE** | 🟠 | Provisional conflict | BIN/OTU shared only with placeholder/provisional names (non-Linnaean: `sp.`, `cf.`, `aff.`, numeric codes, etc.) |
| **RED** | 🔴 | Taxonomic conflict | BIN/OTU shared with at least one formally described species outside this taxon |
| **BLACK** | ⚫ | No coverage | No records for this taxon |

**Note on ORANGE vs RED:** The script tests each external name sharing a BIN/OTU against a set of criteria for a valid Linnaean binomial (two-part name, no placeholder markers, no digits). If all sharing names are provisional/placeholder, the status is ORANGE; if any is a proper species name, the status is RED. Both cases receive BAGS grade **E**.

### Decision Logic

```
┌────────────────────────────────────────────────────────────────┐
│                     GRADING DECISION TREE                      │
├────────────────────────────────────────────────────────────────┤
│                                                                │
│  Has records?                                                  │
│  ├── NO  → Grade F, Status BLACK                              │
│  └── YES                                                       │
│       │                                                        │
│       ▼                                                        │
│  Other names share BIN/OTU?                                    │
│  ├── YES                                                       │
│  │    ├── Any sharer is a Linnaean binomial? → Grade E, RED   │
│  │    └── All sharers are placeholder names? → Grade E, ORANGE│
│  └── NO                                                        │
│       │                                                        │
│       ▼                                                        │
│  Which names are recorded?                                     │
│  ├── Valid + synonym(s) → Status AMBER                        │
│  ├── Valid only         → Status GREEN                        │
│  └── Synonym(s) only   → Status BLUE                         │
│       │                                                        │
│       ▼                                                        │
│  How many BIN/OTUs?                                            │
│  ├── 0 BIN/OTUs (records exist but unassigned) → Grade F      │
│  ├── 1 BIN/OTU                                                │
│  │    ├── ≥11 records → Grade A                               │
│  │    ├── 3–10 records → Grade B                              │
│  │    └── <3 records  → Grade D                               │
│  └── Multiple BIN/OTUs → Grade C                              │
│                                                                │
└────────────────────────────────────────────────────────────────┘
```

## Examples

### Example Scenarios

| Taxon | BIN contains | Status | Grade | Explanation |
|-------|--------------|--------|-------|-------------|
| Valid: "Alpha vulgaris", Syn: "Alpha communis" | Only "Alpha vulgaris" (15 records) | GREEN | A | Clean, good coverage |
| Valid: "Beta marina", Syn: "Beta aquatica" | Both "Beta marina" and "Beta aquatica" (8 records) | AMBER | B | Nomenclatural cleanup needed |
| Valid: "Gamma riparia", Syn: "Gamma fluviatilis" | Only "Gamma fluviatilis" (5 records) | BLUE | B | Valid name not in BOLD |
| Valid: "Delta palustris", Syn: none | "Delta palustris" + "Epsilon montanus" | RED | E | Taxonomic conflict with named species |
| Valid: "Theta lacustris", Syn: none | "Theta lacustris" + "Theta sp. BOLD:AAB1234" | ORANGE | E | Conflict only with placeholder name |
| Valid: "Zeta alpina", Syn: "Zeta montana" | No records | BLACK | F | No coverage |
| Valid: "Eta borealis", Syn: none | "Eta borealis" in 3 different BINs | GREEN | C | Split across clusters |

### Example Output

```
species           synonyms        kingdom   number_records  gb_records  bags_grade  species_status  bin_uris            otu_ids  other_names
Alpha vulgaris    Alpha communis  Animalia  15              12          A           GREEN           BOLD:AAA0001
Beta marina       Beta aquatica   Animalia  8               2           B           AMBER           BOLD:AAB0002
Gamma riparia     Gamma fluv...   Animalia  5               0           B           BLUE            BOLD:AAC0003
Delta palustris                   Animalia  12              8           E           RED             BOLD:AAD0004                 Epsilon montanus
Theta lacustris                   Animalia  3               1           E           ORANGE          BOLD:AAE0005                 Theta sp. BOLD:AAE0005
Zeta alpina       Zeta montana    Animalia  0               0           F           BLACK
```

## BOLD Mode Details

Raw BOLD TSV downloads contain formatting quirks that can break standard CSV parsers:

- Unescaped double-quote characters (e.g. `"Syn.`) in free-text fields
- Embedded newlines and carriage returns within field values

`--bold` mode handles these by:

1. Using `QUOTE_NONE` parsing to avoid misinterpreting quote characters
2. Sanitising each field: stripping embedded newlines/CRs and removing double quotes

In addition, `--bold` mode enables three filters applied **before** indexing (reducing memory usage and noise):

| Filter | Flag | Default | Description |
|--------|------|---------|-------------|
| Species validity | `--filter-species` / `--no-filter-species` | ON | Skip records where `species` is empty, `None`, or other placeholder values |
| Kingdom | `--filter-kingdom` / `--no-filter-kingdom` | ON | Keep only records matching the specified kingdom(s) |
| Marker code | `--marker <CODE>` | `COI-5P` | Keep only records with the specified `marker_code` value |

Each filter logs how many records were skipped, giving a summary at the end of index building.

## Performance

### Benchmarks

Tested on typical HPC hardware (32 cores, 64GB RAM):

| Records | Taxa | Time (32 workers) | Time (1 worker) |
|---------|------|-------------------|-----------------|
| 300K | 50K | ~30 seconds | ~2 minutes |
| 1M | 100K | ~1 minute | ~5 minutes |
| 5M | 300K | ~3 minutes | ~15 minutes |

### Memory Usage

The script builds five in-memory indices during the records pass:
- `name_to_count`: record count per name
- `name_to_bins`: BIN/OTU set per name (primary analysis)
- `bin_to_names`: name set per BIN/OTU
- `name_to_bin_uris`: `bin_uri` set per name (output only)
- `name_to_otu_ids`: `otu_id` set per name (output only)
- `name_to_gb_count`: UK record count per name

For 5M records with 500K unique names and 100K BINs:
- Estimated memory: ~150–200MB for indices + working memory
- Recommended: 8GB+ for safety margin

### Optimization Tips

1. **Use BOLD mode filters** to reduce index size before analysis — filtering by marker and kingdom can cut records by 80%+ for a typical BOLD download
2. **Use parallel processing** for datasets >10K taxa
3. **Adjust batch size** based on memory: larger batches = less overhead, more memory
4. **Single-threaded mode** (`--workers 1`) uses less memory and avoids multiprocessing overhead for smaller datasets

## Troubleshooting

### Common Issues

**"Species list must have a 'taxon_name' or 'species' column"**
- Ensure your input file is tab-separated
- Check that the header row contains `taxon_name` or `species`

**"Records file must have 'bin_uri' or 'otu_id' column"**
- The records file needs a cluster identifier column
- Supported column names: `bin_uri`, `otu_id`, `OTU_ID`, `BIN`

**"Records file must have 'species' or 'organism' column"**
- The records file needs a species name column
- Supported column names: `species` (BOLD format), `organism` (NCBI format)

**Memory errors on large datasets**
- Reduce `--batch-size`
- Use fewer `--workers`
- In BOLD mode, apply stricter filters (`--marker`, `--kingdom-list`) to reduce the number of records indexed
- Ensure sufficient RAM allocation on HPC

**Encoding errors**
- The script automatically tries UTF-8 then Latin-1
- If issues persist, convert files to UTF-8

**Unexpected records in BOLD mode**
- Use `--log-level DEBUG` to see per-row filter decisions
- Check that the `marker_code` and `kingdom` column names match what's in your BOLD download

## SLURM Job Script (sh-run_gap_analysis.sh)

A ready-made SLURM batch submission script is provided for running the gap analysis on HPC clusters.

### SLURM Configuration

```bash
#SBATCH --partition=day
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=ALL
```

### Usage

Edit the path variables in the script to match your environment, then submit:

```bash
sbatch sh-run_gap_analysis.sh
```

### Variables to Configure

| Variable | Description |
|----------|-------------|
| `SCRIPT_DIR` | Directory containing `gap_analysis.py` |
| `SPECIES_LIST` | Path to your species list TSV |
| `RECORDS_FILE` | Path to your records TSV |
| `OUTPUT_FILE` | Desired output path |
| `WORKERS` | Number of parallel workers (should match `--cpus-per-task`) |

Adjust `--mem` and `--cpus-per-task` based on dataset size. See the Performance section above for guidance.

## Integration with Pipeline

This script is designed to work with BOLD data processing pipelines:

```bash
# Step 1: Download/prepare BOLD data
# ... (your data preparation steps)

# Step 2: Run gap analysis
python gap_analysis.py \
    --species-list target_species.tsv \
    --records BOLD_Public.tsv \
    --output gap_analysis.tsv \
    --workers 32 \
    --bold

# Step 3: Review results
# Filter for problem taxa:
# - RED status: taxonomic conflicts with named species — requires investigation
# - ORANGE status: BIN shared with provisional names only — lower priority
# - BLUE status: nomenclatural updates needed in BOLD (valid name absent)
# - BLACK status: sampling gaps to fill
# - AMBER status: synonym cleanup needed

# Step 4: Visualise results (optional)
# Use gap_analysis_reporting/gap_analysis_figures.py to produce publication-quality figures
# from the gap analysis output (requires pandas, matplotlib, numpy)
```

## Author

Ben Price / Claude
Date: 2025-01-25

## License

[Specify your license here]
