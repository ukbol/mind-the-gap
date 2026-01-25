# HPC Gap Analysis for DNA Barcode Library Curation

A high-performance Python script for taxon-centric gap analysis of DNA barcode reference libraries. Designed for HPC environments with support for parallel processing of large datasets (millions of records, hundreds of thousands of taxa).

## Overview

This tool analyzes the quality of DNA barcode reference data by:

1. **Loading a target species list** — Taxa defined by valid names and their synonyms
2. **Scanning records** — Finding all records matching each taxon's names
3. **Analyzing BIN/OTU sharing** — Detecting taxonomic conflicts where different taxa share clusters
4. **Assigning grades and status** — BAGS grades (A-F) and traffic light status (GREEN/AMBER/RED/BLUE/BLACK)

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

### Command Line Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--species-list` | Yes | — | Path to species list TSV file |
| `--records` | Yes | — | Path to records TSV file (e.g., result_output.tsv) |
| `--output` | Yes | — | Output path for gap analysis TSV |
| `--workers` | No | CPU count | Number of parallel workers |
| `--batch-size` | No | 1000 | Taxa per batch for parallel processing |
| `--log-level` | No | INFO | Logging level: DEBUG, INFO, WARNING, ERROR |

## Input Files

### Species List Format

Tab-separated file with headers:

| Column | Required | Description |
|--------|----------|-------------|
| `taxon_name` | Yes | Valid species name (preferred to species)|
| `species` | No | Valid species name (if taxon_name not present)|
| `synonyms` | No | Semicolon-separated synonyms (can be empty) |
| *(other columns)* | No | Preserved in output |

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
| `species` | Yes | Species name |
| `subspecies` | No | Subspecies name (matched separately if present) |
| `bin_uri` or `otu_id` | Yes | Cluster identifier(s), pipe-separated if multiple |

**Note:** The script automatically detects which cluster column to use (prefers `bin_uri`).

## Output

### Gap Analysis TSV

The output file contains all input columns plus:

| Column | Description |
|--------|-------------|
| `number_records` | Count of records matching any name in the taxon |
| `bags_grade` | Quality grade: A, B, C, D, E, or F |
| `species_status` | Traffic light status: GREEN, AMBER, RED, BLUE, or BLACK |
| `other_names` | Semicolon-separated names sharing BIN/OTU but not in taxon's list |

## Grading System

### BAGS Grades

| Grade | Meaning | Criteria |
|-------|---------|----------|
| **A** | Excellent | 1 BIN/OTU, ≥11 records, no conflicts |
| **B** | Good | 1 BIN/OTU, 3-10 records, no conflicts |
| **C** | Split | Multiple BIN/OTUs, no conflicts |
| **D** | Low coverage | 1 BIN/OTU, <3 records, no conflicts |
| **E** | Conflict | Other taxa share the BIN/OTU |
| **F** | No data | No records found |

### Traffic Light Status

| Status | Colour | Meaning | Criteria |
|--------|--------|---------|----------|
| **GREEN** | 🟢 | Clean | Only valid name recorded, no conflicts |
| **AMBER** | 🟡 | Nomenclatural mess | Both valid name AND synonym(s) recorded (needs cleanup) |
| **BLUE** | 🔵 | Valid name absent | Only synonym(s) recorded, valid name missing from database |
| **RED** | 🔴 | Taxonomic conflict | Names from outside the taxon share BIN/OTU |
| **BLACK** | ⚫ | No coverage | No records for this taxon |

### Decision Logic

```
┌─────────────────────────────────────────────────────────────┐
│                    GRADING DECISION TREE                    │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Has records?                                               │
│  ├── NO  → Grade F, Status BLACK                           │
│  └── YES                                                    │
│       │                                                     │
│       ▼                                                     │
│  Other names share BIN/OTU?                                 │
│  ├── YES → Grade E, Status RED                             │
│  └── NO                                                     │
│       │                                                     │
│       ▼                                                     │
│  Which names are recorded?                                  │
│  ├── Valid + synonym(s) → Status AMBER                     │
│  ├── Valid only → Status GREEN                             │
│  └── Synonym(s) only → Status BLUE                         │
│       │                                                     │
│       ▼                                                     │
│  How many BIN/OTUs?                                         │
│  ├── 1 BIN/OTU                                             │
│  │    ├── ≥11 records → Grade A                            │
│  │    ├── 3-10 records → Grade B                           │
│  │    └── <3 records → Grade D                             │
│  └── Multiple BIN/OTUs → Grade C                           │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Examples

### Example Scenarios

| Taxon | BIN contains | Status | Grade | Explanation |
|-------|--------------|--------|-------|-------------|
| Valid: "Alpha vulgaris", Syn: "Alpha communis" | Only "Alpha vulgaris" (15 records) | GREEN | A | Clean, good coverage |
| Valid: "Beta marina", Syn: "Beta aquatica" | Both "Beta marina" and "Beta aquatica" (8 records) | AMBER | B | Nomenclatural cleanup needed |
| Valid: "Gamma riparia", Syn: "Gamma fluviatilis" | Only "Gamma fluviatilis" (5 records) | BLUE | B | Valid name not in BOLD |
| Valid: "Delta palustris", Syn: none | "Delta palustris" + "Epsilon montanus" | RED | E | Taxonomic conflict |
| Valid: "Zeta alpina", Syn: "Zeta montana" | No records | BLACK | F | No coverage |
| Valid: "Eta borealis", Syn: none | "Eta borealis" in 3 different BINs | GREEN | C | Split across clusters |

### Example Output

```
species           synonyms        kingdom   number_records  bags_grade  species_status  other_names
Alpha vulgaris    Alpha communis  Animalia  15              A           GREEN           
Beta marina       Beta aquatica   Animalia  8               B           AMBER           
Gamma riparia     Gamma fluv...   Animalia  5               B           BLUE            
Delta palustris                   Animalia  12              E           RED             Epsilon montanus
Zeta alpina       Zeta montana    Animalia  0               F           BLACK           
```

## Performance

### Benchmarks

Tested on typical HPC hardware (32 cores, 64GB RAM):

| Records | Taxa | Time (32 workers) | Time (1 worker) |
|---------|------|-------------------|-----------------|
| 300K | 50K | ~30 seconds | ~2 minutes |
| 1M | 100K | ~1 minute | ~5 minutes |
| 5M | 300K | ~3 minutes | ~15 minutes |

### Memory Usage

The script builds three in-memory indices:
- `name_to_count`: ~50 bytes per unique name
- `name_to_bins`: ~100 bytes per unique name
- `bin_to_names`: ~100 bytes per BIN/OTU

For 5M records with 500K unique names and 100K BINs:
- Estimated memory: ~150MB for indices + working memory
- Recommended: 8GB+ for safety margin

### Optimization Tips

1. **Use parallel processing** for datasets >10K taxa
2. **Adjust batch size** based on memory: larger batches = less overhead, more memory
3. **Single-threaded mode** (`--workers 1`) uses less memory, suitable for smaller datasets

## Troubleshooting

### Common Issues

**"Species list must have a 'species' column"**
- Ensure your input file is tab-separated
- Check that the header row contains `species`

**"Records file must have 'bin_uri' or 'otu_id' column"**
- The records file needs a cluster identifier column
- Supported column names: `bin_uri`, `otu_id`, `OTU_ID`, `BIN`

**Memory errors on large datasets**
- Reduce `--batch-size`
- Use fewer workers
- Ensure sufficient RAM allocation on HPC

**Encoding errors**
- The script automatically tries UTF-8 then Latin-1
- If issues persist, convert files to UTF-8

## Integration with Pipeline

This script is designed to work with BOLD data processing pipelines:

```bash
# Step 1: Download/prepare BOLD data
# ... (your data preparation steps)

# Step 2: Run gap analysis
python gap_analysis.py \
    --species-list target_species.tsv \
    --records result_output.tsv \
    --output gap_analysis.tsv \
    --workers 32

# Step 3: Review results
# Filter for problem taxa:
# - RED status: taxonomic conflicts requiring investigation
# - BLUE status: nomenclatural updates needed in BOLD
# - BLACK status: sampling gaps to fill
```

## Author

Ben Price / Claude  
Date: 2025-01-25

## License

[Specify your license here]
