# BAGS Assessment

Calculate BAGS (Barcode, Audit & Grade System) grades for taxonomic barcode reference libraries based on OTU clustering results.

## Overview

BAGS grades assess the quality of species-level DNA barcode reference data by evaluating:
- Whether species are assigned to OTUs (clusters/BINs)
- Whether species share OTUs with other species (indicating potential misidentification or cryptic diversity)
- Whether species are split across multiple OTUs
- The number of specimens available for each species

This script is designed to run after OTU clustering (e.g., using `otu_clustering.py`) and produces a quality assessment for each taxonomic unit.

## BAGS Grades

| Grade | Description | Criteria |
|-------|-------------|----------|
| **A** | High quality | 1 species, 1 unshared OTU, ≥11 specimens |
| **B** | Good quality | 1 species, 1 unshared OTU, 3-10 specimens |
| **C** | Split species | 1 species split across multiple OTUs, each OTU unshared |
| **D** | Low coverage | 1 species, 1 unshared OTU, <3 specimens |
| **E** | BIN sharing | Multiple species share the same OTU (conflict) |
| **F** | No barcode | No OTU_ID assigned for this species |

**Grade Precedence** (highest to lowest): E > C > A > B > D > F

- Grade **E** indicates potential taxonomic or identification problems
- Grade **C** may indicate cryptic species or geographic variation
- Grades **A** and **B** represent well-supported species
- Grade **D** needs more sampling
- Grade **F** has no usable barcode data

## Installation

No additional dependencies required beyond Python 3.7+. The script uses only standard library modules.

```bash
# Clone or copy the script
cd /path/to/mind-the-gap/bags_assessment
```

## Usage

### Basic Usage

```bash
python bags_assessment.py input.tsv
```

### With Verbose Output

```bash
python bags_assessment.py -v input.tsv
```

### Specify Output Directory

```bash
python bags_assessment.py input.tsv -o /path/to/output/
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `input` | Input TSV file with OTU clustering results (required) |
| `-o, --output-dir` | Output directory (default: same as input file) |
| `-v, --verbose` | Print progress and statistics |

## Input Format

Tab-separated file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `species` or `organism` | Yes | Species name |
| `OTU_ID` | Yes | OTU cluster assignment (from clustering step) |
| `accession` or `processid` | Yes | Record identifier |
| `taxid` | No | Taxonomic ID (will be generated if missing) |

**Notes:**
- Rows with blank species names are ignored
- Species matching is case-insensitive
- If `taxid` column is missing, unique taxids are generated for each unique species name

### Example Input

```tsv
processid	species	OTU_ID	...
PROC001	Homo sapiens	OTU_000001	...
PROC002	Homo sapiens	OTU_000001	...
PROC003	Pan troglodytes	OTU_000002	...
```

## Output Files

### 1. assessed_BAGS.tsv

BAGS assessment results with columns:

| Column | Description |
|--------|-------------|
| `taxid` | Taxonomic ID |
| `BAGS` | BAGS grade (A-F) |
| `OTU_ID` | OTU ID(s), pipe-separated if multiple |
| `sharers` | Other species sharing each OTU, comma-separated within OTU, pipe-separated between OTUs |

**Example Output:**
```tsv
taxid	BAGS	OTU_ID	sharers
1	A	OTU_000001	
2	B	OTU_000002	
3	C	OTU_000003|OTU_000004	|
4	E	OTU_000005	Species B,Species C
5	F		
```

### 2. {input}_taxid.tsv

Original input data with `taxid` column appended (if it was missing from the input).

## Pipeline Integration

This script is designed to run after OTU clustering in the mind-the-gap pipeline:

```bash
# Step 1: OTU clustering
python otu_clustering/otu_clustering.py -t 0.99 input.tsv clustered.tsv

# Step 2: BAGS assessment
python bags_assessment/bags_assessment.py -v clustered.tsv
```

### HPC Usage

The script is memory-efficient and processes data in a streaming fashion where possible. For very large datasets:

```bash
# Run with verbose output to monitor progress
python bags_assessment.py -v large_dataset.tsv -o /scratch/output/
```

## Interpreting Results

### Grade Distribution

A healthy reference library should have:
- High proportion of grades A and B
- Low proportion of grades E (conflicts)
- Some grade C may be expected for species complexes

### Sharers Column

For grade E records, the `sharers` column lists other species in the same OTU:
- `Species A,Species B` = two other species share this OTU
- For multiple OTUs: `Species A|Species B,Species C` means OTU1 has Species A, OTU2 has Species B and C

### Actions Based on Grades

| Grade | Recommended Action |
|-------|-------------------|
| A, B | No action needed |
| C | Review - may indicate cryptic diversity |
| D | Add more specimens if possible |
| E | Review identifications - potential errors |
| F | Obtain barcode sequences |

## References

BAGS grading system adapted from:
- Barcode Index Numbers (BINs) - Ratnasingham & Hebert, 2013
- DNA barcode reference library curation workflows

## Author

Ben Price / Claude  
Date: 2025-01-23
