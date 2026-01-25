# UKSI JNCC Annotation

Scripts for annotating UKSI (UK Species Inventory) species data with JNCC (Joint Nature Conservation Committee) conservation designations.

## Overview

These scripts match UKSI species records against JNCC conservation designation data by comparing Taxon Version Keys (TVKs). They support matching via both valid species TVKs and synonym TVKs to maximize coverage.

## Scripts

### uksi_jncc_annotation.py

Annotates UKSI valid species list with JNCC conservation designations by matching on the `Recommended_taxon_version` column.

**Usage:**
```bash
python uksi_jncc_annotation.py --input species.tsv --jncc jncc.tsv [--output-dir /path/to/output] [--stdout]
```

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input`, `-i` | Path to UKSI valid species TSV file |
| `--jncc`, `-j` | Path to JNCC conservation designations TSV file |
| `--output-dir`, `-o` | Output directory (default: same as input file) |
| `--stdout` | Write annotated TSV to stdout instead of file |
| `--log-file` | Custom log file path |

**Matching logic:**
- Checks valid species TVK (`taxon_version_key`) against JNCC `Recommended_taxon_version`
- Checks synonym TVKs (`synonym_tvk_list`) against JNCC TVKs
- Merges designations if multiple TVKs match (combining unique values with semicolons)
- Tracks match status: `valid`, `synonym`, or `valid;synonym`

**Output:**
- `uksi_valid_species_jncc_annotated.tsv` - Annotated species file
- `uksi_jncc_annotation.log` - Processing log

---

### uksi_jncc_annotation_v2.py

Enhanced version that matches against the JNCC `included_tvk_list` column (semicolon-separated TVKs).

**Usage:**
```bash
python uksi_jncc_annotation_v2.py --input uksi_species.tsv --jncc jncc_designations.tsv [--output-dir /path/to/output]
```

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input`, `-i` | Path to UKSI valid species TSV file |
| `--jncc`, `-j` | Path to JNCC conservation designations TSV file |
| `--output-dir`, `-o` | Output directory (default: same as input file) |
| `--log-file` | Custom log file path |

**Key differences from v1:**
- Matches against `included_tvk_list` column (semicolon-separated TVKs) rather than `Recommended_taxon_version`
- Handles cases where a JNCC row covers multiple TVKs
- Tracks which JNCC row indices were matched for reference
- Logs unmatched JNCC rows as errors for debugging

**Output columns added:**
| Column | Description |
|--------|-------------|
| `jncc_matching_tvk` | TVK(s) that matched |
| `tvk_match_status` | `valid`, `synonym`, or `valid;synonym` |
| `jncc_row_index` | JNCC file row number(s) matched |

---

### test_output.py

Utility script for inspecting annotation results. Displays sample species that were matched via different paths (synonym-only vs both valid and synonym).

**Usage:**
```bash
python test_output.py
```

**Note:** This script has hardcoded file paths and is intended for development/debugging purposes.

---

## JNCC Designation Columns

Both annotation scripts add the following conservation designation columns to the output:

| Column | Description |
|--------|-------------|
| A: Bern Convention | Bern Convention listing |
| C: Birds Directive | EU Birds Directive |
| C1: Convention on Migratory Species | CMS listing |
| C2: OSPAR | OSPAR Convention |
| D: Habitats Directive | EU Habitats Directive |
| E: EC Cites | CITES listing |
| F: Global Red list status | IUCN Global Red List |
| Fa: Red Listing based on pre 1994 IUCN guidelines | Pre-1994 Red List |
| Fb: Red Listing based on 1994 IUCN guidelines | 1994 Red List |
| Fc: Red listing based on 2001 IUCN guidelines | 2001 Red List |
| Fd: Red data categories - birds | Bird red data |
| Fe: Red data categories - Spiders | Spider red data |
| Ga: Rare and scarce species | Rare/scarce listing |
| Gb: Rare and scarce species (not IUCN) | Non-IUCN rare/scarce |
| Ha: Biodiversity Action Plan UK | UK BAP priority |
| Hb: Biodiversity Lists - England | England NERC S.41 |
| Hc: Biodiversity Lists - Scotland | Scottish Biodiversity List |
| Hd: Biodiversity Lists - Wales | Wales Environment Act S7 |
| He: Biodiversity Lists - Northern Ireland | NI Priority species |
| I: Wildlife and Countryside Act 1981 | WCA Schedule 5 |
| J: Wildlife (NI) Order 1985 | NI Wildlife Order |
| K: Conservation of Habitats Regulations 2010 | Habitats Regulations |
| L: Conservation Regulations (NI) 1995 | NI Habitats Regulations |
| M: Protection of Badgers Act | Badger protection |

## Input File Requirements

### UKSI Species File
- TSV format with header row
- Must contain `taxon_version_key` column
- Should contain `synonym_tvk_list` column (semicolon-separated TVKs)

### JNCC Designations File
- TSV format with header row
- For v1: `Recommended_taxon_version` column
- For v2: `included_tvk_list` column (semicolon-separated TVKs)
- All designation columns (A through M)

## Example Workflow

```bash
# Step 1: Run annotation with v2 script
python uksi_jncc_annotation_v2.py \
    --input ../uksi_valid_species_output.tsv \
    --jncc ../jncc_mapping/20231206_jncc_conservation_designations_taxon.tsv \
    --output-dir ./output

# Step 2: Check results
python test_output.py
```

## References

- [JNCC Conservation Designations](https://jncc.gov.uk/our-work/conservation-designations-for-uk-taxa/)
- [UKSI Nameserver](https://www.nhm.ac.uk/ukspecies/)
