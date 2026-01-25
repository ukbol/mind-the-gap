# Pantheon to UKSI Name Matcher

Matches taxon names from the Pantheon database against the UK Species Inventory (UKSI) to retrieve taxonomic version keys.

## Purpose

This script takes a Pantheon database export and matches the taxon names against the UKSI to append:
- `TAXON_VERSION_KEY` - The TVK for the matched name
- `RECOMMENDED_TAXON_VERSION_KEY` - The TVK for the currently accepted name
- `NAME_STATUS` - Whether the match was to a Recommended (R) or Synonym (S) name

## Requirements

- Python 3.9+
- pandas

```bash
pip install pandas
```

## Usage

```bash
python pantheon_uksi_matcher.py --uksi <uksi_file> --pantheon <pantheon_file> --output-dir <output_directory>
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--uksi` | Yes | Path to UKSI names TSV file |
| `--pantheon` | Yes | Path to Pantheon input TSV file |
| `--output-dir` | Yes | Output directory for results |

### Example

```bash
python pantheon_uksi_matcher.py \
    --uksi ../uksi_20251203a_input_names.tsv \
    --pantheon pantheon_input_cleaned.tsv \
    --output-dir ./output
```

## Output Files

Output filenames are based on the input Pantheon filename:

| File | Description |
|------|-------------|
| `{input}_matched.tsv` | Pantheon records successfully matched to UKSI, with TVK columns appended |
| `{input}_unmatched.tsv` | Pantheon records not found in UKSI |
| `{input}_log.txt` | Processing log with statistics |

## Matching Logic

1. **Case-insensitive exact matching** - Taxon names are matched exactly (ignoring case)
2. **Fallback matching** - If `Taxon` column doesn't match, tries `original_taxon` column as fallback
3. **Recommended names preferred** - If a name exists as both a recommended name and synonym in UKSI, the recommended version is used
4. **Deduplication** - Duplicate Pantheon rows for the same taxon are combined, with values joined by commas

## Input File Requirements

### UKSI File
- Tab-separated values (TSV)
- Must contain columns: `TAXON_NAME`, `TAXON_VERSION_KEY`, `RECOMMENDED_TAXON_VERSION_KEY`, `NAME_STATUS`

### Pantheon File
- Tab-separated values (TSV)
- Must contain column: `Taxon`
- Optional column: `original_taxon` (used as fallback if `Taxon` doesn't match)

## Output Columns

The matched output file includes all original Pantheon columns plus:

| Column | Description |
|--------|-------------|
| `TAXON_VERSION_KEY` | The TVK for the matched name |
| `RECOMMENDED_TAXON_VERSION_KEY` | The TVK for the currently accepted name |
| `NAME_STATUS` | R=Recommended, S=Synonym, U=Unverified, I=Invalid |
| `MATCHED_ON` | Which column was used for matching: `Taxon` or `original_taxon` |

## Interpreting Results

### NAME_STATUS values

| Status | Meaning |
|--------|---------|
| R | Recommended - the currently accepted name |
| S | Synonym - a name that redirects to another accepted name |
| U | Unverified |
| I | Invalid |

If `NAME_STATUS = 'S'`, the `TAXON_VERSION_KEY` will differ from `RECOMMENDED_TAXON_VERSION_KEY`, indicating the matched name is a synonym.

### MATCHED_ON values

| Value | Meaning |
|-------|---------|
| Taxon | Matched on the primary `Taxon` column |
| original_taxon | Matched on fallback `original_taxon` column (typically subspecies with expanded names) |

## Author

Ben Price / Claude  
Natural History Museum, London
