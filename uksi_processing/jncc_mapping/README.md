# JNCC Conservation Designations to UKSI TVK Mapper

Maps JNCC conservation designation taxa to UK Species Inventory (UKSI) taxon version keys (TVKs), expanding coverage based on taxonomic rank.

## Purpose

Conservation designations in the JNCC file are assigned to specific taxa. This script expands that coverage:

- **Above species level** (e.g., genus, family): The designation applies to all known species within that group, so all descendant TVKs are included
- **Below species level** (e.g., subspecies, variety): The designation should also protect the parent species, so the species-level TVK is included
- **Species level**: Only the taxon itself (and its recommended TVK if different) is included

## Input Files

| File | Description | Key Columns |
|------|-------------|-------------|
| JNCC designations TSV | Conservation designations with taxon references | `Recommended_taxon_version` |
| UKSI names TSV | Taxon names with TVK mappings | `TAXON_VERSION_KEY`, `RECOMMENDED_TAXON_VERSION_KEY`, `RANK` |
| UKSI taxa TSV | Taxonomic hierarchy | `TAXON_VERSION_KEY`, `PARENT_TVK`, `RANK` |

## Output Files

- **`{input_filename}_uksi_mapped.tsv`**: Original JNCC data with added `included_tvk_list` column
- **`jncc_mapping_{timestamp}.log`**: Detailed processing log with warnings and statistics

## Usage

### Default paths (as configured for NHM workflow):


```bash
python jncc_uksi_mapper.py
```

### Custom paths:

```bash
python jncc_uksi_mapper.py \
    --jncc /path/to/jncc_designations.tsv \
    --names /path/to/uksi_names.tsv \
    --taxa /path/to/uksi_taxa.tsv \
    --output-dir /path/to/output/
```

### Command Line Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--jncc` | `-j` | `../20231206_jncc_conservation_designations_taxon.tsv` | JNCC designations file |
| `--names` | `-n` | `../uksi_20251203a_input_names.tsv` | UKSI names file |
| `--taxa` | `-t` | `../uksi_20251203a_input_taxa.tsv` | UKSI taxa hierarchy file |
| `--output-dir` | `-o` | Current directory | Output directory |

## Output Format

The output TSV contains all original columns from the JNCC file plus:

| Column | Description |
|--------|-------------|
| `included_tvk_list` | Semicolon-separated list of all related TVKs |

### What's included in `included_tvk_list`:

1. The original JNCC taxon's TVK (always)
2. The recommended TVK if different from original (always)
3. **If above species**: All descendant TVKs + their recommended TVKs
4. **If below species**: Parent species TVK + its recommended TVK
5. **If species level**: Just items 1 and 2


## Rank Classification

Taxa are classified into three categories based on their RANK value:

### Above Species
Genus, Family, Order, Class, Phylum, Kingdom, and other higher ranks. Also includes aggregates (Species aggregate, Genus aggregate, etc.) and Species sensu lato.

### Species Level
- Species
- Species hybrid
- Unranked (treated as species)

### Below Species
Subspecies, Variety, Form, Cultivar, Race, Breed, and other infraspecific ranks.

### Ambiguous Ranks
When rank is "Unknown" or cannot be determined, the script infers the category by examining the hierarchy:
- If parent is at species level → treat as below species
- If any descendant is at species level → treat as above species
- Otherwise → treat as species level

## Error Handling

- **Unmatched TVKs**: Logged as warnings, included in output with empty `included_tvk_list`
- **Empty TVKs**: Logged as warnings, included in output with empty `included_tvk_list`
- **Encoding issues**: Automatically falls back through UTF-8 → UTF-8-BOM → Latin-1 → Windows-1252

## Dependencies

- Python 3.6+
- No external packages required (uses only standard library)

## Example

Input JNCC row:
```
Taxon_name: Aricia
Recommended_taxon_version: NHMSYS0021366585
(... other columns ...)
```

If Aricia is a genus containing species A, B, C with TVKs TVK_A, TVK_B, TVK_C:

Output row:
```
included_tvk_list: NHMSYS0021366585;TVK_A;TVK_A_REC;TVK_B;TVK_C
```


## Log File Contents

The log file contains:
- Processing parameters and file paths
- Row counts for each input file
- Warnings for unmatched or problematic taxa
- Debug information about rank inference (at DEBUG level)
- Summary statistics

## Author

Generated for NHM London biodiversity genomics work (Principal Curator, Freshwater Entomology & Biodiversity Genomics).

## Version History

- **v1.0** (2025-01): Initial version with hierarchical TVK expansion
