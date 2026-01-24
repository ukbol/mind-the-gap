# UKSI Species Extraction Pipeline

## Overview

This pipeline extracts valid species names from the UK Species Inventory (UKSI) database, along with their higher taxonomy and synonyms. It produces two output files: one for valid species and one for invalid/indeterminate species.

## Input Files

1. **TAXA table** (`uksi_20251203a_input_taxa.tsv`)
   - Contains the hierarchical taxonomy tree
   - Each row represents a taxon (any rank from Kingdom to infraspecific)
   - Key fields:
     - `ORGANISM_KEY`: Unique identifier for the taxon
     - `TAXON_VERSION_KEY` (TVK): Unique identifier for this name version
     - `PARENT_KEY` / `PARENT_TVK`: Links to parent taxon in hierarchy
     - `RANK`: Taxonomic rank (Kingdom, Phylum, Class, Order, Family, Genus, Species, Subspecies, Variety, etc.)
     - `TAXON_NAME`: The scientific name
     - `REDUNDANT_FLAG`: 'Y' = dubious/redundant, '' = valid (NOT used for filtering - see note below)

2. **NAMES table** (`uksi_20251203a_input_names.tsv`)
   - Contains all name versions and their mappings to recommended names
   - Key fields:
     - `TAXON_VERSION_KEY`: TVK for this specific name
     - `TAXON_NAME`: The scientific name
     - `RECOMMENDED_TAXON_VERSION_KEY`: TVK of the accepted/valid name this resolves to
     - `RECOMMENDED_NAME_RANK`: Rank of the recommended name
     - `NAME_STATUS`: 'R' = recommended, 'S' = synonym, 'U' = unverified
     - `DEPRECATED_DATE`: If populated, name is deprecated
     - `LANGUAGE`: 'la' = Latin/scientific, 'en' = English, etc.

## Processing Logic

### Step 1: Load TAXA table
- Index all taxa by `ORGANISM_KEY` and `TAXON_VERSION_KEY`
- Total records: ~124,000

### Step 2: Build child taxa index
- For each taxon, index it under its `PARENT_TVK`
- This allows finding all subspecies/varieties/forms that belong to a species
- Includes ALL children (redundant or not) because even redundant infraspecific names should appear as synonyms of their parent species

### Step 3: Identify valid species
- From TAXA, select all records where `RANK = 'Species'`
- **Note:** We do NOT filter on `REDUNDANT_FLAG` (see "Why no REDUNDANT_FLAG filter?" below)
- Total species: ~79,500

### Step 4: Load NAMES table for synonyms
- Build index: `RECOMMENDED_TAXON_VERSION_KEY` → list of synonym entries
- **Filters applied:**
  - `DEPRECATED_DATE = ''` (not deprecated)
  - `LANGUAGE = 'la'` (scientific names only)
  - `RANK` in species/infraspecific ranks
  - `TAXON_VERSION_KEY ≠ RECOMMENDED_TAXON_VERSION_KEY` (exclude the valid name's own entry)
  - Skip if TVK is itself a valid species (conflict case)

### Step 5: Process each valid species
For each species identified in Step 3:

1. **Filter by Kingdom:** Only include Animalia, Plantae, Fungi, Chromista
2. **Filter invalid names:** Names containing `.`, `?`, `"`, `(Other)`, `(unidentified)`, `indet` go to invalid output
3. **Build higher taxonomy:** Traverse up via `PARENT_KEY` to get Kingdom, Phylum/Division, Class, Order, Family, Genus
4. **Collect synonyms from three sources:**
   - **NAMES table:** All names pointing to this TVK via `RECOMMENDED_TAXON_VERSION_KEY`
   - **TAXA children:** All child taxa (subspecies, varieties, forms) via `PARENT_TVK`
   - **Subgenus expansion:** For names like "Genus (Subgenus) species", generate "Subgenus species" and "Genus species"
5. **Filter synonyms:** Remove any identical to the valid species name
6. **Output:** Write to valid or invalid TSV

## Why no REDUNDANT_FLAG filter?

Analysis revealed that **4,669 species** are marked `REDUNDANT_FLAG = 'Y'` in TAXA but are still referenced as valid recommended names in NAMES (i.e., other names have `RECOMMENDED_TAXON_VERSION_KEY` pointing to them).

This suggests that `REDUNDANT_FLAG` indicates something like "this taxon concept has been superseded or is uncertain" rather than "this is not a valid species name". Since synonyms still resolve to these species, they must be included in the output.

**Evidence:**
- TAXA non-redundant species: 74,853
- TAXA all species (including redundant): 79,522
- NAMES unique species TVKs: 79,522
- Redundant in TAXA but valid in NAMES: 4,669

All species TVKs referenced in NAMES exist in TAXA, so using all species from TAXA captures the complete set.

## Output Files

### Valid species output (`uksi_valid_species_output.tsv`)
| Column | Description |
|--------|-------------|
| organism_key | From TAXA |
| taxon_version_key | TVK of the valid species |
| kingdom | From hierarchy traversal |
| phylum_division | Phylum or Division (combined) |
| class | From hierarchy traversal |
| order | From hierarchy traversal |
| family | From hierarchy traversal |
| genus | From hierarchy traversal |
| species | TAXON_NAME from TAXA |
| synonyms | Semicolon-delimited list of synonym names |
| synonym_tvk_list | Semicolon-delimited list of synonym TVKs (aligned with synonyms) |
| recommended_name_authority | TAXON_AUTHORITY from TAXA |
| non_native_flag | From TAXA |
| terrestrial_freshwater_flag | From TAXA |
| freshwater | From TAXA |
| marine_flag | From TAXA |

### Invalid species output (`uksi_invalid_species_output.tsv`)
Same columns as valid output, but contains species with indeterminate names:
- Names containing `.` (catches sp., spp., cf., aff., n. sp., etc.)
- Names containing `?`
- Names containing `"`
- Names containing `(Other)`
- Names containing `(unidentified)`
- Names containing `indet`

### Log file (`uksi_pipeline_log.txt`)
- Processing statistics
- List of excluded identical synonyms (names matching the valid species that were removed from synonym list)
- List of valid-as-synonym conflicts (for database debugging)
- List of invalid species by reason

## Synonym Sources

Synonyms are collected from three sources:

1. **NAMES table linkage:** All names where `RECOMMENDED_TAXON_VERSION_KEY` = species TVK
   - Excludes deprecated names (`DEPRECATED_DATE` not empty)
   - Excludes vernacular names (only `LANGUAGE = 'la'`)
   - Excludes the valid name's own entry
   - Includes species-level and infraspecific synonyms

2. **TAXA hierarchy (child taxa):** All taxa where `PARENT_TVK` = species TVK
   - Captures subspecies, varieties, forms, etc. that are children of the species
   - Includes redundant children (they should still appear as synonyms)
   - This is essential because infraspecific taxa often point to themselves in NAMES, not to their parent species

3. **Subgenus expansion:** For species names like "Genus (Subgenus) epithet"
   - Generates "Subgenus epithet" as a synonym
   - Generates "Genus epithet" as a synonym
   - These derived synonyms have empty TVKs in the synonym_tvk_list

## File Locations

- Script: `uksi_species_pipeline_v2.py`
- Input TAXA: `uksi_20251203a_input_taxa.tsv`
- Input NAMES: `uksi_20251203a_input_names.tsv`
- Output valid: `uksi_valid_species_output.tsv`
- Output invalid: `uksi_invalid_species_output.tsv`
- Log: `uksi_pipeline_log.txt`

## Reference Documentation

- Guide PDF: `uksi_guide.pdf` - Nameserver structure and SQL examples
- Glossary: `uksi_sandbox_glossary.tsv` - Taxonomy terminology
- NBN Summary: `uksi_summary_nbn.txt` - Overview of UK Species Inventory
