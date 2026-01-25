# UKSI Database Pipeline

A SQLite-based pipeline for processing UK Species Inventory (UKSI) data, linking taxonomic names with Pantheon ecological traits and JNCC conservation designations.

## Overview

This pipeline creates a SQLite database from UKSI export files and generates a comprehensive species checklist with:
- Higher taxonomy (Kingdom → Genus)
- All Latin synonyms (including subspecific taxa and subgenus variants)
- Pantheon invertebrate ecological traits
- JNCC conservation designations (with hierarchical propagation)

## Directory Structure

```
uksi_db/
├── uksi_import.py              # Database creation script
├── uksi_export.py              # Species export script
├── uksi.db                     # SQLite database (~144 MB)
├── uksi_species_export.tsv     # Valid species output (~22 MB)
├── uksi_invalid_species_export.tsv  # Filtered invalid species
├── uksi_import.log             # Import log
├── uksi_export.log             # Export log
└── README.md                   # This file
```

## Input Files

The pipeline requires four input files located in the parent directory:

| File | Description | Source |
|------|-------------|--------|
| `uksi_20251203a_input_names.tsv` | UKSI Nameserver table (~337k names) | NHM UKSI export |
| `uksi_20251203a_input_taxa.tsv` | UKSI Taxa backbone (~124k taxa) | NHM UKSI export |
| `pantheon_mapping/output/pantheon_input_cleaned_matched.tsv` | Pantheon data matched to UKSI (~13k records) | Pre-processed |
| `jncc_mapping/20231206_jncc_conservation_designations_taxon.tsv` | JNCC designations (~14k records) | JNCC export |

## Quick Start

```bash
# 1. Create the database
python uksi_import.py

# 2. Export species checklist
python uksi_export.py
```

## Scripts

### uksi_import.py

Creates a SQLite database with four main tables and supporting indexes/views.

**Tables created:**
- `names` - All UKSI name variants (synonyms, vernaculars, misspellings)
- `taxa` - Backbone taxonomy with hierarchical structure
- `pantheon` - Invertebrate ecological traits
- `jncc` - Conservation designations (original)
- `jncc_resolved` - Conservation designations with TVK resolution

**Key features:**
- Resolves JNCC synonym TVKs through the names table (633 additional matches)
- Creates indexes for efficient lookups
- Creates utility views for common queries

**Output:**
```
Names table: 336,814 total, 124,397 unique recommended TVKs
Taxa table: 124,397 total, 79,522 species
Pantheon table: 12,745 total, 12,745 matched to taxa
JNCC table: 14,395 total
  - Direct match to taxa: 13,762
  - Resolved via names table: 14,395 (633 additional)
```


### uksi_export.py

Exports valid species with comprehensive data to a TSV file.

**Species filtering:**
- Includes: Species, Microspecies, Species hybrid, Species aggregate, Intergeneric hybrid, Species sensu lato/stricto
- Excludes: Names containing `.` (sp., spp., cf., aff.), `?`, `(Other)`, `"`, `(unidentified)`, `indet`, `/` (species aggregates like species1/species2)
- Invalid species are written to a separate file for review

**Synonym generation:**
- Includes all Latin names from the names table pointing to each species
- Includes child taxa from taxa table (subspecies, varieties, forms, etc.)
- For names like `Genus (Subgenus) species`, automatically adds:
  - `Genus species` (without subgenus)
  - `Subgenus species` (subgenus treated as genus)
- Deduplicates by name (each unique synonym appears once)
- Semicolon-separated

**JNCC designation propagation (CRITICAL):**
- **Downward inheritance:** If a higher taxon (e.g., Cetacea) has a designation, all descendant species inherit it
- **Upward inheritance:** If a subspecific taxon has a designation, the parent species inherits it

**Output:**
```
Valid species exported: 76,640
Invalid species filtered: 1,429
  - Contains '.': 1,390
  - Contains '?': 29
  - Contains 'indet': 4
  - Contains '(Other)': 4
  - Contains '(unidentified)': 2
```

## Output Columns

All column headers are lowercase.

### Core identifiers
| Column | Description |
|--------|-------------|
| organism_key | UKSI organism identifier |
| taxon_version_key | UKSI TVK (primary key for linking) |

### Higher taxonomy
| Column | Description |
|--------|-------------|
| kingdom | Kingdom (e.g., Animalia, Plantae) |
| phylum_division | Phylum or Division |
| class | Class |
| order | Order |
| family | Family |
| genus | Genus |

### Synonyms
| Column | Description |
|--------|-------------|
| synonyms | Semicolon-separated list of Latin synonyms |

### Taxon information
| Column | Description |
|--------|-------------|
| taxon_name | Scientific name |
| taxon_authority | Taxonomic authority |
| taxon_rank | Taxonomic rank (Species, Species aggregate, etc.) |

### Flags
| Column | Description |
|--------|-------------|
| non_native_flag | Non-native species indicator |
| terrestrial_freshwater_flag | Terrestrial/freshwater indicator |
| freshwater | Freshwater indicator |
| marine_flag | Marine indicator |


### Pantheon columns
| Column | Source column | Description |
|--------|---------------|-------------|
| pantheon_sqs | SQS | Species Quality Score |
| pantheon_conservation_status | Conservation status | Conservation status |
| pantheon_current_conservation_status | Current conservation status | Current conservation status |
| pantheon_designation_summary | Designation summary | Summary of designations |
| pantheon_larval_feeding_guild | Larval feeding guild | Larval feeding ecology |
| pantheon_adult_feeding_guild | Adult feeding guild | Adult feeding ecology |
| pantheon_broad_biotope | Broad biotope | Broad habitat type |
| pantheon_habitat | Habitat | Specific habitat |
| pantheon_resources | Resources | Resource requirements |
| pantheon_specific_assemblage_type | Specific assemblage type | Assemblage classification |
| pantheon_habitat_score | Habitat score | Habitat quality score |
| pantheon_associations | Associations | Species associations |

### JNCC designation columns (all lowercase)
| Column | Description |
|--------|-------------|
| jncc_a: bern convention | Bern Convention listing |
| jncc_c: birds directive | EU Birds Directive |
| jncc_c1: convention on migratory species | CMS listing |
| jncc_c2: ospar | OSPAR Convention |
| jncc_d: habitats directive | EU Habitats Directive |
| jncc_e: ec cites | CITES listing |
| jncc_f: global red list status | IUCN Global Red List |
| jncc_fa: red listing based on pre 1994 iucn guidelines | Pre-1994 Red List |
| jncc_fb: red listing based on 1994 iucn guidelines | 1994 Red List |
| jncc_fc: red listing based on 2001 iucn guidelines | 2001 Red List |
| jncc_fd: red data categories  - birds (not based on iucn criteria) | Bird red data |
| jncc_fe: red data categories - spiders (not based on iucn criteria) | Spider red data |
| jncc_ga: rare and scarce species | Rare/scarce listing |
| jncc_gb: rare and scarce species (not based on iucn criteria) | Non-IUCN rare/scarce |
| jncc_ha: biodiversity action plan uk list of priority species | UK BAP priority |
| jncc_hb: biodiversity lists - england | England NERC S.41 |
| jncc_hc: biodiversity lists - scotland | Scottish Biodiversity List |
| jncc_hd: biodiversity lists - wales | Wales Environment Act S7 |
| jncc_he: biodiversity lists - northern ireland | NI Priority species |
| jncc_i: wildlife and countryside act 1981 | WCA Schedule 5 |
| jncc_j: the wildlife (northern ireland) order 1985 | NI Wildlife Order |
| jncc_k: the conservation of habitats and species regulations 2010 | Habitats Regulations |
| jncc_l: the conservation (nature habitats, etc_) regulations (ni) 199 | NI Habitats Regulations |
| jncc_m: protection of badgers act | Badger protection |

## Database Schema

### Key relationships

```
names.RECOMMENDED_TAXON_VERSION_KEY → taxa.TAXON_VERSION_KEY
pantheon.RECOMMENDED_TAXON_VERSION_KEY → taxa.TAXON_VERSION_KEY  
jncc_resolved.resolved_tvk → taxa.TAXON_VERSION_KEY
taxa.PARENT_KEY → taxa.ORGANISM_KEY (hierarchy)
```

### Indexes
- `idx_names_recommended_tvk` - Fast synonym lookups
- `idx_names_language` - Filter by language
- `idx_taxa_parent_tvk` - Hierarchy traversal
- `idx_taxa_lineage` - Lineage-based queries
- `idx_pantheon_rec_tvk` - Pantheon joins
- `idx_jncc_resolved_tvk` - JNCC joins


## Example Queries

### Find all synonyms for a species
```sql
SELECT n.TAXON_NAME, n.TAXON_AUTHORITY
FROM names n
WHERE n.RECOMMENDED_TAXON_VERSION_KEY = 'NBNSYS0000007302'
  AND n.LANGUAGE = 'la'
ORDER BY n.TAXON_NAME;
```

### Get species with JNCC designations
```sql
SELECT t.TAXON_NAME, jr.*
FROM taxa t
JOIN jncc_resolved jr ON jr.resolved_tvk = t.TAXON_VERSION_KEY
WHERE t.RANK = 'Species'
  AND jr."Ha: Biodiversity Action Plan UK list of priority species" IS NOT NULL;
```

### Find Pantheon data for a species
```sql
SELECT t.TAXON_NAME, p.*
FROM taxa t
JOIN pantheon p ON p.RECOMMENDED_TAXON_VERSION_KEY = t.TAXON_VERSION_KEY
WHERE t.TAXON_NAME LIKE 'Bombus%';
```

### Traverse taxonomy hierarchy
```sql
WITH RECURSIVE ancestors AS (
    SELECT ORGANISM_KEY, TAXON_NAME, RANK, PARENT_KEY, 0 as level
    FROM taxa 
    WHERE TAXON_NAME = 'Bombus terrestris'
    
    UNION ALL
    
    SELECT t.ORGANISM_KEY, t.TAXON_NAME, t.RANK, t.PARENT_KEY, a.level + 1
    FROM taxa t
    JOIN ancestors a ON t.ORGANISM_KEY = a.PARENT_KEY
)
SELECT * FROM ancestors ORDER BY level;
```

## Technical Notes

### JNCC TVK Resolution
The JNCC file uses some older/synonym TVKs that don't directly match the taxa table. The import script resolves these through the names table:
1. First attempts direct match: `jncc.Recommended_taxon_version → taxa.TAXON_VERSION_KEY`
2. If no match, resolves via names: `jncc.Recommended_taxon_version → names.TAXON_VERSION_KEY → names.RECOMMENDED_TAXON_VERSION_KEY → taxa.TAXON_VERSION_KEY`

This resolves an additional 633 JNCC records (from 13,762 to 14,395 total).

### JNCC Designation Inheritance
Designations cascade through the taxonomic hierarchy:
- **Example:** Cetacea (Infraorder) has designations for Habitats Directive, EC CITES, Wildlife Act, etc.
- All whale species (Balaenoptera, Megaptera, etc.) automatically inherit these designations
- Species may also have their own additional direct designations

### Subgenus Synonym Generation
For names with subgenus notation like `Acartia (Acartiura) clausi`:
- The export automatically generates `Acartia clausi` and `Acartiura clausi` as synonyms
- This ensures records can be matched regardless of whether the subgenus was recorded

### Invalid Species Filtering
Species are filtered to the invalid output file if they contain:
- `.` - Catches sp., spp., cf., aff., n. sp., etc.
- `?` - Uncertain identifications
- `(Other)`, `(unidentified)`, `indet` - Indeterminate names
- `"` - Informal names

## Version History

- **2025-01-25 v2.1** - Reordered columns, added taxon_rank, all lowercase headers, exclude species aggregates with "/" pattern, removed redundant species column
- **2025-01-25 v2.0** - Added invalid species filtering, semicolon synonym separator, subgenus synonym generation, synonym deduplication
- **2025-01-25 v1.0** - Initial database pipeline with JNCC resolution and hierarchical propagation

## Author

Generated for Ben Price, Principal Curator, Natural History Museum London
UKBOL / Biodiversity Genomics Europe project

## References

- [UKSI Nameserver Guide](https://www.nhm.ac.uk/ukspecies/) - Chris Raper, NHM
- [NBN Atlas](https://nbnatlas.org/) - National Biodiversity Network
- [JNCC Conservation Designations](https://jncc.gov.uk/our-work/conservation-designations-for-uk-taxa/)
- [Pantheon Database](https://www.brc.ac.uk/pantheon/) - Biological Records Centre
