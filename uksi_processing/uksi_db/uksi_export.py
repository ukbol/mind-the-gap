#!/usr/bin/env python3
"""
UKSI Database Export Script
===========================
Exports valid species with:
- Higher taxonomy (Kingdom -> Genus)
- All Latin synonyms (including subspecific taxa, subgenus variants)
- Pantheon ecological traits
- JNCC conservation designations (propagated through hierarchy)
- FreshBase and UKCEH freshwater list membership flags

Author: Generated for Ben Price, NHM London
Date: 2025-01-25
Version: 2.2

Changes in v2.2:
- Add FreshBase and UKCEH freshwater list presence columns

Changes in v2.1:
- Filter to only Animalia, Plantae, Chromista, and Fungi kingdoms
- Exclude Species aggregate, Species sensu lato, and Species hybrid ranks from valid output

Changes in v2:
- Filter invalid species (sp., spp., cf., ?, etc.)
- Use semicolon separator for synonyms
- Add subgenus-derived synonyms for Genus (Subgenus) species names
- Deduplicate synonyms while being inclusive

CRITICAL: JNCC designations are propagated:
- Downward: If a higher taxon (e.g., Cetacea) has a designation, all child species inherit it
- Upward: If a subspecific taxon has a designation, the parent species inherits it
"""

import sqlite3
import csv
import sys
import re
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# Configuration
BASE_DIR = Path(r"C:\GitHub\mind-the-gap\uksi_processing")
DB_PATH = BASE_DIR / "uksi_db" / "uksi.db"
OUTPUT_PATH = BASE_DIR / "uksi_db" / "uksi_species_export.tsv"
INVALID_OUTPUT_PATH = BASE_DIR / "uksi_db" / "uksi_invalid_species_export.tsv"
LOG_PATH = BASE_DIR / "uksi_db" / "uksi_export.log"

# Taxonomic ranks for higher taxonomy extraction (in order)
HIGHER_RANKS = ['Kingdom', 'Phylum', 'Division', 'Class', 'Order', 'Family', 'Genus']

# Species-level ranks to include in output
SPECIES_RANKS = [
    'Species', 'Microspecies', 'Intergeneric hybrid', 'Species sensu stricto'
]

# Kingdoms to include in output
VALID_KINGDOMS = ['Animalia', 'Plantae', 'Chromista', 'Fungi']

# Ranks to include in synonyms (species + all infraspecific)
SYNONYM_RANKS = {
    'Species', 'Species aggregate', 'Species group', 'Species hybrid',
    'Species pro parte', 'Species sensu lato', 'Species sensu stricto',
    'Subspecies', 'Subspecies aggregate', 'Subspecies hybrid',
    'Variety', 'Varietal hybrid', 'Subvariety',
    'Form', 'Subform', 'Nothosubspecies', 'Nothovariety',
    'Microspecies', 'Praespecies', 'Cultivar', 'Convariety',
    'Abberation', 'Breed', 'Facies', 'Morphotype', 'Race',
    'Pathovar', 'Forma specialis', 'ecad',
}

# Patterns that indicate an invalid/indeterminate species name
INVALID_NAME_PATTERNS = [
    r'\.',                  # Any period (catches sp., spp., cf., aff., n. sp, etc.)
    r'\?',                  # Question mark
    r'\(Other\)',           # (Other)
    r'"',                   # Double quotes
    r'\(unidentified\)',    # (unidentified)
    r'\bindet\b',           # indet (indeterminate)
    r'/',                   # Slash (catches species aggregates like "species1/species2")
]
INVALID_NAME_REGEX = re.compile('|'.join(INVALID_NAME_PATTERNS), re.IGNORECASE)

# Regex to detect subgenus pattern: Genus (Subgenus) species [subspecies...]
SUBGENUS_PATTERN = re.compile(r'^(\w+)\s+\((\w+)\)\s+(.+)$')


# JNCC columns to export
JNCC_COLUMNS = [
    ("A: Bern Convention", "jncc_A: Bern Convention"),
    ("C: Birds Directive", "jncc_C: Birds Directive"),
    ("C1: Convention on Migratory Species", "jncc_C1: Convention on Migratory Species"),
    ("C2: OSPAR", "jncc_C2: OSPAR"),
    ("D: Habitats Directive", "jncc_D: Habitats Directive"),
    ("E: EC Cites", "jncc_E: EC Cites"),
    ("F: Global Red list status", "jncc_F: Global Red list status"),
    ("Fa: Red Listing based on pre 1994 IUCN guidelines", "jncc_Fa: Red Listing based on pre 1994 IUCN guidelines"),
    ("Fb: Red Listing based on 1994 IUCN guidelines", "jncc_Fb: Red Listing based on 1994 IUCN guidelines"),
    ("Fc: Red listing based on 2001 IUCN guidelines", "jncc_Fc: Red listing based on 2001 IUCN guidelines"),
    ("Fd: Red data categories  - birds (not based on IUCN criteria)", "jncc_Fd: Red data categories  - birds (not based on IUCN criteria)"),
    ("Fe: Red data categories - Spiders (not based on IUCN criteria)", "jncc_Fe: Red data categories - Spiders (not based on IUCN criteria)"),
    ("Ga: Rare and scarce species", "jncc_Ga: Rare and scarce species"),
    ("Gb: Rare and scarce species (not based on IUCN criteria)", "jncc_Gb: Rare and scarce species (not based on IUCN criteria)"),
    ("Ha: Biodiversity Action Plan UK list of priority species", "jncc_Ha: Biodiversity Action Plan UK list of priority species"),
    ("Hb: Biodiversity Lists - England", "jncc_Hb: Biodiversity Lists - England"),
    ("Hc: Biodiversity Lists - Scotland", "jncc_Hc: Biodiversity Lists - Scotland"),
    ("Hd: Biodiversity Lists - Wales", "jncc_Hd: Biodiversity Lists - Wales"),
    ("He: Biodiversity Lists - Northern Ireland", "jncc_He: Biodiversity Lists - Northern Ireland"),
    ("I: Wildlife and Countryside Act 1981", "jncc_I: Wildlife and Countryside Act 1981"),
    ("J: The Wildlife (Northern Ireland) Order 1985", "jncc_J: The Wildlife (Northern Ireland) Order 1985"),
    ("K: The Conservation of Habitats and Species Regulations 2010", "jncc_K: The Conservation of Habitats and Species Regulations 2010"),
    ("L: The Conservation (Nature Habitats, etc_) Regulations (NI) 199", "jncc_L: The Conservation (Nature Habitats, etc_) Regulations (NI) 199"),
    ("M: Protection of Badgers Act", "jncc_M: Protection of Badgers Act"),
]

# Pantheon columns to export
PANTHEON_COLUMNS = [
    ("SQS", "pantheon_sqs"),
    ("Conservation status", "pantheon_conservation_status"),
    ("Current conservation status", "pantheon_current_conservation_status"),
    ("Designation summary", "pantheon_designation_summary"),
    ("Larval feeding guild", "pantheon_larval_feeding_guild"),
    ("Adult feeding guild", "pantheon_adult_feeding_guild"),
    ("Broad biotope", "pantheon_broad_biotope"),
    ("Habitat", "pantheon_habitat"),
    ("Resources", "pantheon_resources"),
    ("Specific assemblage type", "pantheon_specific_assemblage_type"),
    ("Habitat score", "pantheon_habitat_score"),
    ("Associations", "pantheon_associations"),
]


def log(message: str):
    """Log message to both console and file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] {message}"
    print(log_line)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(log_line + "\n")


def is_invalid_species_name(name: str) -> tuple:
    """
    Check if a species name contains patterns indicating it's indeterminate.
    Returns (is_invalid, reason) tuple.
    """
    if INVALID_NAME_REGEX.search(name):
        match = INVALID_NAME_REGEX.search(name)
        return True, f"Contains '{match.group()}'"
    return False, None


def extract_subgenus_synonyms(species_name: str) -> list:
    """
    For names like "Genus (Subgenus) species", extract additional synonyms:
    - "Genus species" (without subgenus)
    - "Subgenus species" (subgenus as genus)
    
    Returns list of additional synonym strings, or empty list if not applicable.
    """
    match = SUBGENUS_PATTERN.match(species_name)
    if match:
        genus = match.group(1)
        subgenus = match.group(2)
        epithet = match.group(3)  # May include subspecies etc.
        
        return [
            f"{genus} {epithet}",       # Genus species (without subgenus)
            f"{subgenus} {epithet}",    # Subgenus species
        ]
    return []


def build_lineage_lookup(conn: sqlite3.Connection) -> dict:
    """
    Build a lookup from ORGANISM_KEY to taxon info for hierarchy traversal.
    Returns dict: organism_key -> {tvk, name, authority, rank, parent_key, lineage}
    """
    log("Building lineage lookup table...")
    cur = conn.cursor()
    cur.execute("""
        SELECT ORGANISM_KEY, TAXON_VERSION_KEY, TAXON_NAME, TAXON_AUTHORITY, 
               RANK, PARENT_KEY, LINEAGE
        FROM taxa
    """)
    
    lookup = {}
    for row in cur.fetchall():
        lookup[row[0]] = {
            'tvk': row[1],
            'name': row[2],
            'authority': row[3],
            'rank': row[4],
            'parent_key': row[5],
            'lineage': row[6]
        }
    
    log(f"  Built lookup with {len(lookup):,} taxa")
    return lookup


def get_higher_taxonomy(organism_key: str, lineage_lookup: dict) -> dict:
    """
    Traverse up the taxonomy to extract higher ranks.
    Returns dict with keys: kingdom, phylum_division, class, order, family, genus
    """
    result = {
        'kingdom': '',
        'phylum_division': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': ''
    }
    
    # Walk up the tree via parent_key
    current_key = organism_key
    visited = set()
    
    while current_key and current_key in lineage_lookup and current_key not in visited:
        visited.add(current_key)
        taxon = lineage_lookup[current_key]
        rank = taxon['rank']
        name = taxon['name']
        
        if rank == 'Kingdom':
            result['kingdom'] = name
        elif rank in ('Phylum', 'Division'):
            result['phylum_division'] = name
        elif rank == 'Class':
            result['class'] = name
        elif rank == 'Order':
            result['order'] = name
        elif rank == 'Family':
            result['family'] = name
        elif rank == 'Genus':
            result['genus'] = name
        
        current_key = taxon['parent_key']
    
    return result


def build_jncc_designation_maps(conn: sqlite3.Connection, lineage_lookup: dict) -> dict:
    """
    Build JNCC designation maps that propagate designations through the hierarchy.
    
    CRITICAL LOGIC:
    1. If a higher taxon (e.g., Cetacea order) has a designation, all descendant species inherit it
    2. If a subspecific taxon has a designation, the parent species inherits it
    
    Returns: dict[tvk] -> dict[jncc_column] -> value
    """
    log("Building JNCC designation maps with hierarchy propagation...")
    cur = conn.cursor()
    
    # Get all JNCC designations with their resolved TVKs
    jncc_col_names = [col[0] for col in JNCC_COLUMNS]
    jncc_col_sql = ', '.join([f'"{col}"' for col in jncc_col_names])
    
    cur.execute(f"""
        SELECT resolved_tvk, {jncc_col_sql}
        FROM jncc_resolved
        WHERE resolved_tvk IS NOT NULL
    """)
    
    # Direct designations: tvk -> {column -> value}
    direct_designations = defaultdict(dict)
    for row in cur.fetchall():
        tvk = row[0]
        for i, col_name in enumerate(jncc_col_names):
            val = row[i + 1]
            if val and val.strip():
                direct_designations[tvk][col_name] = val.strip()
    
    log(f"  Loaded {len(direct_designations):,} taxa with direct JNCC designations")
    
    # Build TVK -> organism_key mapping
    cur.execute("SELECT TAXON_VERSION_KEY, ORGANISM_KEY FROM taxa")
    tvk_to_org = {}
    org_to_tvk = {}
    for row in cur.fetchall():
        tvk_to_org[row[0]] = row[1]
        org_to_tvk[row[1]] = row[0]
    
    # Get all species TVKs
    cur.execute(f"""
        SELECT TAXON_VERSION_KEY, LINEAGE, ORGANISM_KEY 
        FROM taxa 
        WHERE RANK IN ({','.join(['?' for _ in SPECIES_RANKS])})
    """, SPECIES_RANKS)
    species_data = cur.fetchall()
    
    log(f"  Processing {len(species_data):,} species for designation inheritance...")
    
    # Final designations map (includes inherited)
    final_designations = defaultdict(dict)
    
    # For each species, walk up parent chain and collect designations
    inheritance_count = 0
    for species_tvk, species_lineage, species_org in species_data:
        # Start with direct designations
        if species_tvk in direct_designations:
            final_designations[species_tvk].update(direct_designations[species_tvk])
        
        # Walk up parent chain
        current_org = species_org
        visited = set()
        while current_org and current_org in lineage_lookup and current_org not in visited:
            visited.add(current_org)
            parent_org = lineage_lookup[current_org]['parent_key']
            
            if parent_org and parent_org in org_to_tvk:
                parent_tvk = org_to_tvk[parent_org]
                if parent_tvk in direct_designations:
                    # Inherit designations from ancestor
                    for col, val in direct_designations[parent_tvk].items():
                        if col not in final_designations[species_tvk]:
                            final_designations[species_tvk][col] = val
                            inheritance_count += 1
            
            current_org = parent_org
    
    log(f"  Inherited {inheritance_count:,} designations from higher taxa")
    
    # Now handle subspecific taxa propagating UP to species
    cur.execute("""
        SELECT TAXON_VERSION_KEY, PARENT_KEY, ORGANISM_KEY
        FROM taxa 
        WHERE RANK IN ('Subspecies', 'Variety', 'Form', 'Subvariety', 'Breed')
    """)
    subspecific_data = cur.fetchall()
    
    upward_count = 0
    for sub_tvk, parent_key, sub_org in subspecific_data:
        if sub_tvk in direct_designations and parent_key:
            if parent_key in org_to_tvk:
                parent_tvk = org_to_tvk[parent_key]
                for col, val in direct_designations[sub_tvk].items():
                    if col not in final_designations[parent_tvk]:
                        final_designations[parent_tvk][col] = val
                        upward_count += 1
    
    log(f"  Propagated {upward_count:,} designations upward from subspecific taxa")
    log(f"  Final: {len(final_designations):,} taxa with JNCC designations")
    
    return final_designations


def get_latin_synonyms(conn: sqlite3.Connection, lineage_lookup: dict) -> dict:
    """
    Get all Latin synonyms for each species TVK.
    
    Includes:
    - All Latin names from names table pointing to the recommended TVK
    - Child taxa from taxa table (subspecies, varieties, etc.)
    - Subgenus-derived synonyms (Genus species, Subgenus species)
    
    Deduplicates by name (keeping unique names only).
    
    Returns: dict[recommended_tvk] -> list of unique synonym name strings
    """
    log("Loading Latin synonyms (comprehensive)...")
    cur = conn.cursor()
    
    # Build org_to_tvk mapping
    cur.execute("SELECT TAXON_VERSION_KEY, ORGANISM_KEY, TAXON_NAME, TAXON_AUTHORITY FROM taxa")
    org_to_tvk = {}
    tvk_to_info = {}
    for row in cur.fetchall():
        org_to_tvk[row[1]] = row[0]
        tvk_to_info[row[0]] = {'name': row[2], 'authority': row[3]}
    
    # Get all Latin names from names table
    cur.execute("""
        SELECT 
            n.RECOMMENDED_TAXON_VERSION_KEY,
            n.TAXON_NAME,
            n.TAXON_AUTHORITY,
            n.TAXON_VERSION_KEY,
            n.RANK
        FROM names n
        WHERE n.LANGUAGE = 'la'
    """)
    
    # Collect all names per recommended TVK
    synonyms_raw = defaultdict(list)
    for row in cur.fetchall():
        rec_tvk = row[0]
        name = row[1]
        authority = row[2] if row[2] else ''
        tvk = row[3]
        rank = row[4]
        
        # Include name with authority
        if authority:
            full_name = f"{name} {authority}"
        else:
            full_name = name
        
        synonyms_raw[rec_tvk].append({
            'name': name,
            'full_name': full_name,
            'tvk': tvk,
            'rank': rank
        })
    
    log(f"  Loaded names for {len(synonyms_raw):,} taxa from names table")
    
    # Build child taxa index by PARENT_TVK for subspecies/varieties/forms
    cur.execute("""
        SELECT TAXON_VERSION_KEY, PARENT_TVK, TAXON_NAME, TAXON_AUTHORITY, RANK
        FROM taxa
        WHERE PARENT_TVK IS NOT NULL AND PARENT_TVK != ''
    """)
    
    children_by_parent = defaultdict(list)
    for row in cur.fetchall():
        child_tvk = row[0]
        parent_tvk = row[1]
        name = row[2]
        authority = row[3] if row[3] else ''
        rank = row[4]
        
        if authority:
            full_name = f"{name} {authority}"
        else:
            full_name = name
        
        children_by_parent[parent_tvk].append({
            'name': name,
            'full_name': full_name,
            'tvk': child_tvk,
            'rank': rank
        })
    
    log(f"  Built child index for {len(children_by_parent):,} parent taxa")
    
    # Now build final synonym lists with deduplication
    final_synonyms = {}
    
    for rec_tvk, entries in synonyms_raw.items():
        # Get the valid name info
        valid_info = tvk_to_info.get(rec_tvk, {})
        valid_name = valid_info.get('name', '')
        valid_authority = valid_info.get('authority', '')
        if valid_authority:
            valid_full = f"{valid_name} {valid_authority}"
        else:
            valid_full = valid_name
        
        # Collect all synonym names (just the name part, not full with authority for dedup)
        seen_names = set()
        synonym_names = []
        
        # Add names from names table
        for entry in entries:
            name = entry['name']
            # Skip if it's the valid name itself
            if name == valid_name:
                continue
            if name not in seen_names:
                seen_names.add(name)
                synonym_names.append(name)
        
        # Add child taxa from taxa table
        for child in children_by_parent.get(rec_tvk, []):
            name = child['name']
            if name != valid_name and name not in seen_names:
                seen_names.add(name)
                synonym_names.append(name)
        
        # Add subgenus-derived synonyms if valid name has subgenus pattern
        subgenus_syns = extract_subgenus_synonyms(valid_name)
        for syn in subgenus_syns:
            if syn != valid_name and syn not in seen_names:
                seen_names.add(syn)
                synonym_names.append(syn)
        
        # Sort alphabetically
        synonym_names.sort()
        
        final_synonyms[rec_tvk] = synonym_names
    
    total_syns = sum(len(v) for v in final_synonyms.values())
    log(f"  Final: {len(final_synonyms):,} taxa with {total_syns:,} unique synonyms")
    
    return final_synonyms


def get_pantheon_data(conn: sqlite3.Connection) -> dict:
    """
    Load Pantheon ecological trait data.
    Returns: dict[recommended_tvk] -> dict of column values
    """
    log("Loading Pantheon data...")
    cur = conn.cursor()
    
    col_names = [col[0] for col in PANTHEON_COLUMNS]
    col_sql = ', '.join([f'"{col}"' for col in col_names])
    
    cur.execute(f"""
        SELECT RECOMMENDED_TAXON_VERSION_KEY, {col_sql}
        FROM pantheon
    """)
    
    pantheon = {}
    for row in cur.fetchall():
        tvk = row[0]
        data = {}
        for i, col_name in enumerate(col_names):
            val = row[i + 1]
            data[col_name] = val if val else ''
        pantheon[tvk] = data
    
    log(f"  Loaded Pantheon data for {len(pantheon):,} taxa")
    return pantheon


def get_freshwater_presence(conn: sqlite3.Connection) -> tuple:
    """
    Build sets of TVKs that appear in each freshwater species list.
    Uses the resolved tables (which map synonym TVKs to recommended TVKs).

    Returns: (freshbase_tvks, ukceh_freshwater_tvks) as sets of TVK strings
    """
    log("Loading freshwater list data...")
    cur = conn.cursor()

    cur.execute("SELECT DISTINCT resolved_tvk FROM freshbase_resolved WHERE resolved_tvk IS NOT NULL")
    freshbase_tvks = {row[0] for row in cur.fetchall()}
    log(f"  FreshBase: {len(freshbase_tvks):,} resolved TVKs")

    cur.execute("SELECT DISTINCT resolved_tvk FROM ukceh_freshwater_resolved WHERE resolved_tvk IS NOT NULL")
    ukceh_tvks = {row[0] for row in cur.fetchall()}
    log(f"  UKCEH freshwater list: {len(ukceh_tvks):,} resolved TVKs")

    return freshbase_tvks, ukceh_tvks


def export_species(conn: sqlite3.Connection):
    """Main export function."""
    log("\n=== Starting Species Export ===")
    
    # Build lookup tables
    lineage_lookup = build_lineage_lookup(conn)
    jncc_designations = build_jncc_designation_maps(conn, lineage_lookup)
    synonyms = get_latin_synonyms(conn, lineage_lookup)
    pantheon = get_pantheon_data(conn)
    freshbase_tvks, ukceh_tvks = get_freshwater_presence(conn)
    
    # Get all valid species
    log("\nQuerying valid species...")
    cur = conn.cursor()
    cur.execute(f"""
        SELECT 
            ORGANISM_KEY,
            TAXON_VERSION_KEY,
            TAXON_NAME,
            TAXON_AUTHORITY,
            RANK,
            NON_NATIVE_FLAG,
            TERRESTRIAL_FRESHWATER_FLAG,
            FRESHWATER,
            MARINE_FLAG
        FROM taxa
        WHERE RANK IN ({','.join(['?' for _ in SPECIES_RANKS])})
        AND (REDUNDANT_FLAG IS NULL OR REDUNDANT_FLAG = '')
        ORDER BY TAXON_NAME
    """, SPECIES_RANKS)
    
    species_rows = cur.fetchall()
    log(f"  Found {len(species_rows):,} species (before filtering)")
    
    # Build output headers (all lowercase, reordered as specified)
    headers = [
        # Core identifiers
        'organism_key',
        'taxon_version_key',
        # Higher taxonomy
        'kingdom',
        'phylum_division',
        'class',
        'order',
        'family',
        'genus',
        # Taxon info
        'taxon_name',
        'taxon_authority',
        'taxon_rank',
        # Synonyms
        'synonyms',
        # Flags
        'non_native_flag',
        'terrestrial_freshwater_flag',
        'freshwater',
        'marine_flag',
    ]
    
    # Add freshwater list columns
    headers.append('freshbase')
    headers.append('ukceh_freshwater_list')

    # Add Pantheon columns
    for _, output_name in PANTHEON_COLUMNS:
        headers.append(output_name)

    # Add JNCC columns (lowercase)
    for _, output_name in JNCC_COLUMNS:
        headers.append(output_name.lower())
    
    # Write output
    log(f"\nWriting valid species to {OUTPUT_PATH}...")
    log(f"Writing invalid species to {INVALID_OUTPUT_PATH}...")
    
    valid_count = 0
    invalid_count = 0
    kingdom_filtered_count = 0
    invalid_reasons = defaultdict(int)
    
    with open(OUTPUT_PATH, 'w', newline='', encoding='utf-8') as f_valid, \
         open(INVALID_OUTPUT_PATH, 'w', newline='', encoding='utf-8') as f_invalid:
        
        writer_valid = csv.writer(f_valid, delimiter='\t')
        writer_invalid = csv.writer(f_invalid, delimiter='\t')
        
        writer_valid.writerow(headers)
        writer_invalid.writerow(headers)
        
        for species in species_rows:
            org_key = species[0]
            tvk = species[1]
            taxon_name = species[2]
            taxon_authority = species[3] or ''
            rank = species[4]
            non_native = species[5] or ''
            terr_fw = species[6] or ''
            freshwater = species[7] or ''
            marine = species[8] or ''
            
            # Get higher taxonomy first (needed for kingdom filter)
            higher_tax = get_higher_taxonomy(org_key, lineage_lookup)
            
            # Filter by kingdom - skip if not in valid kingdoms
            if higher_tax['kingdom'] not in VALID_KINGDOMS:
                kingdom_filtered_count += 1
                continue
            
            # Check if name is invalid
            is_invalid, reason = is_invalid_species_name(taxon_name)
            
            # Get synonyms (semicolon separated)
            syn_list = synonyms.get(tvk, [])
            syn_str = ';'.join(syn_list)
            
            # Get Pantheon data
            panth = pantheon.get(tvk, {})
            
            # Get JNCC designations
            jncc = jncc_designations.get(tvk, {})
            
            # Build output row (matching new column order)
            row = [
                # Core identifiers
                org_key,
                tvk,
                # Higher taxonomy
                higher_tax['kingdom'],
                higher_tax['phylum_division'],
                higher_tax['class'],
                higher_tax['order'],
                higher_tax['family'],
                higher_tax['genus'],
                # Taxon info
                taxon_name,
                taxon_authority,
                rank,
                # Synonyms
                syn_str,
                # Flags
                non_native,
                terr_fw,
                freshwater,
                marine,
            ]
            
            # Add freshwater list presence
            row.append('Y' if tvk in freshbase_tvks else '')
            row.append('Y' if tvk in ukceh_tvks else '')

            # Add Pantheon columns
            for input_col, _ in PANTHEON_COLUMNS:
                row.append(panth.get(input_col, ''))

            # Add JNCC columns
            for input_col, _ in JNCC_COLUMNS:
                row.append(jncc.get(input_col, ''))
            
            # Write to appropriate file
            if is_invalid:
                writer_invalid.writerow(row)
                invalid_count += 1
                invalid_reasons[reason] += 1
            else:
                writer_valid.writerow(row)
                valid_count += 1
            
            if (valid_count + invalid_count) % 10000 == 0:
                log(f"  Processed {valid_count + invalid_count:,} species...")
    
    log(f"\n  Valid species exported: {valid_count:,}")
    log(f"  Invalid species filtered: {invalid_count:,}")
    log(f"  Kingdom filtered (not Animalia/Plantae/Chromista/Fungi): {kingdom_filtered_count:,}")
    for reason, count in sorted(invalid_reasons.items(), key=lambda x: -x[1]):
        log(f"    - {reason}: {count:,}")
    
    return valid_count, invalid_count


def validate_export(conn: sqlite3.Connection, valid_count: int, invalid_count: int):
    """Validate the export with summary statistics."""
    log("\n=== Export Validation ===")
    
    # Read back the valid export file and check stats
    with open(OUTPUT_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    
    log(f"Total valid rows in export: {len(rows):,}")
    
    # Count records with various data
    with_synonyms = sum(1 for r in rows if r.get('synonyms', ''))
    with_pantheon = sum(1 for r in rows if r.get('pantheon_sqs', ''))
    with_jncc = sum(1 for r in rows if any(r.get(col[1].lower(), '') for col in JNCC_COLUMNS))
    with_freshbase = sum(1 for r in rows if r.get('freshbase', ''))
    with_ukceh = sum(1 for r in rows if r.get('ukceh_freshwater_list', ''))

    log(f"  - With synonyms: {with_synonyms:,}")
    log(f"  - With Pantheon data: {with_pantheon:,}")
    log(f"  - With JNCC designations: {with_jncc:,}")
    log(f"  - On FreshBase list: {with_freshbase:,}")
    log(f"  - On UKCEH freshwater list: {with_ukceh:,}")
    
    # Check higher taxonomy coverage
    with_kingdom = sum(1 for r in rows if r.get('kingdom', ''))
    with_phylum = sum(1 for r in rows if r.get('phylum_division', ''))
    with_class = sum(1 for r in rows if r.get('class', ''))
    with_order = sum(1 for r in rows if r.get('order', ''))
    with_family = sum(1 for r in rows if r.get('family', ''))
    with_genus = sum(1 for r in rows if r.get('genus', ''))
    
    log(f"  Higher taxonomy coverage:")
    log(f"    - Kingdom: {with_kingdom:,}")
    log(f"    - Phylum/Division: {with_phylum:,}")
    log(f"    - Class: {with_class:,}")
    log(f"    - Order: {with_order:,}")
    log(f"    - Family: {with_family:,}")
    log(f"    - Genus: {with_genus:,}")
    
    # Check subgenus synonym generation
    subgenus_count = 0
    for r in rows:
        name = r.get('taxon_name', '')
        syns = r.get('synonyms', '')
        if SUBGENUS_PATTERN.match(name) and syns:
            subgenus_count += 1
    log(f"\n  Species with subgenus pattern: {subgenus_count:,}")
    
    # Sample a subgenus species to verify synonym generation
    log("\n  Sample subgenus species with synonyms:")
    count = 0
    for r in rows:
        name = r.get('taxon_name', '')
        syns = r.get('synonyms', '')
        if SUBGENUS_PATTERN.match(name) and syns:
            log(f"    {name}")
            log(f"      Synonyms: {syns[:100]}...")
            count += 1
            if count >= 3:
                break


def main():
    """Main entry point for the export script."""
    log("=" * 60)
    log("UKSI Database Export Script v2.2")
    log("=" * 60)
    
    # Clear log file
    if LOG_PATH.exists():
        LOG_PATH.unlink()
    
    # Check database exists
    if not DB_PATH.exists():
        log(f"ERROR: Database not found: {DB_PATH}")
        log("Run uksi_import.py first to create the database.")
        sys.exit(1)
    
    log(f"\nDatabase: {DB_PATH}")
    log(f"Valid output: {OUTPUT_PATH}")
    log(f"Invalid output: {INVALID_OUTPUT_PATH}")
    
    # Connect to database
    conn = sqlite3.connect(str(DB_PATH))
    
    try:
        # Export species
        valid_count, invalid_count = export_species(conn)
        
        # Validate
        validate_export(conn, valid_count, invalid_count)
        
        log("\n" + "=" * 60)
        log("Export completed successfully!")
        log(f"Valid species: {OUTPUT_PATH}")
        log(f"Invalid species: {INVALID_OUTPUT_PATH}")
        log("=" * 60)
        
    except Exception as e:
        log(f"\nERROR: {e}")
        import traceback
        log(traceback.format_exc())
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
