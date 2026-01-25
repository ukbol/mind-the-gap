"""
UKSI Species Extraction Pipeline v2
====================================
Extracts valid species names from UKSI TAXA and NAMES tables,
with higher taxonomy and synonyms.

Author: Claude (for Ben Price, NHM)
Date: 2026-01-24
Version: 2.0

Changes in v2:
- Exclude synonyms identical to the valid name (logged)
- Fix broken lines in authority field (handle embedded newlines)
- Strip whitespace from synonym list, use semicolon without space
- Filter invalid species (sp., cf., aff., etc.) to separate output
- Add subgenus-derived synonyms for names like Genus (Subgenus) species
- Check for names that are both valid and synonym (logged for DB fix)
"""

import csv
import sys
import re
from collections import defaultdict
from datetime import datetime

# Configuration
TAXA_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_20251203a_input_taxa.tsv'
NAMES_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_20251203a_input_names.tsv'
OUTPUT_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_valid_species_output.tsv'
INVALID_OUTPUT_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_invalid_species_output.tsv'
LOG_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_pipeline_log.txt'

# Kingdoms to include
TARGET_KINGDOMS = {'Animalia', 'Plantae', 'Fungi', 'Chromista'}

# Ranks that count as "species level" for the main output
SPECIES_RANKS = {'Species', 'Microspecies', 'Praespecies'}

# Ranks to include in synonyms (species + all infraspecific)
SYNONYM_RANKS = {
    'Species',
    'Species aggregate',
    'Species group', 
    'Species hybrid',
    'Species pro parte',
    'Species sensu lato',
    'Subspecies',
    'Subspecies aggregate',
    'Subspecies hybrid',
    'Variety',
    'Varietal hybrid',
    'Subvariety',
    'Form',
    'Subform',
    'Nothosubspecies',
    'Nothovariety',
    'Microspecies',
    'Praespecies',
    'Cultivar',
    'Convariety',
    'Abberation',
    'Breed',
    'Facies',
    'Morphotype',
    'Race',
    'Pathovar',
    'Forma specialis',
    'ecad',
}

# Higher taxonomy ranks we want to extract (in order of hierarchy)
HIGHER_RANKS = ['Kingdom', 'Phylum', 'Division', 'Class', 'Order', 'Family', 'Genus']

# Patterns that indicate an invalid/indeterminate species name
# Note: Any name containing a period (.) is considered invalid, as this typically
# indicates abbreviations like sp., spp., cf., aff., n. sp, etc.
INVALID_NAME_PATTERNS = [
    r'\.',                  # Any period in name (catches sp., spp., cf., aff., n. sp, etc.)
    r'\?',                  # ?
    r'\(Other\)',           # (Other)
    r'"',                   # "
    r'\(unidentified\)',    # (unidentified)
    r'\bindet\b',           # indet (indeterminate)
]

# Compile into single regex
INVALID_NAME_REGEX = re.compile('|'.join(INVALID_NAME_PATTERNS), re.IGNORECASE)

# Regex to detect subgenus pattern: Genus (Subgenus) species [subspecies...]
SUBGENUS_PATTERN = re.compile(r'^(\w+)\s+\((\w+)\)\s+(.+)$')


class PipelineLogger:
    """Logger for pipeline events."""
    
    def __init__(self, log_path):
        self.log_path = log_path
        self.log_entries = []
        self.excluded_identical_synonyms = []
        self.valid_as_synonym_conflicts = []
        self.invalid_species = []
        
    def log(self, message):
        """Add a log message."""
        self.log_entries.append(message)
        print(message)
    
    def log_excluded_identical_synonym(self, species_name, tvk, synonym_name, synonym_tvk):
        """Log when a synonym is excluded because it matches the valid name."""
        self.excluded_identical_synonyms.append({
            'species': species_name,
            'species_tvk': tvk,
            'synonym': synonym_name,
            'synonym_tvk': synonym_tvk,
        })
    
    def log_valid_as_synonym_conflict(self, valid_name, valid_tvk, synonym_tvk, rec_tvk):
        """Log when a name appears as both valid and synonym."""
        self.valid_as_synonym_conflicts.append({
            'valid_name': valid_name,
            'valid_tvk': valid_tvk,
            'synonym_tvk': synonym_tvk,
            'recommended_tvk': rec_tvk,
        })
    
    def log_invalid_species(self, species_name, tvk, reason):
        """Log when a species is filtered as invalid."""
        self.invalid_species.append({
            'species': species_name,
            'tvk': tvk,
            'reason': reason,
        })
    
    def write_log(self):
        """Write all log entries to file."""
        with open(self.log_path, 'w', encoding='utf-8') as f:
            f.write("UKSI Species Extraction Pipeline Log\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write("=" * 80 + "\n\n")
            
            # General log entries
            f.write("PROCESSING LOG\n")
            f.write("-" * 40 + "\n")
            for entry in self.log_entries:
                f.write(entry + "\n")
            f.write("\n")
            
            # Excluded identical synonyms
            f.write("=" * 80 + "\n")
            f.write(f"EXCLUDED IDENTICAL SYNONYMS ({len(self.excluded_identical_synonyms)} entries)\n")
            f.write("-" * 40 + "\n")
            f.write("These synonyms were excluded because the name matches the valid species name.\n\n")
            for entry in self.excluded_identical_synonyms[:100]:  # Limit to first 100
                f.write(f"Species: {entry['species']} (TVK: {entry['species_tvk']})\n")
                f.write(f"  Excluded synonym: {entry['synonym']} (TVK: {entry['synonym_tvk']})\n")
            if len(self.excluded_identical_synonyms) > 100:
                f.write(f"\n... and {len(self.excluded_identical_synonyms) - 100} more\n")
            f.write("\n")
            
            # Valid as synonym conflicts
            f.write("=" * 80 + "\n")
            f.write(f"VALID NAME ALSO LISTED AS SYNONYM ({len(self.valid_as_synonym_conflicts)} entries)\n")
            f.write("-" * 40 + "\n")
            f.write("These names appear as both a valid species AND as a synonym pointing to another taxon.\n")
            f.write("The valid form is kept, and the synonym reference is removed.\n")
            f.write("This may indicate a database issue that needs manual review.\n\n")
            for entry in self.valid_as_synonym_conflicts:
                f.write(f"Valid name: {entry['valid_name']}\n")
                f.write(f"  Valid TVK: {entry['valid_tvk']}\n")
                f.write(f"  Also appears as synonym TVK: {entry['synonym_tvk']}\n")
                f.write(f"  Points to recommended TVK: {entry['recommended_tvk']}\n")
                f.write("\n")
            
            # Invalid species summary
            f.write("=" * 80 + "\n")
            f.write(f"INVALID SPECIES FILTERED ({len(self.invalid_species)} entries)\n")
            f.write("-" * 40 + "\n")
            f.write("These species were moved to the invalid output file due to indeterminate names.\n\n")
            # Group by reason
            by_reason = defaultdict(list)
            for entry in self.invalid_species:
                by_reason[entry['reason']].append(entry['species'])
            for reason, species_list in sorted(by_reason.items()):
                f.write(f"Reason: {reason} ({len(species_list)} species)\n")
                for sp in species_list[:10]:
                    f.write(f"  - {sp}\n")
                if len(species_list) > 10:
                    f.write(f"  ... and {len(species_list) - 10} more\n")
                f.write("\n")


def clean_field(value):
    """Clean a field value - remove embedded newlines and strip whitespace."""
    if value:
        # Replace newlines with spaces
        value = value.replace('\n', ' ').replace('\r', ' ')
        # Collapse multiple spaces
        value = ' '.join(value.split())
        return value.strip()
    return value


def load_taxa_with_clean_fields(filepath):
    """
    Load TAXA table with field cleaning to handle embedded newlines.
    Returns dicts indexed by ORGANISM_KEY and TAXON_VERSION_KEY.
    """
    taxa_by_org_key = {}
    taxa_by_tvk = {}
    
    with open(filepath, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Clean all fields
            cleaned_row = {k: clean_field(v) for k, v in row.items()}
            taxa_by_org_key[cleaned_row['ORGANISM_KEY']] = cleaned_row
            taxa_by_tvk[cleaned_row['TAXON_VERSION_KEY']] = cleaned_row
    
    return taxa_by_org_key, taxa_by_tvk


def load_synonyms_index(filepath, valid_species_tvks, logger):
    """
    Load NAMES table and build index of synonyms by RECOMMENDED_TAXON_VERSION_KEY.
    
    Includes ALL Latin names linking to species, regardless of deprecation status.
    
    Also builds a set of valid species names (from TAXA) to check for conflicts.
    
    Returns:
        synonyms_by_rec_tvk: dict mapping rec_tvk -> list of synonym names
        valid_name_tvks: set of TVKs that are valid species names
    """
    synonyms_by_rec_tvk = defaultdict(list)
    
    # First pass: collect all name entries
    all_name_entries = []
    with open(filepath, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Clean fields
            cleaned_row = {k: clean_field(v) for k, v in row.items()}
            all_name_entries.append(cleaned_row)
    
    # Build set of valid TVKs (from TAXA table)
    valid_tvk_set = set(valid_species_tvks.keys())
    
    # Check for conflicts: names that are both valid AND appear as synonyms elsewhere
    for row in all_name_entries:
        tvk = row['TAXON_VERSION_KEY']
        rec_tvk = row['RECOMMENDED_TAXON_VERSION_KEY']
        
        # If this TVK is a valid species but points to a different recommended TVK
        if tvk in valid_tvk_set and rec_tvk and tvk != rec_tvk:
            valid_name = valid_species_tvks[tvk]['TAXON_NAME']
            logger.log_valid_as_synonym_conflict(valid_name, tvk, tvk, rec_tvk)
    
    # Second pass: build synonym index
    for row in all_name_entries:
        # Filter criteria
        # NOTE: DEPRECATED_DATE filter removed - include all names regardless of deprecation
        if row['LANGUAGE'] != 'la':
            continue  # Skip non-scientific names
        if row['RANK'] not in SYNONYM_RANKS:
            continue  # Skip non-species/infraspecific ranks
        if row['TAXON_VERSION_KEY'] == row['RECOMMENDED_TAXON_VERSION_KEY']:
            continue  # Skip the valid name's own entry
        
        # Skip if this TVK is itself a valid species (conflict case)
        if row['TAXON_VERSION_KEY'] in valid_tvk_set:
            continue
        
        rec_tvk = row['RECOMMENDED_TAXON_VERSION_KEY']
        if rec_tvk:
            synonyms_by_rec_tvk[rec_tvk].append({
                'name': row['TAXON_NAME'],
                'tvk': row['TAXON_VERSION_KEY'],
            })
    
    return synonyms_by_rec_tvk


def get_higher_taxonomy(species_row, taxa_by_org_key):
    """
    Traverse up the hierarchy from a species to extract higher taxonomy.
    Returns dict with Kingdom, Phylum, Division, Class, Order, Family, Genus.
    """
    hierarchy = {rank: '' for rank in HIGHER_RANKS}
    current = species_row
    visited = set()
    
    while current and current['ORGANISM_KEY'] not in visited:
        visited.add(current['ORGANISM_KEY'])
        rank = current['RANK']
        
        if rank in hierarchy:
            hierarchy[rank] = current['TAXON_NAME']
        
        parent_key = current.get('PARENT_KEY', '')
        if parent_key and parent_key in taxa_by_org_key:
            current = taxa_by_org_key[parent_key]
        else:
            break
    
    return hierarchy


def is_invalid_species_name(name):
    """
    Check if a species name contains patterns indicating it's indeterminate.
    Returns (is_invalid, reason) tuple.
    """
    if INVALID_NAME_REGEX.search(name):
        match = INVALID_NAME_REGEX.search(name)
        return True, f"Contains '{match.group()}'"
    return False, None


def extract_subgenus_synonyms(species_name):
    """
    For names like "Genus (Subgenus) species", extract additional synonyms:
    - "Subgenus species" 
    - "Genus species"
    
    Returns list of additional synonym strings, or empty list if not applicable.
    """
    match = SUBGENUS_PATTERN.match(species_name)
    if match:
        genus = match.group(1)
        subgenus = match.group(2)
        epithet = match.group(3)  # May include subspecies etc.
        
        additional = [
            f"{subgenus} {epithet}",   # Subgenus species
            f"{genus} {epithet}",       # Genus species
        ]
        return additional
    return []


def process_species():
    """Main processing function."""
    
    logger = PipelineLogger(LOG_FILE)
    logger.log(f"Pipeline started: {datetime.now().isoformat()}")
    logger.log(f"Input TAXA: {TAXA_FILE}")
    logger.log(f"Input NAMES: {NAMES_FILE}")
    logger.log(f"Output valid: {OUTPUT_FILE}")
    logger.log(f"Output invalid: {INVALID_OUTPUT_FILE}")
    logger.log("")
    
    # Load TAXA data
    logger.log("Loading TAXA table...")
    taxa_by_org_key, taxa_by_tvk = load_taxa_with_clean_fields(TAXA_FILE)
    logger.log(f"  Loaded {len(taxa_by_org_key)} taxa records")
    
    # Build index of child taxa by PARENT_TVK (for subspecies/varieties/forms)
    # Note: We include redundant taxa here because even if an infraspecific taxon
    # is marked redundant, it should still appear as a synonym of the valid parent species
    logger.log("Building child taxa index by PARENT_TVK...")
    children_by_parent_tvk = defaultdict(list)
    for org_key, taxon in taxa_by_org_key.items():
        parent_tvk = taxon.get('PARENT_TVK', '')
        if parent_tvk:
            children_by_parent_tvk[parent_tvk].append({
                'name': taxon['TAXON_NAME'],
                'tvk': taxon['TAXON_VERSION_KEY'],
                'rank': taxon['RANK'],
            })
    logger.log(f"  Built child index for {len(children_by_parent_tvk)} parent taxa")
    
    # First pass: identify all valid species (for conflict checking)
    # Note: We do NOT filter on REDUNDANT_FLAG here because species marked redundant
    # in TAXA may still be the recommended name that synonyms point to in NAMES.
    # Analysis showed 4,669 species are redundant in TAXA but valid in NAMES.
    logger.log("Identifying valid species for conflict checking...")
    valid_species_tvks = {}
    for org_key, taxon in taxa_by_org_key.items():
        if taxon['RANK'] in SPECIES_RANKS:
            valid_species_tvks[taxon['TAXON_VERSION_KEY']] = taxon
    logger.log(f"  Found {len(valid_species_tvks)} valid species in TAXA")
    
    # Load NAMES and build synonym index
    logger.log("Loading NAMES table for synonyms...")
    synonyms_by_rec_tvk = load_synonyms_index(NAMES_FILE, valid_species_tvks, logger)
    logger.log(f"  Built synonym index for {len(synonyms_by_rec_tvk)} valid taxa")
    logger.log("")
    
    # Output columns
    output_columns = [
        'organism_key',
        'taxon_version_key',
        'kingdom',
        'phylum_division',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'synonyms',
        'synonym_tvk_list',
        'recommended_name_authority',
        'non_native_flag',
        'terrestrial_freshwater_flag',
        'freshwater',
        'marine_flag',
    ]
    
    logger.log(f"Processing species...")
    
    species_count = 0
    invalid_count = 0
    skipped_rank = 0
    skipped_kingdom = 0
    excluded_identical_count = 0
    identical_tvks_captured = 0  # TVKs from identical-name synonyms
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8', newline='') as outfile, \
         open(INVALID_OUTPUT_FILE, 'w', encoding='utf-8', newline='') as invalidfile:
        
        writer = csv.DictWriter(outfile, fieldnames=output_columns, delimiter='\t')
        writer.writeheader()
        
        invalid_writer = csv.DictWriter(invalidfile, fieldnames=output_columns, delimiter='\t')
        invalid_writer.writeheader()
        
        for org_key, taxon in taxa_by_org_key.items():
            # Filter: species-level rank only (no REDUNDANT_FLAG filter - see README)
            if taxon['RANK'] not in SPECIES_RANKS:
                skipped_rank += 1
                continue
            
            # Get higher taxonomy
            hierarchy = get_higher_taxonomy(taxon, taxa_by_org_key)
            
            # Filter: target kingdoms only
            if hierarchy['Kingdom'] not in TARGET_KINGDOMS:
                skipped_kingdom += 1
                continue
            
            # Get valid species name
            species_name = taxon['TAXON_NAME']
            tvk = taxon['TAXON_VERSION_KEY']
            
            # Check if name is invalid/indeterminate
            is_invalid, reason = is_invalid_species_name(species_name)
            
            # Get synonyms from NAMES table
            synonym_entries = synonyms_by_rec_tvk.get(tvk, [])
            
            # Process synonyms: 
            # - Names different from valid name go to synonyms list
            # - Names identical to valid name: capture their TVKs separately
            # Store as list of (name, tvk) tuples
            filtered_synonyms = []
            identical_name_tvks = []  # TVKs for names identical to valid species name
            
            for syn_entry in synonym_entries:
                syn_name = syn_entry['name']
                syn_tvk = syn_entry['tvk']
                if syn_name == species_name:
                    # Name is identical - don't add to synonyms list but capture the TVK
                    logger.log_excluded_identical_synonym(species_name, tvk, syn_name, syn_tvk)
                    excluded_identical_count += 1
                    if syn_tvk:  # Only add if TVK exists
                        identical_name_tvks.append(syn_tvk)
                else:
                    filtered_synonyms.append((syn_name, syn_tvk))
            
            # Add child taxa from TAXA hierarchy (subspecies, varieties, forms, etc.)
            # These are linked via PARENT_TVK in the TAXA table
            child_taxa = children_by_parent_tvk.get(tvk, [])
            for child in child_taxa:
                child_name = child['name']
                child_tvk = child['tvk']
                # Don't add if name is identical to species or already in list
                existing_names = [s[0] for s in filtered_synonyms]
                if child_name != species_name and child_name not in existing_names:
                    filtered_synonyms.append((child_name, child_tvk))
            
            # Add subgenus-derived synonyms if applicable (these have no TVK)
            subgenus_synonyms = extract_subgenus_synonyms(species_name)
            for sub_syn in subgenus_synonyms:
                # Check if this name is already in the list
                existing_names = [s[0] for s in filtered_synonyms]
                if sub_syn not in existing_names and sub_syn != species_name:
                    filtered_synonyms.append((sub_syn, ''))  # No TVK for derived synonyms
            
            # Remove duplicates by name, keeping first occurrence (which has TVK if available)
            seen_names = set()
            unique_synonyms = []
            for syn_name, syn_tvk in filtered_synonyms:
                if syn_name not in seen_names:
                    seen_names.add(syn_name)
                    unique_synonyms.append((syn_name, syn_tvk))
            
            # Sort by name
            unique_synonyms.sort(key=lambda x: x[0])
            
            # Format synonyms: semicolon-delimited, no spaces around semicolon
            synonyms_str = ';'.join(s[0] for s in unique_synonyms)
            
            # Build TVK list: include synonym TVKs + identical-name TVKs
            # The identical-name TVKs are added at the end (they have no corresponding synonym name)
            all_tvks = [s[1] for s in unique_synonyms] + identical_name_tvks
            synonym_tvks_str = ';'.join(all_tvks)
            identical_tvks_captured += len(identical_name_tvks)
            
            # Combine Phylum/Division
            phylum_division = hierarchy['Phylum'] if hierarchy['Phylum'] else hierarchy['Division']
            
            # Build output row
            output_row = {
                'organism_key': taxon['ORGANISM_KEY'],
                'taxon_version_key': tvk,
                'kingdom': hierarchy['Kingdom'],
                'phylum_division': phylum_division,
                'class': hierarchy['Class'],
                'order': hierarchy['Order'],
                'family': hierarchy['Family'],
                'genus': hierarchy['Genus'],
                'species': species_name,
                'synonyms': synonyms_str,
                'synonym_tvk_list': synonym_tvks_str,
                'recommended_name_authority': taxon['TAXON_AUTHORITY'],
                'non_native_flag': taxon['NON_NATIVE_FLAG'],
                'terrestrial_freshwater_flag': taxon['TERRESTRIAL_FRESHWATER_FLAG'],
                'freshwater': taxon['FRESHWATER'],
                'marine_flag': taxon['MARINE_FLAG'],
            }
            
            if is_invalid:
                invalid_writer.writerow(output_row)
                invalid_count += 1
                logger.log_invalid_species(species_name, tvk, reason)
            else:
                writer.writerow(output_row)
                species_count += 1
            
            if (species_count + invalid_count) % 10000 == 0:
                logger.log(f"  Processed {species_count + invalid_count} species...")
    
    logger.log("")
    logger.log("=" * 60)
    logger.log("PROCESSING COMPLETE")
    logger.log("=" * 60)
    logger.log(f"Valid species written: {species_count}")
    logger.log(f"Invalid species written: {invalid_count}")
    logger.log(f"Skipped (not species rank): {skipped_rank}")
    logger.log(f"Skipped (not target kingdom): {skipped_kingdom}")
    logger.log(f"Identical-name synonyms (name excluded, TVK captured): {excluded_identical_count}")
    logger.log(f"Identical-name TVKs added to synonym_tvk_list: {identical_tvks_captured}")
    logger.log(f"Valid-as-synonym conflicts: {len(logger.valid_as_synonym_conflicts)}")
    logger.log(f"Output valid file: {OUTPUT_FILE}")
    logger.log(f"Output invalid file: {INVALID_OUTPUT_FILE}")
    logger.log(f"Log file: {LOG_FILE}")
    
    # Write log file
    logger.write_log()
    print(f"\nLog written to: {LOG_FILE}")


if __name__ == '__main__':
    process_species()
