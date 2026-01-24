"""
UKSI Species Extraction Pipeline v2.2
======================================
Extracts valid species names from UKSI TAXA and NAMES tables,
with higher taxonomy and synonyms.

Author: Claude (for Ben Price, NHM)
Date: 2026-01-24
Version: 2.2

Changes in v2.2:
- FIXED DEDUPLICATION: Now allows duplicate synonym names if they have different TVKs.
  (Fixes issue where deprecated TVKs were discarded because a valid TVK for the same name existed).
- Relaxed synonym RANK filtering (v2.1 change retained).
- Debug tracing enabled for specific TVKs.
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

# Debugging: Add TVKs here to see exactly why they are accepted/rejected in the log
DEBUG_TVKS = ['BMSSYS0000013986'] 

# Kingdoms to include
TARGET_KINGDOMS = {'Animalia', 'Plantae', 'Fungi', 'Chromista'}

# Ranks that count as "species level" for the main output
SPECIES_RANK = 'Species'

# Ranks to include in synonyms
SYNONYM_RANKS = {
    'Species', 'Species aggregate', 'Species group', 'Species hybrid',
    'Species pro parte', 'Species sensu lato', 'Subspecies',
    'Subspecies aggregate', 'Subspecies hybrid', 'Variety', 'Varietal hybrid',
    'Subvariety', 'Form', 'Subform', 'Nothosubspecies', 'Nothovariety',
    'Microspecies', 'Praespecies', 'Cultivar', 'Convariety', 'Abberation',
    'Breed', 'Facies', 'Morphotype', 'Race', 'Pathovar', 'Forma specialis', 'ecad',
}

HIGHER_RANKS = ['Kingdom', 'Phylum', 'Division', 'Class', 'Order', 'Family', 'Genus']

INVALID_NAME_PATTERNS = [
    r'\.', r'\?', r'\(Other\)', r'"', r'\(unidentified\)', r'\bindet\b'
]
INVALID_NAME_REGEX = re.compile('|'.join(INVALID_NAME_PATTERNS), re.IGNORECASE)
SUBGENUS_PATTERN = re.compile(r'^(\w+)\s+\((\w+)\)\s+(.+)$')


class PipelineLogger:
    def __init__(self, log_path):
        self.log_path = log_path
        self.log_entries = []
        self.excluded_identical_synonyms = []
        self.valid_as_synonym_conflicts = []
        self.invalid_species = []
        self.debug_entries = []
        
    def log(self, message):
        self.log_entries.append(message)
        print(message)

    def log_debug(self, tvk, message):
        entry = f"[DEBUG TRACE] {tvk}: {message}"
        self.debug_entries.append(entry)
        print(entry)
    
    def log_excluded_identical_synonym(self, species_name, tvk, synonym_name, synonym_tvk):
        self.excluded_identical_synonyms.append({
            'species': species_name, 'species_tvk': tvk,
            'synonym': synonym_name, 'synonym_tvk': synonym_tvk,
        })
    
    def log_valid_as_synonym_conflict(self, valid_name, valid_tvk, synonym_tvk, rec_tvk):
        self.valid_as_synonym_conflicts.append({
            'valid_name': valid_name, 'valid_tvk': valid_tvk,
            'synonym_tvk': synonym_tvk, 'recommended_tvk': rec_tvk,
        })
    
    def log_invalid_species(self, species_name, tvk, reason):
        self.invalid_species.append({'species': species_name, 'tvk': tvk, 'reason': reason})
    
    def write_log(self):
        with open(self.log_path, 'w', encoding='utf-8') as f:
            f.write(f"UKSI Pipeline Log - {datetime.now().isoformat()}\n")
            f.write("=" * 80 + "\n\n")
            
            if self.debug_entries:
                f.write("DEBUG TRACES\n" + "-" * 40 + "\n")
                for entry in self.debug_entries: f.write(entry + "\n")
                f.write("\n")

            f.write("PROCESSING LOG\n" + "-" * 40 + "\n")
            for entry in self.log_entries: f.write(entry + "\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write(f"EXCLUDED IDENTICAL SYNONYMS ({len(self.excluded_identical_synonyms)})\n")
            for entry in self.excluded_identical_synonyms[:50]:
                f.write(f"{entry['species']} <- {entry['synonym']} ({entry['synonym_tvk']})\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write(f"INVALID SPECIES ({len(self.invalid_species)})\n")
            for entry in self.invalid_species[:50]:
                f.write(f"{entry['species']} ({entry['reason']})\n")

def clean_field(value):
    if value:
        value = value.replace('\n', ' ').replace('\r', ' ')
        value = ' '.join(value.split())
        return value.strip()
    return value

def load_taxa_with_clean_fields(filepath):
    taxa_by_org_key = {}
    taxa_by_tvk = {}
    with open(filepath, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cleaned_row = {k: clean_field(v) for k, v in row.items()}
            taxa_by_org_key[cleaned_row['ORGANISM_KEY']] = cleaned_row
            taxa_by_tvk[cleaned_row['TAXON_VERSION_KEY']] = cleaned_row
    return taxa_by_org_key, taxa_by_tvk

def load_synonyms_index(filepath, valid_species_tvks, logger):
    synonyms_by_rec_tvk = defaultdict(list)
    valid_tvk_set = set(valid_species_tvks.keys())
    
    logger.log("  Scanning names for valid-as-synonym conflicts...")
    
    with open(filepath, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tvk = clean_field(row.get('TAXON_VERSION_KEY', ''))
            rec_tvk = clean_field(row.get('RECOMMENDED_TAXON_VERSION_KEY', ''))
            
            if tvk in DEBUG_TVKS:
                logger.log_debug(tvk, f"Found in NAMES. Rank: '{row.get('RANK','')}', Rec_TVK: '{rec_tvk}'")

            if tvk in valid_tvk_set and rec_tvk and tvk != rec_tvk:
                valid_name = valid_species_tvks[tvk]['TAXON_NAME']
                logger.log_valid_as_synonym_conflict(valid_name, tvk, tvk, rec_tvk)

            if clean_field(row.get('LANGUAGE', '')) != 'la': continue
            if tvk == rec_tvk: continue
            if tvk in valid_tvk_set: continue

            points_to_valid_species = (rec_tvk in valid_tvk_set)
            rank = clean_field(row.get('RANK', ''))
            
            if not points_to_valid_species and rank not in SYNONYM_RANKS:
                if tvk in DEBUG_TVKS: logger.log_debug(tvk, f"REJECTED: Rank '{rank}' not in permitted list")
                continue

            if rec_tvk:
                if tvk in DEBUG_TVKS: logger.log_debug(tvk, f"ACCEPTED: Added as synonym for {rec_tvk}")
                synonyms_by_rec_tvk[rec_tvk].append({
                    'name': clean_field(row.get('TAXON_NAME', '')),
                    'tvk': tvk,
                })
    return synonyms_by_rec_tvk

def get_higher_taxonomy(species_row, taxa_by_org_key):
    hierarchy = {rank: '' for rank in HIGHER_RANKS}
    current = species_row
    visited = set()
    while current and current['ORGANISM_KEY'] not in visited:
        visited.add(current['ORGANISM_KEY'])
        if current['RANK'] in hierarchy: hierarchy[current['RANK']] = current['TAXON_NAME']
        parent_key = current.get('PARENT_KEY', '')
        if parent_key in taxa_by_org_key: current = taxa_by_org_key[parent_key]
        else: break
    return hierarchy

def is_invalid_species_name(name):
    match = INVALID_NAME_REGEX.search(name)
    if match: return True, f"Contains '{match.group()}'"
    return False, None

def extract_subgenus_synonyms(species_name):
    match = SUBGENUS_PATTERN.match(species_name)
    if match:
        return [f"{match.group(2)} {match.group(3)}", f"{match.group(1)} {match.group(3)}"]
    return []

def process_species():
    logger = PipelineLogger(LOG_FILE)
    logger.log(f"Pipeline started: {datetime.now().isoformat()} (v2.2)")
    
    logger.log("Loading TAXA table...")
    taxa_by_org_key, taxa_by_tvk = load_taxa_with_clean_fields(TAXA_FILE)
    
    logger.log("Building child taxa index...")
    children_by_parent_tvk = defaultdict(list)
    for org_key, taxon in taxa_by_org_key.items():
        if taxon.get('PARENT_TVK'):
            children_by_parent_tvk[taxon['PARENT_TVK']].append({
                'name': taxon['TAXON_NAME'], 'tvk': taxon['TAXON_VERSION_KEY']
            })
            
    valid_species_tvks = {t['TAXON_VERSION_KEY']: t for t in taxa_by_org_key.values() if t['RANK'] == SPECIES_RANK}
    logger.log(f"  Found {len(valid_species_tvks)} valid species")
    
    logger.log("Loading NAMES table...")
    synonyms_by_rec_tvk = load_synonyms_index(NAMES_FILE, valid_species_tvks, logger)
    
    output_cols = ['organism_key', 'taxon_version_key', 'kingdom', 'phylum_division', 
                   'class', 'order', 'family', 'genus', 'species', 'synonyms', 'synonym_tvk_list',
                   'recommended_name_authority', 'non_native_flag', 'terrestrial_freshwater_flag',
                   'freshwater', 'marine_flag']
    
    species_count = 0
    with open(OUTPUT_FILE, 'w', encoding='utf-8', newline='') as outfile, \
         open(INVALID_OUTPUT_FILE, 'w', encoding='utf-8', newline='') as invalidfile:
        
        writer = csv.DictWriter(outfile, fieldnames=output_cols, delimiter='\t')
        writer.writeheader()
        invalid_writer = csv.DictWriter(invalidfile, fieldnames=output_cols, delimiter='\t')
        invalid_writer.writeheader()
        
        for org_key, taxon in taxa_by_org_key.items():
            if taxon['RANK'] != SPECIES_RANK: continue
            
            hierarchy = get_higher_taxonomy(taxon, taxa_by_org_key)
            if hierarchy['Kingdom'] not in TARGET_KINGDOMS: continue
            
            species_name = taxon['TAXON_NAME']
            tvk = taxon['TAXON_VERSION_KEY']
            is_invalid, reason = is_invalid_species_name(species_name)
            
            synonym_entries = synonyms_by_rec_tvk.get(tvk, [])
            
            # 1. Gather all potential synonyms
            filtered_synonyms = []
            for syn_entry in synonym_entries:
                if syn_entry['name'] == species_name:
                    logger.log_excluded_identical_synonym(species_name, tvk, syn_entry['name'], syn_entry['tvk'])
                else:
                    filtered_synonyms.append((syn_entry['name'], syn_entry['tvk']))
            
            # 2. Add child taxa
            for child in children_by_parent_tvk.get(tvk, []):
                if child['name'] != species_name:
                    filtered_synonyms.append((child['name'], child['tvk']))
            
            # 3. Add subgenus forms (no TVK)
            for sub_syn in extract_subgenus_synonyms(species_name):
                if sub_syn != species_name:
                    filtered_synonyms.append((sub_syn, ''))
            
            # 4. Deduplicate by (NAME + TVK) pair to preserve multiple TVKs for same name
            seen_entries = set()
            unique_synonyms = []
            for syn_name, syn_tvk in filtered_synonyms:
                entry_key = (syn_name, syn_tvk)
                if entry_key not in seen_entries:
                    seen_entries.add(entry_key)
                    unique_synonyms.append((syn_name, syn_tvk))
            
            unique_synonyms.sort(key=lambda x: x[0])
            
            row_data = {
                'organism_key': taxon['ORGANISM_KEY'],
                'taxon_version_key': tvk,
                'kingdom': hierarchy['Kingdom'],
                'phylum_division': hierarchy.get('Phylum') or hierarchy.get('Division'),
                'class': hierarchy['Class'],
                'order': hierarchy['Order'],
                'family': hierarchy['Family'],
                'genus': hierarchy['Genus'],
                'species': species_name,
                'synonyms': ';'.join(s[0] for s in unique_synonyms),
                'synonym_tvk_list': ';'.join(s[1] for s in unique_synonyms),
                'recommended_name_authority': taxon['TAXON_AUTHORITY'],
                'non_native_flag': taxon['NON_NATIVE_FLAG'],
                'terrestrial_freshwater_flag': taxon['TERRESTRIAL_FRESHWATER_FLAG'],
                'freshwater': taxon['FRESHWATER'],
                'marine_flag': taxon['MARINE_FLAG'],
            }
            
            if is_invalid:
                invalid_writer.writerow(row_data)
                logger.log_invalid_species(species_name, tvk, reason)
            else:
                writer.writerow(row_data)
                species_count += 1
            
            if species_count % 10000 == 0: logger.log(f"  Processed {species_count} species...")

    logger.log("PROCESSING COMPLETE")
    logger.write_log()

if __name__ == '__main__':
    process_species()