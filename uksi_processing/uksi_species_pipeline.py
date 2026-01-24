"""
UKSI Species Extraction Pipeline
================================
Extracts valid species names from UKSI TAXA and NAMES tables,
with higher taxonomy and synonyms.

Author: Claude (for Ben Price, NHM)
Date: 2026-01-24
"""

import csv
import sys
from collections import defaultdict

# Configuration
TAXA_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_20251203a_input_taxa.tsv'
NAMES_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_20251203a_input_names.tsv'
OUTPUT_FILE = r'C:\_claude_files\projects\ukbol_gaplist\uksi\uksi_valid_species_output.tsv'

# Kingdoms to include
TARGET_KINGDOMS = {'Animalia', 'Plantae', 'Fungi', 'Chromista'}

# Ranks that count as "species level" for the main output
SPECIES_RANK = 'Species'

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


def load_taxa_index():
    """Load TAXA table and index by ORGANISM_KEY for hierarchy traversal."""
    print("Loading TAXA table...")
    taxa_by_org_key = {}
    taxa_by_tvk = {}
    
    with open(TAXA_FILE, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            taxa_by_org_key[row['ORGANISM_KEY']] = row
            taxa_by_tvk[row['TAXON_VERSION_KEY']] = row
    
    print(f"  Loaded {len(taxa_by_org_key)} taxa records")
    return taxa_by_org_key, taxa_by_tvk


def load_synonyms_index():
    """
    Load NAMES table and build index of synonyms by RECOMMENDED_TAXON_VERSION_KEY.
    
    Includes all names where:
    - DEPRECATED_DATE is empty (not deprecated)
    - LANGUAGE = 'la' (scientific names only)
    - TAXON_VERSION_KEY != RECOMMENDED_TAXON_VERSION_KEY (excludes the valid name itself)
    - RANK is species-level or infraspecific
    """
    print("Loading NAMES table for synonyms...")
    synonyms_by_rec_tvk = defaultdict(list)
    
    with open(NAMES_FILE, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Filter criteria
            if row['DEPRECATED_DATE'] != '':
                continue  # Skip deprecated names
            if row['LANGUAGE'] != 'la':
                continue  # Skip non-scientific names
            if row['RANK'] not in SYNONYM_RANKS:
                continue  # Skip non-species/infraspecific ranks
            if row['TAXON_VERSION_KEY'] == row['RECOMMENDED_TAXON_VERSION_KEY']:
                continue  # Skip the valid name itself
            
            rec_tvk = row['RECOMMENDED_TAXON_VERSION_KEY']
            if rec_tvk:
                synonyms_by_rec_tvk[rec_tvk].append(row['TAXON_NAME'])
    
    print(f"  Built synonym index for {len(synonyms_by_rec_tvk)} valid taxa")
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


def process_species():
    """Main processing function."""
    
    # Load data
    taxa_by_org_key, taxa_by_tvk = load_taxa_index()
    synonyms_by_rec_tvk = load_synonyms_index()
    
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
        'recommended_name_authority',
        'non_native_flag',
        'terrestrial_freshwater_flag',
        'freshwater',
        'marine_flag',
    ]
    
    print(f"Processing species and writing to {OUTPUT_FILE}...")
    
    species_count = 0
    skipped_redundant = 0
    skipped_rank = 0
    skipped_kingdom = 0
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_columns, delimiter='\t')
        writer.writeheader()
        
        for org_key, taxon in taxa_by_org_key.items():
            # Filter: not redundant
            if taxon['REDUNDANT_FLAG'] != '':
                skipped_redundant += 1
                continue
            
            # Filter: species rank only
            if taxon['RANK'] != SPECIES_RANK:
                skipped_rank += 1
                continue
            
            # Get higher taxonomy
            hierarchy = get_higher_taxonomy(taxon, taxa_by_org_key)
            
            # Filter: target kingdoms only
            if hierarchy['Kingdom'] not in TARGET_KINGDOMS:
                skipped_kingdom += 1
                continue
            
            # Get synonyms
            tvk = taxon['TAXON_VERSION_KEY']
            synonym_list = synonyms_by_rec_tvk.get(tvk, [])
            # Remove duplicates and sort
            synonym_list = sorted(set(synonym_list))
            synonyms_str = '; '.join(synonym_list)
            
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
                'species': taxon['TAXON_NAME'],
                'synonyms': synonyms_str,
                'recommended_name_authority': taxon['TAXON_AUTHORITY'],
                'non_native_flag': taxon['NON_NATIVE_FLAG'],
                'terrestrial_freshwater_flag': taxon['TERRESTRIAL_FRESHWATER_FLAG'],
                'freshwater': taxon['FRESHWATER'],
                'marine_flag': taxon['MARINE_FLAG'],
            }
            
            writer.writerow(output_row)
            species_count += 1
            
            if species_count % 10000 == 0:
                print(f"  Processed {species_count} species...")
    
    print()
    print("=" * 60)
    print("PROCESSING COMPLETE")
    print("=" * 60)
    print(f"Valid species written: {species_count}")
    print(f"Skipped (redundant): {skipped_redundant}")
    print(f"Skipped (not species rank): {skipped_rank}")
    print(f"Skipped (not target kingdom): {skipped_kingdom}")
    print(f"Output file: {OUTPUT_FILE}")


if __name__ == '__main__':
    process_species()
