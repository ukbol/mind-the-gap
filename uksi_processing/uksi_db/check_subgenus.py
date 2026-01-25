#!/usr/bin/env python3
"""Check subgenus synonym generation."""
import csv

with open(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi_species_export.tsv', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    count = 0
    for row in reader:
        name = row['TAXON_NAME']
        if '(' in name and ')' in name and count < 10:
            print(f"Name: {name}")
            print(f"  Synonyms: {row['synonyms']}")
            print()
            count += 1
