#!/usr/bin/env python3
"""Verify JNCC propagation is working correctly."""
import sqlite3
import csv

conn = sqlite3.connect(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi.db')
cur = conn.cursor()

# Check what JNCC designations Cetacea has
print("=== JNCC designations for Cetacea (Infraorder) ===")
cur.execute("""
    SELECT j.* 
    FROM jncc_resolved j
    JOIN taxa t ON t.TAXON_VERSION_KEY = j.resolved_tvk
    WHERE t.TAXON_NAME = 'Cetacea'
""")
cols = [desc[0] for desc in cur.description]
for row in cur.fetchall():
    for i, val in enumerate(row):
        if val and val.strip() and cols[i] not in ('informal_group', 'Recommended_taxon_name', 'Recommended_authority', 'Recommended_qualifier', 'Recommended_taxon_version', 'resolved_tvk'):
            print(f"  {cols[i]}: {val}")

# Now check if whale species inherited these
print("\n=== Checking export file for whale species ===")
with open(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi_species_export.tsv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if 'Balaenoptera' in row.get('TAXON_NAME', ''):
            print(f"\n{row['TAXON_NAME']}:")
            # Show JNCC columns that have values
            for key, val in row.items():
                if key.startswith('jncc_') and val:
                    print(f"  {key}: {val}")

conn.close()
