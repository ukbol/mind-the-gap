#!/usr/bin/env python3
"""Quick test of the UKSI database."""
import sqlite3

conn = sqlite3.connect(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi.db')
cur = conn.cursor()

# Test: Get a sample species and its synonyms
cur.execute("""
    SELECT t.TAXON_NAME, t.RANK, n.TAXON_NAME, n.NAME_STATUS, n.NAME_FORM
    FROM taxa t
    JOIN names n ON n.RECOMMENDED_TAXON_VERSION_KEY = t.TAXON_VERSION_KEY
    WHERE t.TAXON_NAME = 'Bombus terrestris'
    LIMIT 10
""")
print('Bombus terrestris and its names:')
for row in cur.fetchall():
    print(f'  {row}')

# Test: Check JNCC linkage
cur.execute("""
    SELECT t.TAXON_NAME, j."Ha: Biodiversity Action Plan UK list of priority species"
    FROM taxa t
    JOIN jncc j ON j.Recommended_taxon_version = t.TAXON_VERSION_KEY
    WHERE j."Ha: Biodiversity Action Plan UK list of priority species" IS NOT NULL 
    AND j."Ha: Biodiversity Action Plan UK list of priority species" != ''
    LIMIT 5
""")
print('\nSample JNCC BAP species:')
for row in cur.fetchall():
    print(f'  {row}')

# Test: Check Pantheon linkage
cur.execute("""
    SELECT t.TAXON_NAME, p.SQS, p."Conservation status", p."Broad biotope"
    FROM taxa t
    JOIN pantheon p ON p.RECOMMENDED_TAXON_VERSION_KEY = t.TAXON_VERSION_KEY
    WHERE p.SQS IS NOT NULL AND p.SQS != ''
    LIMIT 5
""")
print('\nSample Pantheon data:')
for row in cur.fetchall():
    print(f'  {row}')

conn.close()
