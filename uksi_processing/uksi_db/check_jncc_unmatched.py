#!/usr/bin/env python3
"""Investigate unmatched JNCC records."""
import sqlite3

conn = sqlite3.connect(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi.db')
cur = conn.cursor()

# Find JNCC records that don't match taxa directly
cur.execute("""
    SELECT 
        j.informal_group,
        j.Recommended_taxon_name,
        j.Recommended_taxon_version,
        CASE WHEN EXISTS (SELECT 1 FROM names n WHERE n.TAXON_VERSION_KEY = j.Recommended_taxon_version) 
             THEN 'YES' ELSE 'NO' END as in_names,
        CASE WHEN EXISTS (SELECT 1 FROM names n WHERE n.RECOMMENDED_TAXON_VERSION_KEY = j.Recommended_taxon_version) 
             THEN 'YES' ELSE 'NO' END as is_recommended_in_names
    FROM jncc j
    WHERE NOT EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = j.Recommended_taxon_version)
    ORDER BY j.informal_group, j.Recommended_taxon_name
""")

results = cur.fetchall()
print(f"Found {len(results)} unmatched JNCC records\n")

# Count by match status
in_names_count = sum(1 for r in results if r[3] == 'YES')
is_recommended_count = sum(1 for r in results if r[4] == 'YES')
print(f"Of these:")
print(f"  - {in_names_count} have TVK in names table")
print(f"  - {is_recommended_count} have TVK as a RECOMMENDED_TAXON_VERSION_KEY in names table")
print()

# Show sample of unmatched
print("Sample unmatched records:")
print("-" * 100)
for row in results[:50]:
    print(f"  {row[0]:30} | {row[1]:40} | {row[2]} | in_names={row[3]} | is_rec={row[4]}")

# Check if the TVKs exist in names but point to different recommended TVKs
print("\n\nChecking if unmatched TVKs have different recommended names in names table:")
print("-" * 100)
cur.execute("""
    SELECT 
        j.Recommended_taxon_name as jncc_name,
        j.Recommended_taxon_version as jncc_tvk,
        n.TAXON_NAME as names_taxon,
        n.RECOMMENDED_TAXON_VERSION_KEY as names_rec_tvk,
        n.RECOMMENDED_SCIENTIFIC_NAME as names_rec_name,
        t.TAXON_NAME as taxa_name
    FROM jncc j
    JOIN names n ON n.TAXON_VERSION_KEY = j.Recommended_taxon_version
    LEFT JOIN taxa t ON t.TAXON_VERSION_KEY = n.RECOMMENDED_TAXON_VERSION_KEY
    WHERE NOT EXISTS (SELECT 1 FROM taxa t2 WHERE t2.TAXON_VERSION_KEY = j.Recommended_taxon_version)
    LIMIT 30
""")

for row in cur.fetchall():
    print(f"  JNCC: {row[0][:35]:35} -> Names rec: {row[4][:35] if row[4] else 'NULL':35} -> Taxa: {row[5] if row[5] else 'NOT IN TAXA'}")

conn.close()
