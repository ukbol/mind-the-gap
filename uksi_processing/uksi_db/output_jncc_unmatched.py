#!/usr/bin/env python3
"""
Output unmatched JNCC records with their recommended names for reference.
"""
import sqlite3
import csv

conn = sqlite3.connect(r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\uksi.db')
cur = conn.cursor()

# Get all unmatched JNCC records with their resolution path
cur.execute("""
    SELECT 
        j.informal_group,
        j.Recommended_taxon_name as jncc_name,
        j.Recommended_taxon_version as jncc_tvk,
        n.RECOMMENDED_SCIENTIFIC_NAME as resolved_name,
        n.RECOMMENDED_TAXON_VERSION_KEY as resolved_tvk,
        CASE WHEN EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = n.RECOMMENDED_TAXON_VERSION_KEY) 
             THEN 'YES' ELSE 'NO' END as resolved_in_taxa
    FROM jncc j
    LEFT JOIN names n ON n.TAXON_VERSION_KEY = j.Recommended_taxon_version
    WHERE NOT EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = j.Recommended_taxon_version)
    ORDER BY j.informal_group, j.Recommended_taxon_name
""")

results = cur.fetchall()

# Write to TSV
output_path = r'C:\GitHub\mind-the-gap\uksi_processing\uksi_db\jncc_unmatched_analysis.tsv'
with open(output_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['informal_group', 'jncc_name', 'jncc_tvk', 'resolved_name', 'resolved_tvk', 'resolved_in_taxa'])
    writer.writerows(results)

print(f"Written {len(results)} unmatched records to: {output_path}")

# Summary stats
resolved_yes = sum(1 for r in results if r[5] == 'YES')
resolved_no = sum(1 for r in results if r[5] == 'NO')
no_names_match = sum(1 for r in results if r[3] is None)

print(f"\nSummary:")
print(f"  - Total unmatched: {len(results)}")
print(f"  - Resolvable via names table: {resolved_yes}")
print(f"  - NOT resolvable (no taxa match): {resolved_no}")
print(f"  - No match in names table: {no_names_match}")

conn.close()
