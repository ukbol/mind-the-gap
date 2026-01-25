import csv

print("=== Species matched via synonym TVK only ===")
with open(r'C:\GitHub\mind-the-gap\uksi_processing\jncc_annotation\uksi_valid_species_jncc_annotated.tsv', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    count = 0
    for row in reader:
        if row.get('tvk_match_status') == 'synonym':
            print(f"Species: {row['species']}")
            print(f"Valid TVK: {row['taxon_version_key']}")
            print(f"Synonym TVKs: {row.get('synonym_tvk_list', '')[:100]}...")
            print(f"Matching TVK: {row['jncc_matching_tvk']}")
            print(f"Match status: {row['tvk_match_status']}")
            print('---')
            count += 1
            if count >= 5:
                break

print("\n=== Species matched via both valid and synonym TVK ===")
with open(r'C:\GitHub\mind-the-gap\uksi_processing\jncc_annotation\uksi_valid_species_jncc_annotated.tsv', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    count = 0
    for row in reader:
        if row.get('tvk_match_status') == 'valid;synonym':
            print(f"Species: {row['species']}")
            print(f"Valid TVK: {row['taxon_version_key']}")
            print(f"Matching TVKs: {row['jncc_matching_tvk']}")
            print(f"Match status: {row['tvk_match_status']}")
            print(f"BAP: {row.get('Ha: Biodiversity Action Plan UK list of priority species', '')}")
            print('---')
            count += 1
            if count >= 5:
                break
