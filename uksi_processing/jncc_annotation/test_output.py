import csv

with open(r'C:\GitHub\mind-the-gap\uksi_processing\jncc_annotation\uksi_valid_species_jncc_annotated.tsv', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    count = 0
    for row in reader:
        if row.get('jncc_matching_tvk') and any(row.get(c) for c in ['A: Bern Convention', 'I: Wildlife and Countryside Act 1981', 'Ha: Biodiversity Action Plan UK list of priority species']):
            print(f"Species: {row['species']}")
            print(f"TVK: {row['taxon_version_key']}")
            print(f"Bern: {row.get('A: Bern Convention', '')}")
            print(f"WCA: {row.get('I: Wildlife and Countryside Act 1981', '')}")
            print(f"BAP: {row.get('Ha: Biodiversity Action Plan UK list of priority species', '')}")
            print('---')
            count += 1
            if count >= 5:
                break
