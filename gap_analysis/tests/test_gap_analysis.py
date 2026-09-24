"""Tests for gap_analysis.py status/issue assignment and the legacy converter."""

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pytest

GAP_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GAP_DIR))

import gap_analysis as ga  # noqa: E402
import convert_legacy_status as conv  # noqa: E402


def build_indexes(records):
    """Build analyze_taxon indexes from (species, cluster_or_None, count) tuples."""
    name_to_count = defaultdict(int)
    name_to_bins = defaultdict(set)
    bin_to_names = defaultdict(set)
    for species, cluster, count in records:
        name = ga.normalize_species_name(species)
        name_to_count[name] += count
        if cluster:
            name_to_bins[name].add(cluster)
            bin_to_names[cluster].add(name)
    return name_to_count, name_to_bins, bin_to_names, name_to_bins, {}


def analyze(valid, synonyms, records):
    taxon = ga.Taxon(row_index=0, valid_name=valid, synonyms=synonyms, input_data={})
    return ga.analyze_taxon(taxon, *build_indexes(records))


@pytest.mark.parametrize('records, status, grade, issues', [
    ([], 'no_records', 'F', set()),
    ([('Aus bus', 'B1', 12)], 'valid_name', 'A', set()),
    ([('Aus bus', 'B1', 5)], 'valid_name', 'B', set()),
    ([('Aus bus', 'B1', 2)], 'valid_name', 'D', {'few_records'}),
    ([('Aus bus', 'B1', 5), ('Aus bus', 'B2', 5)], 'valid_name', 'C', {'split_bins'}),
    ([('Aus bus', None, 4)], 'valid_name', 'F', {'no_cluster'}),
    ([('Aus bus', 'B1', 5), ('Aus old', 'B1', 5)], 'valid_and_synonym', 'B', {'synonym_records'}),
    ([('Aus old', 'B1', 5)], 'synonym_only', 'B', {'synonym_records', 'valid_name_absent'}),
    ([('Aus bus', 'B1', 5), ('Aus cus', 'B1', 1)], 'shared_bin_species', 'E', {'shared_bin_species'}),
    ([('Aus bus', 'B1', 5), ('Aus sp. X1', 'B1', 1)], 'shared_bin_interim', 'E', {'shared_bin_interim'}),
])
def test_status_grade_and_issues(records, status, grade, issues):
    result = analyze('Aus bus', ['Aus old'], records)
    assert result.species_status == status
    assert result.bags_grade == grade
    assert result.issues == issues


def test_shared_bin_keeps_all_other_issues():
    # Previously RED hid the synonym and split problems; now all are flagged
    result = analyze('Aus bus', ['Aus old'], [
        ('Aus old', 'B1', 1),
        ('Aus cus', 'B1', 3),
        ('Aus cf. bus', 'B2', 1),
        ('Aus old', 'B2', 1),
    ])
    assert result.species_status == 'shared_bin_species'
    assert result.bags_grade == 'E'
    assert result.issues == {
        'shared_bin_species', 'shared_bin_interim', 'synonym_records',
        'valid_name_absent', 'split_bins', 'few_records',
    }
    assert ga.format_issues(result.issues) == (
        'shared_bin_species;shared_bin_interim;synonym_records;'
        'valid_name_absent;split_bins;few_records'
    )


def test_legacy_map_covers_every_status():
    assert set(ga.LEGACY_STATUS_MAP.values()) == set(ga.STATUS_ORDER)


def write_tsv(path, rows):
    with open(path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def test_end_to_end_output_columns(tmp_path):
    species = tmp_path / 'species.tsv'
    records = tmp_path / 'records.tsv'
    output = tmp_path / 'out.tsv'
    write_tsv(species, [
        {'taxon_name': 'Aus bus', 'synonyms': 'Aus old'},
        {'taxon_name': 'Aus cus', 'synonyms': ''},
        {'taxon_name': 'Dus eus', 'synonyms': ''},
    ])
    write_tsv(records, [
        {'species': 'Aus bus', 'bin_uri': 'BOLD:1'},
        {'species': 'Aus old', 'bin_uri': 'BOLD:2'},
        {'species': 'Aus cus', 'bin_uri': 'BOLD:3'},
    ])
    subprocess.run(
        [sys.executable, str(GAP_DIR / 'gap_analysis.py'), '--species-list', str(species),
         '--records', str(records), '--output', str(output), '--workers', '1'],
        check=True, capture_output=True,
    )
    rows = read_tsv(output)
    header = list(rows[0])
    assert header.index('issues') == header.index('species_status') + 1
    assert [(r['species_status'], r['issues']) for r in rows] == [
        ('valid_and_synonym', 'synonym_records;split_bins;few_records'),
        ('valid_name', 'few_records'),
        ('no_records', ''),
    ]


def test_convert_gene_file(tmp_path):
    path = tmp_path / 'gene.tsv'
    base = {'taxon_name': '', 'number_records': '0', 'bags_grade': 'F',
            'species_status': 'BLACK', 'bin_uris': '', 'otu_ids': '', 'other_names': ''}
    write_tsv(path, [
        dict(base),
        dict(base, number_records='20', bags_grade='A', species_status='GREEN'),
        dict(base, number_records='2', bags_grade='D', species_status='AMBER'),
        dict(base, number_records='5', bags_grade='C', species_status='BLUE'),
        dict(base, number_records='5', bags_grade='E', species_status='RED',
             bin_uris='BOLD:1;BOLD:2', other_names='Aus cus;Aus sp. 1'),
        # Pre-ORANGE file: RED with only an interim sharer becomes shared_bin_interim
        dict(base, number_records='5', bags_grade='E', species_status='RED',
             bin_uris='BOLD:1', other_names='Aus cf. bus'),
        dict(base, number_records='3', bags_grade='E', species_status='ORANGE',
             bin_uris='BOLD:1', other_names='Aus sp.'),
    ])
    assert conv.convert_file(path, path)
    rows = read_tsv(path)
    header = list(rows[0])
    assert header.index('issues') == header.index('species_status') + 1
    assert [(r['species_status'], r['issues']) for r in rows] == [
        ('no_records', ''),
        ('valid_name', ''),
        ('valid_and_synonym', 'synonym_records;few_records'),
        ('synonym_only', 'synonym_records;valid_name_absent;split_bins'),
        ('shared_bin_species', 'shared_bin_species;shared_bin_interim;split_bins'),
        ('shared_bin_interim', 'shared_bin_interim'),
        ('shared_bin_interim', 'shared_bin_interim'),
    ]

    # Idempotent: a second run changes nothing
    before = path.read_bytes()
    assert not conv.convert_file(path, path)
    assert path.read_bytes() == before


@pytest.mark.parametrize('extra_col, legacy, expected', [
    ('dtol_status', ['GREEN', 'BLUE', 'AMBER', 'RED', 'BLACK'],
     ['annotation_complete', 'assembly_submitted', 'raw_data_submitted',
      'biosample_submitted', 'not_in_dtol']),
    ('mitogenome_count', ['GREEN', 'BLACK'], ['mitogenome_present', 'no_mitogenome']),
])
def test_convert_dtol_and_mitogenome(tmp_path, extra_col, legacy, expected):
    path = tmp_path / 'genome.tsv'
    write_tsv(path, [{'taxon_name': 'x', extra_col: '', 'species_status': s} for s in legacy])
    conv.convert_file(path, path)
    rows = read_tsv(path)
    assert 'issues' not in rows[0]
    assert [r['species_status'] for r in rows] == expected
