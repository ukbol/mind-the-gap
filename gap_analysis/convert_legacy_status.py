#!/usr/bin/env python3
"""
Convert legacy traffic-light gap analysis outputs to descriptive status codes.

Older gap_analysis.py, dtol_status.py and mito_status.py outputs wrote
species_status as a colour (GREEN/BLUE/AMBER/ORANGE/RED/BLACK). This script
rewrites those values to the descriptive codes now produced by the scripts,
so existing results can be used without re-running the analysis.

The dataset type is detected from the columns:
  - dtol_status present       -> DToL genome output
  - mitogenome_count present  -> ENA mitogenome output
  - otherwise                 -> barcode gene gap analysis output

For gene outputs an `issues` column is also added and rebuilt as far as the
existing columns allow:
  - shared_bin_species / shared_bin_interim  from other_names (this also moves
    pre-ORANGE RED rows that only share with interim names to
    shared_bin_interim)
  - synonym_records / valid_name_absent      from AMBER / BLUE
  - split_bins / no_cluster                  from bags_grade (for grade E rows,
    split_bins is estimated from the bin_uris/otu_ids lists)
  - few_records                              from number_records

Synonym flags cannot be recovered for rows that share a BIN/OTU (RED/ORANGE),
because the old output does not record which names had records; re-run
gap_analysis.py to populate those.

Files already using the new codes are left unchanged, so the script is safe
to run repeatedly.

Usage:
  python convert_legacy_status.py final_result/*.tsv            # in place
  python convert_legacy_status.py old.tsv --output new.tsv

Author: Ben Price / Claude
"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'dtol_processing'))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'ena_processing'))

from gap_analysis import (  # noqa: E402
    LEGACY_STATUS_MAP as GENE_LEGACY_STATUS_MAP,
    MIN_RECORDS_FOR_BAGS_B,
    ISSUE_SHARED_BIN_SPECIES, ISSUE_SHARED_BIN_INTERIM, ISSUE_SYNONYM_RECORDS,
    ISSUE_VALID_NAME_ABSENT, ISSUE_SPLIT_BINS, ISSUE_NO_CLUSTER, ISSUE_FEW_RECORDS,
    STATUS_SHARED_BIN_SPECIES, STATUS_SHARED_BIN_INTERIM,
    format_issues, is_linnaean_name, normalize_species_name,
)
from dtol_status import LEGACY_STATUS_MAP as DTOL_LEGACY_STATUS_MAP  # noqa: E402
from mito_status import LEGACY_STATUS_MAP as MITO_LEGACY_STATUS_MAP  # noqa: E402

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)  # Windows compatibility


def detect_dataset_type(columns: List[str]) -> str:
    """Return 'dtol', 'mitogenome' or 'gene' based on the output columns."""
    if 'dtol_status' in columns:
        return 'dtol'
    if 'mitogenome_count' in columns:
        return 'mitogenome'
    return 'gene'


def split_list(value: Optional[str]) -> List[str]:
    """Split a ';'-joined output field into its non-empty parts."""
    return [v.strip() for v in (value or '').split(';') if v.strip()]


def to_int(value: Optional[str]) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def convert_gene_row(row: Dict[str, str], legacy: str) -> None:
    """Set species_status and issues on a legacy gene output row, in place."""
    issues: Set[str] = set()
    status = GENE_LEGACY_STATUS_MAP[legacy]
    number_records = to_int(row.get('number_records'))

    if legacy != 'BLACK':
        other_names = [normalize_species_name(n) for n in split_list(row.get('other_names'))]
        if any(is_linnaean_name(n) for n in other_names):
            issues.add(ISSUE_SHARED_BIN_SPECIES)
        if any(not is_linnaean_name(n) for n in other_names):
            issues.add(ISSUE_SHARED_BIN_INTERIM)

        if legacy in ('RED', 'ORANGE') and other_names:
            status = (STATUS_SHARED_BIN_SPECIES if ISSUE_SHARED_BIN_SPECIES in issues
                      else STATUS_SHARED_BIN_INTERIM)

        if legacy == 'AMBER':
            issues.add(ISSUE_SYNONYM_RECORDS)
        elif legacy == 'BLUE':
            issues.update({ISSUE_SYNONYM_RECORDS, ISSUE_VALID_NAME_ABSENT})

        grade = (row.get('bags_grade') or '').strip().upper()
        if grade == 'C':
            issues.add(ISSUE_SPLIT_BINS)
        elif grade == 'F':
            issues.add(ISSUE_NO_CLUSTER)
        elif grade == 'E':
            n_clusters = max(len(set(split_list(row.get('bin_uris')))),
                             len(set(split_list(row.get('otu_ids')))))
            if n_clusters > 1:
                issues.add(ISSUE_SPLIT_BINS)

        if 0 < number_records < MIN_RECORDS_FOR_BAGS_B:
            issues.add(ISSUE_FEW_RECORDS)

    row['species_status'] = status
    row['issues'] = format_issues(issues)


def convert_file(input_path: Path, output_path: Path) -> bool:
    """Convert one file. Returns True if anything changed."""
    with open(input_path, encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        columns = list(reader.fieldnames or [])
        rows = list(reader)

    if 'species_status' not in columns:
        logging.warning(f"{input_path}: no species_status column, skipping")
        return False

    dataset_type = detect_dataset_type(columns)
    legacy_map = {
        'gene': GENE_LEGACY_STATUS_MAP,
        'dtol': DTOL_LEGACY_STATUS_MAP,
        'mitogenome': MITO_LEGACY_STATUS_MAP,
    }[dataset_type]

    add_issues = dataset_type == 'gene' and 'issues' not in columns
    if add_issues:
        columns.insert(columns.index('species_status') + 1, 'issues')

    changed = add_issues
    for row in rows:
        legacy = (row.get('species_status') or '').strip().upper()
        if legacy not in legacy_map:
            if add_issues:
                row['issues'] = ''
            continue
        changed = True
        if dataset_type == 'gene':
            convert_gene_row(row, legacy)
        else:
            row['species_status'] = legacy_map[legacy]

    if not changed and output_path == input_path:
        logging.info(f"{input_path}: already uses status codes, unchanged")
        return False

    tmp_path = output_path.with_name(output_path.name + '.tmp')
    with open(tmp_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    tmp_path.replace(output_path)

    counts = Counter(r.get('species_status', '') for r in rows)
    summary = ', '.join(f"{k}={v:,}" for k, v in counts.most_common())
    logging.info(f"{input_path} -> {output_path} ({dataset_type}): {summary}")
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert legacy colour species_status values to descriptive codes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage:')[0].split('\n', 2)[2],
    )
    parser.add_argument('inputs', nargs='+', type=Path, help='Gap analysis TSV file(s)')
    parser.add_argument('--output', type=Path,
                        help='Output path (only with a single input; default: overwrite input)')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')

    if args.output and len(args.inputs) != 1:
        parser.error('--output can only be used with a single input file')

    for input_path in args.inputs:
        convert_file(input_path, args.output or input_path)


if __name__ == '__main__':
    main()
