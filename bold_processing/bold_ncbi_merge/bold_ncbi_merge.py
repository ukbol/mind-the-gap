#!/usr/bin/env python3
"""
bold_ncbi_merge.py — Merge GenBank and BOLD TSV records, preferring BOLD.

Deduplicates on GenBank accession number:
  - BOLD records use the 'insdc_acs' column
  - GenBank records use the 'locus_name' column (= accession)

If a record exists in both files, the BOLD version is kept.
GenBank-only records are mapped to BOLD columns and appended.

Usage:
    python bold_ncbi_merge.py <genbank_tsv> <bold_tsv> <output_tsv>

    python bold_ncbi_merge.py sequences.tsv bold_plants.tsv merged.tsv

Options:
    -v, --verbose    Print progress and statistics to stderr

Author: Ben Price / Claude
Date: 2026-03-28
"""

import argparse
import csv
import sys
from pathlib import Path

# Mapping from GenBank column names to BOLD column names.
# Only fields with a clear equivalent are mapped; the rest are left empty.
GENBANK_TO_BOLD = {
    'locus_name':             'insdc_acs',
    'organism':               'identification',
    'source_collected_by':    'collectors',
    'source_collection_date': 'collection_date_start',
    'source_geo_loc_name':    'site',
    'source_identified_by':   'identified_by',
    'source_specimen_voucher':'museumid',
    'nucleotide_sequence':    'nuc',
    'locus_length':           'nuc_basecount',
}

# Source tag written into the 'notes' column so merged records are traceable
GENBANK_SOURCE_TAG = 'source:GenBank'


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge GenBank and BOLD TSV records, preferring BOLD on duplicate accessions.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('genbank', help='GenBank TSV file (locus_name = accession)')
    parser.add_argument('bold',    help='BOLD TSV file (insdc_acs = accession)')
    parser.add_argument('output',  help='Output merged TSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress and statistics to stderr')
    return parser.parse_args()


def load_bold(bold_path: Path, verbose: bool) -> tuple[list[dict], list[str], set[str]]:
    """
    Load all BOLD records into memory.

    Returns:
        records   : list of row dicts
        fieldnames: ordered list of column names
        accessions: set of insdc_acs values present (empty strings excluded)
    """
    records = []
    accessions = set()

    with open(bold_path, encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'insdc_acs' not in reader.fieldnames:
            raise ValueError(f"'insdc_acs' column not found in BOLD file: {bold_path}")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            records.append(row)
            acc = row.get('insdc_acs', '').strip()
            if acc:
                accessions.add(acc)

    if verbose:
        print(f"BOLD: loaded {len(records):,} records, "
              f"{len(accessions):,} with accession numbers", file=sys.stderr)
    return records, fieldnames, accessions


def map_genbank_row(gb_row: dict, bold_fieldnames: list[str]) -> dict:
    """
    Convert a GenBank row dict into a BOLD-schema row dict.
    Unmapped columns are left as empty strings.
    Adds a source tag to the 'notes' field for traceability.
    """
    bold_row = {col: '' for col in bold_fieldnames}

    for gb_col, bold_col in GENBANK_TO_BOLD.items():
        if bold_col in bold_row:
            bold_row[bold_col] = gb_row.get(gb_col, '').strip()

    # Traceability tag
    existing_notes = bold_row.get('notes', '')
    bold_row['notes'] = (
        f"{existing_notes}; {GENBANK_SOURCE_TAG}".lstrip('; ')
        if existing_notes else GENBANK_SOURCE_TAG
    )

    return bold_row


def load_genbank_unique(
    gb_path: Path,
    bold_accessions: set[str],
    bold_fieldnames: list[str],
    verbose: bool
) -> tuple[list[dict], int, int]:
    """
    Stream the GenBank TSV and return only records whose accession
    is NOT already present in the BOLD dataset.

    Returns:
        unique_rows   : list of mapped BOLD-schema dicts
        total_gb      : total GenBank records read
        duplicate_count: records skipped as duplicates
    """
    unique_rows = []
    total_gb = 0
    duplicate_count = 0

    with open(gb_path, encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'locus_name' not in reader.fieldnames:
            raise ValueError(f"'locus_name' column not found in GenBank file: {gb_path}")

        for row in reader:
            total_gb += 1
            acc = row.get('locus_name', '').strip()
            if acc and acc in bold_accessions:
                duplicate_count += 1
                continue
            unique_rows.append(map_genbank_row(row, bold_fieldnames))

    if verbose:
        print(f"GenBank: {total_gb:,} records read, "
              f"{duplicate_count:,} duplicates skipped, "
              f"{len(unique_rows):,} unique records to append", file=sys.stderr)

    return unique_rows, total_gb, duplicate_count


def write_output(
    output_path: Path,
    bold_records: list[dict],
    genbank_unique: list[dict],
    fieldnames: list[str],
    verbose: bool
) -> None:
    """Write all BOLD records followed by unique GenBank records to output TSV."""
    with open(output_path, 'w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter='\t',
            extrasaction='ignore',
            lineterminator='\n'
        )
        writer.writeheader()

        for row in bold_records:
            writer.writerow(row)

        for row in genbank_unique:
            writer.writerow(row)

    total = len(bold_records) + len(genbank_unique)
    if verbose:
        print(f"Output: {total:,} records written to {output_path}", file=sys.stderr)
        print(f"  BOLD records : {len(bold_records):,}", file=sys.stderr)
        print(f"  GenBank-only : {len(genbank_unique):,}", file=sys.stderr)


def main():
    args = parse_args()

    gb_path  = Path(args.genbank)
    bold_path = Path(args.bold)
    out_path  = Path(args.output)

    for p in (gb_path, bold_path):
        if not p.exists():
            print(f"ERROR: File not found: {p}", file=sys.stderr)
            sys.exit(1)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        bold_records, fieldnames, bold_accessions = load_bold(bold_path, args.verbose)
        genbank_unique, _, _ = load_genbank_unique(
            gb_path, bold_accessions, fieldnames, args.verbose
        )
        write_output(out_path, bold_records, genbank_unique, fieldnames, args.verbose)

    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if not args.verbose:
        total = len(bold_records) + len(genbank_unique)
        print(f"Done. {total:,} records written to {out_path}")


if __name__ == '__main__':
    main()
