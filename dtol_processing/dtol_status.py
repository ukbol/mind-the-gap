#!/usr/bin/env python3
"""
DToL (Darwin Tree of Life) Status Processing for UK Species.

This script matches DToL genome sequencing metadata against a UKSI species
list (the same list used by gap_analysis.py) to produce a per-species
status panel showing how far each UK species has progressed through the
DToL genome sequencing pipeline.

Matching uses valid names AND synonyms so that nomenclatural differences
between DToL and UKSI are bridged automatically.

Output is a TSV file that preserves all input species-list columns and
appends DToL-specific columns, following the same pattern as
gap_analysis.py output.

Author: Ben Price / Claude
Date: 2025-02-14
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, field

# Increase CSV field size limit for large files
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)


# =============================================================================
# CONSTANTS
# =============================================================================

# DToL pipeline stages in order of progress (furthest first).
# Used to determine the "best" status when multiple DToL records
# match a single UKSI taxon.
DTOL_STAGE_ORDER = [
    'Annotation Complete',
    'Assemblies - Submitted',
    'Raw Data - Submitted',
    'Submitted to BioSamples',
]

# Mapping from DToL pipeline stage to a traffic-light label,
# giving visual consistency with the barcode gap analysis panels.
DTOL_STATUS_COLOURS = {
    'Annotation Complete':    'GREEN',
    'Assemblies - Submitted': 'AMBER',
    'Raw Data - Submitted':   'BLUE',
    'Submitted to BioSamples':'RED',
}


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Taxon:
    """A UKSI taxon with valid name, synonyms, and original row data."""
    row_index: int
    valid_name: str
    synonyms: List[str]
    input_data: Dict[str, str]

    @property
    def all_names(self) -> Set[str]:
        names = {normalize_name(self.valid_name)}
        names.update(normalize_name(s) for s in self.synonyms)
        return names


@dataclass
class DtolRecord:
    """A single row from the DToL metadata CSV."""
    organism: str
    common_name: str
    current_status: str
    insdc_id: str
    tol_ids: List[str]


@dataclass
class TaxonResult:
    """DToL analysis results for a single UKSI taxon."""
    taxon: Taxon
    dtol_status: str = 'Not in DToL'
    dtol_status_colour: str = 'BLACK'
    dtol_organism_name: str = ''
    dtol_common_name: str = ''
    dtol_insdc_ids: List[str] = field(default_factory=list)
    dtol_tol_ids: List[str] = field(default_factory=list)
    dtol_matched_name: str = ''  # which name (valid/synonym) matched


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def normalize_name(name: str) -> str:
    """Normalize a species name for matching (lowercase, underscores→spaces)."""
    return name.strip().lower().replace('_', ' ')


def setup_logging(log_level: str = "INFO") -> None:
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )


def stage_rank(status: str) -> int:
    """Return a numeric rank for a DToL stage (lower = further along)."""
    try:
        return DTOL_STAGE_ORDER.index(status)
    except ValueError:
        return len(DTOL_STAGE_ORDER)  # unknown stages ranked last


# =============================================================================
# FILE LOADING
# =============================================================================

def load_species_list(species_file: Path) -> Tuple[List[Taxon], List[str]]:
    """
    Load the UKSI species list TSV (same format as gap_analysis.py input).

    Returns (list of Taxon objects, list of input column names).
    """
    logging.info(f"Loading species list from {species_file}")
    taxa = []
    input_columns: List[str] = []

    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(species_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')
                input_columns = list(reader.fieldnames) if reader.fieldnames else []

                species_col = None
                if 'taxon_name' in input_columns:
                    species_col = 'taxon_name'
                elif 'species' in input_columns:
                    species_col = 'species'
                else:
                    logging.error("Species list must have a 'taxon_name' or 'species' column")
                    sys.exit(1)

                logging.info(f"Using '{species_col}' column for valid species names")

                for row_idx, row in enumerate(reader):
                    valid_name = row.get(species_col, '').strip()
                    if not valid_name:
                        continue

                    synonyms_field = row.get('synonyms', '').strip()
                    synonyms = [s.strip() for s in synonyms_field.split(';') if s.strip()]

                    input_data = {col: row.get(col, '').strip() for col in input_columns}

                    taxa.append(Taxon(
                        row_index=row_idx,
                        valid_name=valid_name,
                        synonyms=synonyms,
                        input_data=input_data,
                    ))

            logging.info(f"Loaded {len(taxa):,} taxa ({sum(len(t.synonyms) for t in taxa):,} synonyms)")
            return taxa, input_columns

        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                continue
            raise

    logging.error("Failed to load species list")
    sys.exit(1)


def load_dtol_metadata(dtol_file: Path) -> List[DtolRecord]:
    """
    Load the DToL portal CSV export.

    Expected columns: Organism, Common Name, Current Status, INSDC ID, ToL ID
    """
    logging.info(f"Loading DToL metadata from {dtol_file}")
    records: List[DtolRecord] = []

    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(dtol_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f)

                if reader.fieldnames is None:
                    logging.error("DToL file appears to be empty")
                    sys.exit(1)

                # Validate expected columns
                expected = {'Organism', 'Current Status'}
                missing = expected - set(reader.fieldnames)
                if missing:
                    logging.error(f"DToL file missing required columns: {missing}")
                    sys.exit(1)

                for row in reader:
                    organism = row.get('Organism', '').strip()
                    if not organism:
                        continue

                    common_name = row.get('Common Name', '').strip()
                    if common_name == 'Not specified':
                        common_name = ''

                    tol_id_raw = row.get('ToL ID', '').strip()
                    tol_ids = [t.strip() for t in tol_id_raw.split(',') if t.strip()]

                    records.append(DtolRecord(
                        organism=organism,
                        common_name=common_name,
                        current_status=row.get('Current Status', '').strip(),
                        insdc_id=row.get('INSDC ID', '').strip(),
                        tol_ids=tol_ids,
                    ))

            logging.info(f"Loaded {len(records):,} DToL records")

            # Log status breakdown
            from collections import Counter
            status_counts = Counter(r.current_status for r in records)
            for status, count in status_counts.most_common():
                logging.info(f"  {status}: {count:,}")

            return records

        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                continue
            raise

    logging.error("Failed to load DToL metadata")
    sys.exit(1)


# =============================================================================
# ANALYSIS
# =============================================================================

def build_dtol_index(records: List[DtolRecord]) -> Dict[str, List[DtolRecord]]:
    """
    Build a lookup from normalized organism name → list of DToL records.

    A single organism name should normally appear only once in the DToL
    portal export, but we handle duplicates gracefully.
    """
    index: Dict[str, List[DtolRecord]] = {}
    for rec in records:
        key = normalize_name(rec.organism)
        index.setdefault(key, []).append(rec)
    logging.info(f"Built DToL index with {len(index):,} unique organism names")
    return index


def analyze_taxa(
    taxa: List[Taxon],
    dtol_index: Dict[str, List[DtolRecord]],
) -> Tuple[List[TaxonResult], List[DtolRecord]]:
    """
    For each UKSI taxon, find the best-matching DToL record(s).

    Matching logic:
    1. Check valid name against DToL organism names
    2. Check each synonym against DToL organism names
    3. If multiple DToL records match (via different names), keep the one
       at the furthest pipeline stage.

    Returns:
        Tuple of (list of TaxonResult, list of unmatched DtolRecords)
    """
    results: List[TaxonResult] = []
    matched_count = 0

    # Build set of all UKSI names (valid + synonyms) for unmatched detection
    all_uksi_names: Set[str] = set()
    for taxon in taxa:
        all_uksi_names.add(normalize_name(taxon.valid_name))
        for syn in taxon.synonyms:
            all_uksi_names.add(normalize_name(syn))

    for taxon in taxa:
        result = TaxonResult(taxon=taxon)

        # Collect all matching DToL records across all names
        matches: List[Tuple[str, DtolRecord]] = []  # (matched_name, record)

        # Check valid name first, then synonyms
        names_to_check = [taxon.valid_name] + taxon.synonyms
        for name in names_to_check:
            key = normalize_name(name)
            if key in dtol_index:
                for rec in dtol_index[key]:
                    matches.append((name, rec))

        if matches:
            matched_count += 1

            # Pick the best record (furthest pipeline stage)
            matches.sort(key=lambda m: stage_rank(m[1].current_status))
            best_name, best_rec = matches[0]

            result.dtol_status = best_rec.current_status
            result.dtol_status_colour = DTOL_STATUS_COLOURS.get(
                best_rec.current_status, 'BLACK'
            )
            result.dtol_organism_name = best_rec.organism
            result.dtol_common_name = best_rec.common_name
            result.dtol_matched_name = best_name

            # Aggregate INSDC IDs and ToL IDs from ALL matches
            seen_insdc: Set[str] = set()
            seen_tol: Set[str] = set()
            for _, rec in matches:
                if rec.insdc_id and rec.insdc_id not in seen_insdc:
                    result.dtol_insdc_ids.append(rec.insdc_id)
                    seen_insdc.add(rec.insdc_id)
                for tid in rec.tol_ids:
                    if tid not in seen_tol:
                        result.dtol_tol_ids.append(tid)
                        seen_tol.add(tid)

        results.append(result)

    logging.info(f"Matched {matched_count:,} / {len(taxa):,} UKSI taxa to DToL records")

    # Collect unmatched DToL records
    unmatched: List[DtolRecord] = []
    for norm_name, recs in dtol_index.items():
        if norm_name not in all_uksi_names:
            unmatched.extend(recs)

    if unmatched:
        logging.info(f"{len(unmatched):,} DToL organisms did not match any UKSI taxon")

    return results, unmatched


# =============================================================================
# OUTPUT
# =============================================================================

def write_results(
    results: List[TaxonResult],
    output_file: Path,
    input_columns: List[str],
) -> None:
    """Write results as TSV, preserving input columns + appending DToL columns."""
    logging.info(f"Writing results to {output_file}")

    dtol_columns = [
        'dtol_status',
        'dtol_status_colour',
        'dtol_organism_name',
        'dtol_common_name',
        'dtol_insdc_ids',
        'dtol_tol_ids',
        'dtol_matched_name',
    ]

    output_columns = list(input_columns) + dtol_columns

    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(
                f,
                fieldnames=output_columns,
                delimiter='\t',
                extrasaction='ignore',
            )
            writer.writeheader()

            for result in results:
                row = dict(result.taxon.input_data)
                row['dtol_status'] = result.dtol_status
                row['dtol_status_colour'] = result.dtol_status_colour
                row['dtol_organism_name'] = result.dtol_organism_name
                row['dtol_common_name'] = result.dtol_common_name
                row['dtol_insdc_ids'] = ';'.join(result.dtol_insdc_ids)
                row['dtol_tol_ids'] = ';'.join(result.dtol_tol_ids)
                row['dtol_matched_name'] = result.dtol_matched_name
                writer.writerow(row)

        logging.info(f"Successfully wrote {len(results):,} results")

    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


def write_unmatched(
    unmatched: List[DtolRecord],
    output_file: Path,
) -> None:
    """Write unmatched DToL organisms to a TSV for manual follow-up."""
    if not unmatched:
        logging.info("No unmatched DToL organisms to write")
        return

    logging.info(f"Writing {len(unmatched):,} unmatched DToL organisms to {output_file}")

    columns = ['organism', 'common_name', 'current_status', 'insdc_id', 'tol_ids']

    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            writer.writeheader()

            for rec in sorted(unmatched, key=lambda r: r.organism):
                writer.writerow({
                    'organism': rec.organism,
                    'common_name': rec.common_name,
                    'current_status': rec.current_status,
                    'insdc_id': rec.insdc_id,
                    'tol_ids': ';'.join(rec.tol_ids),
                })

        logging.info(f"Successfully wrote {len(unmatched):,} unmatched records")

    except Exception as e:
        logging.error(f"Failed to write unmatched file: {e}")
        sys.exit(1)


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def print_summary(results: List[TaxonResult]) -> None:
    """Print a summary of the analysis to the log."""
    from collections import Counter

    total = len(results)
    colour_counts = Counter(r.dtol_status_colour for r in results)
    status_counts = Counter(r.dtol_status for r in results)

    logging.info("=" * 60)
    logging.info("DToL STATUS SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Total UKSI taxa: {total:,}")
    logging.info("")
    logging.info("By pipeline stage:")
    for status in DTOL_STAGE_ORDER:
        count = status_counts.get(status, 0)
        pct = (count / total * 100) if total else 0
        logging.info(f"  {status:30s}  {count:>6,}  ({pct:5.1f}%)")
    not_in = status_counts.get('Not in DToL', 0)
    pct = (not_in / total * 100) if total else 0
    logging.info(f"  {'Not in DToL':30s}  {not_in:>6,}  ({pct:5.1f}%)")
    logging.info("")
    logging.info("By status colour:")
    for colour in ['GREEN', 'AMBER', 'BLUE', 'RED', 'BLACK']:
        count = colour_counts.get(colour, 0)
        pct = (count / total * 100) if total else 0
        logging.info(f"  {colour:10s}  {count:>6,}  ({pct:5.1f}%)")
    logging.info("=" * 60)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Process DToL metadata against a UKSI species list.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dtol_status.py \\
      --species-list uksi_species_export.tsv \\
      --dtol-metadata download.csv \\
      --output dtol_gap_analysis.tsv

Status colour mapping:
  GREEN  = Annotation Complete (genome published)
  AMBER  = Assemblies Submitted (assembly in progress)
  BLUE   = Raw Data Submitted (sequencing done)
  RED    = Submitted to BioSamples (sample registered)
  BLACK  = Not in DToL pipeline
        """,
    )
    parser.add_argument(
        '--species-list', required=True, type=Path,
        help='Path to UKSI species list TSV (with taxon_name/species + synonyms)',
    )
    parser.add_argument(
        '--dtol-metadata', required=True, type=Path,
        help='Path to DToL portal CSV export',
    )
    parser.add_argument(
        '--output', required=True, type=Path,
        help='Output path for DToL status TSV',
    )
    parser.add_argument(
        '--log-level', default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)',
    )

    args = parser.parse_args()
    setup_logging(args.log_level)

    # Validate input files exist
    if not args.species_list.exists():
        logging.error(f"Species list not found: {args.species_list}")
        sys.exit(1)
    if not args.dtol_metadata.exists():
        logging.error(f"DToL metadata not found: {args.dtol_metadata}")
        sys.exit(1)

    # Load data
    taxa, input_columns = load_species_list(args.species_list)
    dtol_records = load_dtol_metadata(args.dtol_metadata)

    # Build index and analyze
    dtol_index = build_dtol_index(dtol_records)
    results, unmatched = analyze_taxa(taxa, dtol_index)

    # Output
    write_results(results, args.output, input_columns)

    # Write unmatched DToL organisms alongside the main output
    unmatched_file = args.output.parent / (args.output.stem + '_unmatched.tsv')
    write_unmatched(unmatched, unmatched_file)

    print_summary(results)

    logging.info("Done.")


if __name__ == '__main__':
    main()
