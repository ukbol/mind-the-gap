#!/usr/bin/env python3
"""
ENA Mitogenome Status Processing for UK Species.

This script matches ENA mitogenome accession data against a species list
(the same list used by gap_analysis.py) to produce a per-species count
showing how many mitogenome accessions are available for each taxon.

Matching uses valid names AND synonyms so that nomenclatural differences
between ENA and the species list are bridged automatically.

Output is a TSV file that preserves all input species-list columns and
appends mitogenome-specific columns, following the same pattern as
dtol_status.py and gap_analysis.py output.

Author: Ben Price / Claude
Date: 2025-02-17
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass, field
from collections import Counter, defaultdict

# Increase CSV field size limit for large files
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Taxon:
    """A taxon with valid name, synonyms, and original row data."""
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
class MitoRecord:
    """A single row from the ENA mitogenome TSV."""
    accession: str
    scientific_name: str
    tax_id: str


@dataclass
class TaxonResult:
    """Mitogenome analysis results for a single taxon."""
    taxon: Taxon
    mitogenome_count: int = 0
    species_status: str = 'BLACK'
    mitogenome_accessions: List[str] = field(default_factory=list)
    matched_names: Set[str] = field(default_factory=set)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def normalize_name(name: str) -> str:
    """Normalize a species name for matching (lowercase, underscores to spaces)."""
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


# =============================================================================
# FILE LOADING
# =============================================================================

def load_species_list(species_file: Path) -> Tuple[List[Taxon], List[str]]:
    """
    Load species list from TSV file.

    Expected columns:
    - 'taxon_name' or 'species': Valid name (required, prefers taxon_name)
    - 'synonyms': Semicolon-separated synonyms (optional)

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


def load_mito_metadata(mito_file: Path) -> List[MitoRecord]:
    """
    Load ENA mitogenome TSV export.

    Expected columns: accession, scientific_name, tax_id
    """
    logging.info(f"Loading mitogenome metadata from {mito_file}")
    records: List[MitoRecord] = []

    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(mito_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')

                if reader.fieldnames is None:
                    logging.error("Mitogenome file appears to be empty")
                    sys.exit(1)

                expected = {'accession', 'scientific_name'}
                missing = expected - set(reader.fieldnames)
                if missing:
                    logging.error(f"Mitogenome file missing required columns: {missing}")
                    sys.exit(1)

                for row in reader:
                    sci_name = row.get('scientific_name', '').strip()
                    if not sci_name:
                        continue

                    records.append(MitoRecord(
                        accession=row.get('accession', '').strip(),
                        scientific_name=sci_name,
                        tax_id=row.get('tax_id', '').strip(),
                    ))

            logging.info(f"Loaded {len(records):,} mitogenome records")

            # Log basic stats
            unique_species = len(set(r.scientific_name for r in records))
            logging.info(f"  Unique species names: {unique_species:,}")

            return records

        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                continue
            raise

    logging.error("Failed to load mitogenome metadata")
    sys.exit(1)


# =============================================================================
# ANALYSIS
# =============================================================================

def build_mito_index(records: List[MitoRecord]) -> Dict[str, List[MitoRecord]]:
    """
    Build a lookup from normalized scientific name to list of MitoRecords.

    Multiple accessions may exist for a single species.
    """
    index: Dict[str, List[MitoRecord]] = defaultdict(list)
    for rec in records:
        key = normalize_name(rec.scientific_name)
        index[key].append(rec)
    logging.info(f"Built mitogenome index with {len(index):,} unique species names")
    return index


def analyze_taxa(
    taxa: List[Taxon],
    mito_index: Dict[str, List[MitoRecord]],
) -> Tuple[List[TaxonResult], List[MitoRecord]]:
    """
    For each taxon, find all matching mitogenome records.

    Matching logic:
    1. Check valid name against mitogenome species names
    2. Check each synonym against mitogenome species names
    3. Aggregate all matching accessions across all names

    Returns:
        Tuple of (list of TaxonResult, list of unmatched MitoRecords)
    """
    results: List[TaxonResult] = []
    matched_count = 0

    # Build set of all species-list names for unmatched detection
    all_taxa_names: Set[str] = set()
    for taxon in taxa:
        all_taxa_names.update(taxon.all_names)

    for taxon in taxa:
        result = TaxonResult(taxon=taxon)

        # Check all names (valid + synonyms) against mito index
        for name in taxon.all_names:
            if name in mito_index:
                result.matched_names.add(name)
                for rec in mito_index[name]:
                    result.mitogenome_accessions.append(rec.accession)

        if result.mitogenome_accessions:
            matched_count += 1
            result.mitogenome_count = len(result.mitogenome_accessions)
            result.species_status = 'GREEN'

        results.append(result)

    logging.info(f"Matched {matched_count:,} / {len(taxa):,} taxa to mitogenome records")

    # Collect unmatched mitogenome records
    unmatched: List[MitoRecord] = []
    for norm_name, recs in mito_index.items():
        if norm_name not in all_taxa_names:
            unmatched.extend(recs)

    if unmatched:
        unique_unmatched = len(set(r.scientific_name for r in unmatched))
        logging.info(
            f"{len(unmatched):,} mitogenome accessions ({unique_unmatched:,} species) "
            f"did not match any taxon in the species list"
        )

    return results, unmatched


# =============================================================================
# OUTPUT
# =============================================================================

def write_results(
    results: List[TaxonResult],
    output_file: Path,
    input_columns: List[str],
) -> None:
    """Write results as TSV, preserving input columns + appending mito columns."""
    logging.info(f"Writing results to {output_file}")

    mito_columns = [
        'mitogenome_count',
        'species_status',
        'mitogenome_accessions',
    ]

    output_columns = list(input_columns) + mito_columns

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
                row['mitogenome_count'] = result.mitogenome_count
                row['species_status'] = result.species_status
                row['mitogenome_accessions'] = ';'.join(sorted(result.mitogenome_accessions))
                writer.writerow(row)

        logging.info(f"Successfully wrote {len(results):,} results")

    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


def write_unmatched(
    unmatched: List[MitoRecord],
    output_file: Path,
) -> None:
    """Write unmatched mitogenome records to a TSV for manual follow-up."""
    if not unmatched:
        logging.info("No unmatched mitogenome records to write")
        return

    logging.info(f"Writing {len(unmatched):,} unmatched mitogenome records to {output_file}")

    columns = ['accession', 'scientific_name', 'tax_id']

    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            writer.writeheader()

            for rec in sorted(unmatched, key=lambda r: r.scientific_name):
                writer.writerow({
                    'accession': rec.accession,
                    'scientific_name': rec.scientific_name,
                    'tax_id': rec.tax_id,
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
    total = len(results)
    with_mito = sum(1 for r in results if r.mitogenome_count > 0)
    without_mito = total - with_mito
    total_accessions = sum(r.mitogenome_count for r in results)

    logging.info("=" * 60)
    logging.info("MITOGENOME STATUS SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Total taxa: {total:,}")
    logging.info(f"Taxa with mitogenomes:    {with_mito:>6,}  ({100 * with_mito / total:.1f}%)")
    logging.info(f"Taxa without mitogenomes: {without_mito:>6,}  ({100 * without_mito / total:.1f}%)")
    logging.info(f"Total accessions matched: {total_accessions:,}")
    logging.info("")

    # Distribution of accession counts among taxa that have mitogenomes
    if with_mito > 0:
        counts = [r.mitogenome_count for r in results if r.mitogenome_count > 0]
        counts.sort()
        logging.info("Accession count distribution (taxa with mitogenomes):")
        logging.info(f"  Min: {counts[0]}")
        logging.info(f"  Median: {counts[len(counts) // 2]}")
        logging.info(f"  Max: {counts[-1]}")

        # Bracket breakdown
        brackets = [(1, 1), (2, 5), (6, 10), (11, 50), (51, None)]
        for lo, hi in brackets:
            if hi is None:
                n = sum(1 for c in counts if c >= lo)
                label = f"  {lo}+ accessions:"
            else:
                n = sum(1 for c in counts if lo <= c <= hi)
                label = f"  {lo}-{hi} accessions:" if lo != hi else f"  {lo} accession:"
            logging.info(f"{label:30s} {n:>5,}")

    logging.info("=" * 60)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Process ENA mitogenome metadata against a species list.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mito_status.py \\
      --species-list uksi_species_export.tsv \\
      --mito-metadata mito_species_ena.tsv \\
      --output mito_gap_analysis.tsv

Status colour mapping:
  GREEN = Mitogenome(s) available in ENA
  BLACK = No mitogenome found in ENA
        """,
    )
    parser.add_argument(
        '--species-list', required=True, type=Path,
        help='Path to species list TSV (with taxon_name/species + synonyms columns)',
    )
    parser.add_argument(
        '--mito-metadata', required=True, type=Path,
        help='Path to ENA mitogenome TSV (with accession, scientific_name, tax_id columns)',
    )
    parser.add_argument(
        '--output', required=True, type=Path,
        help='Output path for mitogenome status TSV',
    )
    parser.add_argument(
        '--log-level', default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)',
    )

    args = parser.parse_args()
    setup_logging(args.log_level)

    logging.info("=" * 60)
    logging.info("ENA Mitogenome Status Processing")
    logging.info("=" * 60)

    # Validate input files exist
    if not args.species_list.exists():
        logging.error(f"Species list not found: {args.species_list}")
        sys.exit(1)
    if not args.mito_metadata.exists():
        logging.error(f"Mitogenome metadata not found: {args.mito_metadata}")
        sys.exit(1)

    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Load data
    taxa, input_columns = load_species_list(args.species_list)
    mito_records = load_mito_metadata(args.mito_metadata)

    # Build index and analyze
    mito_index = build_mito_index(mito_records)
    results, unmatched = analyze_taxa(taxa, mito_index)

    # Output
    write_results(results, args.output, input_columns)

    # Write unmatched mitogenome records alongside the main output
    unmatched_file = args.output.parent / (args.output.stem + '_unmatched.tsv')
    write_unmatched(unmatched, unmatched_file)

    print_summary(results)

    logging.info("Done.")


if __name__ == '__main__':
    main()
