#!/usr/bin/env python3
"""
HPC-optimized Gap Analysis for DNA Barcode Library Curation.

This script performs taxon-centric gap analysis by:
1. Loading a species list (taxa with valid names and synonyms)
2. Scanning a records file to find all records matching each taxon's names
3. Analyzing BIN/OTU sharing to detect taxonomic conflicts
4. Assigning BAGS grades (A-F) and traffic light status (GREEN/AMBER/RED/BLUE/BLACK)

Optimized for HPC environments with parallel processing support.

Author: Ben Price / Claude
Date: 2025-01-25
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
import multiprocessing as mp
import time

# Increase CSV field size limit for large files
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2147483647)  # Windows compatibility


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Taxon:
    """Represents a taxon with valid name and synonyms."""
    row_index: int
    valid_name: str
    synonyms: List[str]
    input_data: Dict[str, str]  # All columns from input file
    
    @property
    def all_names(self) -> Set[str]:
        """All names (valid + synonyms) as lowercase set."""
        names = {self.valid_name.lower()}
        names.update(s.lower() for s in self.synonyms)
        return names
    
    @property
    def all_names_list(self) -> List[str]:
        """All names as list, valid name first."""
        return [self.valid_name] + self.synonyms


@dataclass
class TaxonResult:
    """Analysis results for a single taxon."""
    taxon: Taxon
    number_records: int = 0
    bags_grade: str = 'F'
    species_status: str = 'BLACK'
    other_names: List[str] = field(default_factory=list)
    bins_found: Set[str] = field(default_factory=set)
    names_recorded: Set[str] = field(default_factory=set)  # Which of taxon's names have records


# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure logging with timestamp."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')
    
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# =============================================================================
# FILE LOADING
# =============================================================================

def load_species_list(species_file: Path) -> Tuple[List[Taxon], List[str]]:
    """
    Load species list from TSV file.
    
    Expected columns:
    - 'species': Valid name (required)
    - 'synonyms': Semicolon-separated synonyms (optional, can be empty)
    - Additional columns are preserved
    
    Returns:
        Tuple of (list of Taxon objects, list of input column names)
    """
    logging.info(f"Loading species list from {species_file}")
    taxa = []
    input_columns = []
    
    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(species_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')
                input_columns = list(reader.fieldnames) if reader.fieldnames else []
                
                if 'species' not in input_columns:
                    logging.error("Species list must have a 'species' column")
                    sys.exit(1)
                
                for row_idx, row in enumerate(reader):
                    valid_name = row.get('species', '').strip()
                    if not valid_name:
                        continue
                    
                    # Parse synonyms
                    synonyms_field = row.get('synonyms', '').strip()
                    synonyms = [s.strip() for s in synonyms_field.split(';') if s.strip()]
                    
                    # Store all input data
                    input_data = {col: row.get(col, '').strip() for col in input_columns}
                    
                    taxa.append(Taxon(
                        row_index=row_idx,
                        valid_name=valid_name,
                        synonyms=synonyms,
                        input_data=input_data
                    ))
            
            logging.info(f"Loaded {len(taxa)} taxa from species list")
            total_synonyms = sum(len(t.synonyms) for t in taxa)
            logging.info(f"Total synonyms: {total_synonyms}")
            total_names = len(taxa) + total_synonyms
            logging.info(f"Total names to match: {total_names}")
            return taxa, input_columns
            
        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                continue
            raise
    
    logging.error("Failed to load species list")
    sys.exit(1)


def detect_cluster_column(fieldnames: List[str]) -> str:
    """
    Detect which column to use for cluster identity (bin_uri or otu_id).
    
    Returns column name to use.
    """
    if 'bin_uri' in fieldnames:
        return 'bin_uri'
    elif 'otu_id' in fieldnames:
        return 'otu_id'
    elif 'OTU_ID' in fieldnames:
        return 'OTU_ID'
    elif 'BIN' in fieldnames:
        return 'BIN'
    else:
        logging.error("Records file must have 'bin_uri' or 'otu_id' column")
        sys.exit(1)


def build_indices_from_records(
    records_file: Path,
    chunk_size: int = 100000
) -> Tuple[Dict[str, int], Dict[str, Set[str]], Dict[str, Set[str]], str]:
    """
    Build indices from records file in a single pass.
    
    Returns:
        - name_to_count: species_name (lowercase) -> record count
        - name_to_bins: species_name (lowercase) -> set of BIN/OTU IDs
        - bin_to_names: BIN/OTU ID -> set of species_names (lowercase)
        - cluster_column: name of column used for clustering
    """
    logging.info(f"Building indices from {records_file}")
    start_time = time.time()
    
    name_to_count: Dict[str, int] = defaultdict(int)
    name_to_bins: Dict[str, Set[str]] = defaultdict(set)
    bin_to_names: Dict[str, Set[str]] = defaultdict(set)
    
    total_records = 0
    records_with_cluster = 0
    unique_names = set()
    unique_bins = set()
    
    cluster_column = None
    has_subspecies = False
    
    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(records_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')
                fieldnames = reader.fieldnames or []
                
                # Detect cluster column (file-level decision)
                cluster_column = detect_cluster_column(fieldnames)
                logging.info(f"Using cluster column: {cluster_column}")
                
                # Check for subspecies column
                has_subspecies = 'subspecies' in fieldnames
                if has_subspecies:
                    logging.info("Subspecies column detected - will match against both species and subspecies")
                
                for row in reader:
                    total_records += 1
                    
                    # Get species name
                    species = row.get('species', '').strip()
                    if not species:
                        continue
                    
                    # Get cluster IDs (may be pipe-separated)
                    cluster_field = row.get(cluster_column, '').strip()
                    cluster_ids = [c.strip() for c in cluster_field.split('|') if c.strip()]
                    
                    # Process species name
                    species_lower = species.lower()
                    name_to_count[species_lower] += 1
                    unique_names.add(species_lower)
                    
                    if cluster_ids:
                        records_with_cluster += 1
                        for cid in cluster_ids:
                            name_to_bins[species_lower].add(cid)
                            bin_to_names[cid].add(species_lower)
                            unique_bins.add(cid)
                    
                    # Also process subspecies if present
                    if has_subspecies:
                        subspecies = row.get('subspecies', '').strip()
                        if subspecies and subspecies.lower() not in ['', 'none', 'null']:
                            # Build trinomial
                            subspecies_lower = subspecies.lower()
                            # Count under subspecies name too
                            name_to_count[subspecies_lower] += 1
                            unique_names.add(subspecies_lower)
                            
                            if cluster_ids:
                                for cid in cluster_ids:
                                    name_to_bins[subspecies_lower].add(cid)
                                    bin_to_names[cid].add(subspecies_lower)
                    
                    # Progress logging
                    if total_records % 500000 == 0:
                        logging.info(f"  Processed {total_records:,} records...")
            
            elapsed = time.time() - start_time
            logging.info(f"Index building complete in {elapsed:.1f} seconds")
            logging.info(f"  Total records: {total_records:,}")
            logging.info(f"  Records with cluster ID: {records_with_cluster:,}")
            logging.info(f"  Unique species names: {len(unique_names):,}")
            logging.info(f"  Unique BIN/OTU IDs: {len(unique_bins):,}")
            
            return dict(name_to_count), dict(name_to_bins), dict(bin_to_names), cluster_column
            
        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                # Reset
                name_to_count = defaultdict(int)
                name_to_bins = defaultdict(set)
                bin_to_names = defaultdict(set)
                total_records = 0
                records_with_cluster = 0
                unique_names = set()
                unique_bins = set()
                continue
            raise
    
    logging.error("Failed to read records file")
    sys.exit(1)


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def analyze_taxon(
    taxon: Taxon,
    name_to_count: Dict[str, int],
    name_to_bins: Dict[str, Set[str]],
    bin_to_names: Dict[str, Set[str]]
) -> TaxonResult:
    """
    Analyze a single taxon and determine its BAGS grade and status.
    
    Args:
        taxon: The taxon to analyze
        name_to_count: Mapping of species name -> record count
        name_to_bins: Mapping of species name -> set of BIN/OTU IDs
        bin_to_names: Mapping of BIN/OTU ID -> set of species names
    
    Returns:
        TaxonResult with grade, status, and other metrics
    """
    result = TaxonResult(taxon=taxon)
    taxon_names = taxon.all_names  # lowercase set
    valid_name_lower = taxon.valid_name.lower()
    
    # Step 1: Count records and find which names are recorded
    for name in taxon_names:
        count = name_to_count.get(name, 0)
        if count > 0:
            result.number_records += count
            result.names_recorded.add(name)
    
    # No records = Grade F, Status BLACK
    if result.number_records == 0:
        result.bags_grade = 'F'
        result.species_status = 'BLACK'
        return result
    
    # Step 2: Collect all BIN/OTUs for this taxon
    for name in taxon_names:
        bins = name_to_bins.get(name, set())
        result.bins_found.update(bins)
    
    # Step 3: Find all names in those BIN/OTUs
    all_names_in_bins: Set[str] = set()
    for bin_id in result.bins_found:
        names_in_bin = bin_to_names.get(bin_id, set())
        all_names_in_bins.update(names_in_bin)
    
    # Identify "other names" (not in focal taxon's list)
    other_names = all_names_in_bins - taxon_names
    result.other_names = sorted(other_names)
    
    # Step 4: Determine Status
    if other_names:
        # Taxonomic conflict - external names share BIN/OTU
        result.species_status = 'RED'
        result.bags_grade = 'E'
        return result
    
    # No external names - check nomenclatural status
    valid_recorded = valid_name_lower in result.names_recorded
    synonym_recorded = bool(result.names_recorded - {valid_name_lower})
    
    if valid_recorded and synonym_recorded:
        # Both valid and synonym(s) recorded - nomenclatural mess
        result.species_status = 'AMBER'
    elif valid_recorded:
        # Only valid name recorded - clean
        result.species_status = 'GREEN'
    else:
        # Only synonym(s) recorded - valid name absent
        result.species_status = 'BLUE'
    
    # Step 5: Determine Grade based on BIN/OTU count and records
    num_bins = len(result.bins_found)
    
    if num_bins == 0:
        # Records exist but no cluster assignment
        result.bags_grade = 'F'
    elif num_bins == 1:
        # Single cluster
        if result.number_records >= 11:
            result.bags_grade = 'A'
        elif result.number_records >= 3:
            result.bags_grade = 'B'
        else:
            result.bags_grade = 'D'
    else:
        # Multiple clusters (split)
        result.bags_grade = 'C'
    
    return result


def analyze_taxa_batch(
    taxa_batch: List[Taxon],
    name_to_count: Dict[str, int],
    name_to_bins: Dict[str, Set[str]],
    bin_to_names: Dict[str, Set[str]]
) -> List[TaxonResult]:
    """Analyze a batch of taxa (for parallel processing)."""
    return [analyze_taxon(t, name_to_count, name_to_bins, bin_to_names) for t in taxa_batch]


def analyze_all_taxa_parallel(
    taxa: List[Taxon],
    name_to_count: Dict[str, int],
    name_to_bins: Dict[str, Set[str]],
    bin_to_names: Dict[str, Set[str]],
    num_workers: int = None,
    batch_size: int = 1000
) -> List[TaxonResult]:
    """
    Analyze all taxa using parallel processing.
    
    Args:
        taxa: List of taxa to analyze
        name_to_count: Mapping of species name -> record count
        name_to_bins: Mapping of species name -> set of BIN/OTU IDs
        bin_to_names: Mapping of BIN/OTU ID -> set of species names
        num_workers: Number of parallel workers (default: CPU count)
        batch_size: Taxa per batch
    
    Returns:
        List of TaxonResult in same order as input taxa
    """
    if num_workers is None:
        num_workers = mp.cpu_count()
    
    logging.info(f"Analyzing {len(taxa):,} taxa using {num_workers} workers")
    start_time = time.time()
    
    # For small datasets or single worker, don't use multiprocessing
    if len(taxa) < batch_size or num_workers == 1:
        logging.info("Using single-threaded analysis")
        results = [analyze_taxon(t, name_to_count, name_to_bins, bin_to_names) for t in taxa]
        elapsed = time.time() - start_time
        logging.info(f"Analysis complete in {elapsed:.1f} seconds")
        return results
    
    # Split into batches
    batches = []
    for i in range(0, len(taxa), batch_size):
        batches.append(taxa[i:i + batch_size])
    
    logging.info(f"Split into {len(batches)} batches of up to {batch_size} taxa each")
    
    # Process in parallel
    # Note: For true HPC with shared memory, consider using multiprocessing.Manager
    # or memory-mapped structures. For simplicity, we pass the dicts to each worker.
    results = []
    completed = 0
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(analyze_taxa_batch, batch, name_to_count, name_to_bins, bin_to_names): i
            for i, batch in enumerate(batches)
        }
        
        # Collect results in order
        batch_results = [None] * len(batches)
        for future in as_completed(futures):
            batch_idx = futures[future]
            batch_results[batch_idx] = future.result()
            completed += 1
            if completed % 10 == 0 or completed == len(batches):
                logging.info(f"  Completed {completed}/{len(batches)} batches")
    
    # Flatten results
    for batch_result in batch_results:
        results.extend(batch_result)
    
    elapsed = time.time() - start_time
    logging.info(f"Analysis complete in {elapsed:.1f} seconds")
    
    return results


def analyze_all_taxa_serial(
    taxa: List[Taxon],
    name_to_count: Dict[str, int],
    name_to_bins: Dict[str, Set[str]],
    bin_to_names: Dict[str, Set[str]]
) -> List[TaxonResult]:
    """
    Analyze all taxa using single-threaded processing.
    
    More memory efficient for smaller datasets or when multiprocessing
    overhead isn't worth it.
    """
    logging.info(f"Analyzing {len(taxa):,} taxa (single-threaded)")
    start_time = time.time()
    
    results = []
    for i, taxon in enumerate(taxa):
        results.append(analyze_taxon(taxon, name_to_count, name_to_bins, bin_to_names))
        
        if (i + 1) % 10000 == 0:
            logging.info(f"  Processed {i + 1:,}/{len(taxa):,} taxa")
    
    elapsed = time.time() - start_time
    logging.info(f"Analysis complete in {elapsed:.1f} seconds")
    
    return results


# =============================================================================
# OUTPUT
# =============================================================================

def format_species_name(name: str) -> str:
    """Format species name with capitalized genus."""
    parts = name.split()
    if len(parts) >= 2:
        return parts[0].capitalize() + ' ' + ' '.join(parts[1:])
    elif len(parts) == 1:
        return parts[0].capitalize()
    return name


def write_results(
    results: List[TaxonResult],
    output_file: Path,
    input_columns: List[str]
) -> None:
    """
    Write analysis results to TSV file.
    
    Preserves all input columns and adds analysis columns.
    """
    logging.info(f"Writing results to {output_file}")
    
    # Define output columns: input columns first, then analysis columns
    analysis_columns = ['number_records', 'bags_grade', 'species_status', 'other_names']
    
    # Build fieldnames - input columns (excluding duplicates with analysis) + analysis
    output_columns = list(input_columns)
    for col in analysis_columns:
        if col not in output_columns:
            output_columns.append(col)
    
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=output_columns, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            
            for result in results:
                row = dict(result.taxon.input_data)
                row['number_records'] = result.number_records
                row['bags_grade'] = result.bags_grade
                row['species_status'] = result.species_status
                row['other_names'] = ';'.join(format_species_name(n) for n in result.other_names)
                
                writer.writerow(row)
        
        logging.info(f"Successfully wrote {len(results):,} results")
        
    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


def print_summary(results: List[TaxonResult]) -> None:
    """Print summary statistics."""
    logging.info("=" * 60)
    logging.info("SUMMARY")
    logging.info("=" * 60)
    
    # Grade distribution
    grade_counts = defaultdict(int)
    for r in results:
        grade_counts[r.bags_grade] += 1
    
    logging.info("BAGS Grade Distribution:")
    for grade in ['A', 'B', 'C', 'D', 'E', 'F']:
        count = grade_counts.get(grade, 0)
        pct = 100 * count / len(results) if results else 0
        logging.info(f"  Grade {grade}: {count:,} ({pct:.1f}%)")
    
    # Status distribution
    status_counts = defaultdict(int)
    for r in results:
        status_counts[r.species_status] += 1
    
    logging.info("Status Distribution:")
    for status in ['GREEN', 'AMBER', 'BLUE', 'RED', 'BLACK']:
        count = status_counts.get(status, 0)
        pct = 100 * count / len(results) if results else 0
        emoji = {'GREEN': '🟢', 'AMBER': '🟡', 'BLUE': '🔵', 'RED': '🔴', 'BLACK': '⚫'}.get(status, '')
        logging.info(f"  {status} {emoji}: {count:,} ({pct:.1f}%)")
    
    # Record coverage
    with_records = sum(1 for r in results if r.number_records > 0)
    total_records = sum(r.number_records for r in results)
    logging.info(f"Coverage:")
    logging.info(f"  Taxa with records: {with_records:,}/{len(results):,} ({100*with_records/len(results):.1f}%)")
    logging.info(f"  Total records matched: {total_records:,}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='HPC-optimized gap analysis for DNA barcode library curation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python gap_analysis.py \\
      --species-list species.tsv \\
      --records results/result_output.tsv \\
      --output results/gap_analysis.tsv

  # With parallel processing options
  python gap_analysis.py \\
      --species-list species.tsv \\
      --records results/result_output.tsv \\
      --output results/gap_analysis.tsv \\
      --workers 16 \\
      --batch-size 2000

  # Single-threaded mode (lower memory)
  python gap_analysis.py \\
      --species-list species.tsv \\
      --records results/result_output.tsv \\
      --output results/gap_analysis.tsv \\
      --workers 1
        """
    )
    
    parser.add_argument(
        '--species-list',
        type=Path,
        required=True,
        help='Path to species list TSV file (must have species and synonyms columns)'
    )
    
    parser.add_argument(
        '--records',
        type=Path,
        required=True,
        help='Path to records TSV file (e.g., result_output.tsv)'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Path to output gap analysis TSV file'
    )
    
    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        help='Number of parallel workers (default: number of CPUs)'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=1000,
        help='Taxa per batch for parallel processing (default: 1000)'
    )
    
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    logging.info("=" * 60)
    logging.info("HPC Gap Analysis for DNA Barcode Library Curation")
    logging.info("=" * 60)
    
    # Validate input files
    if not args.species_list.exists():
        logging.error(f"Species list not found: {args.species_list}")
        sys.exit(1)
    
    if not args.records.exists():
        logging.error(f"Records file not found: {args.records}")
        sys.exit(1)
    
    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # Load species list
    taxa, input_columns = load_species_list(args.species_list)
    
    # Build indices from records file
    name_to_count, name_to_bins, bin_to_names, cluster_col = build_indices_from_records(args.records)
    
    # Analyze taxa
    if args.workers == 1:
        results = analyze_all_taxa_serial(taxa, name_to_count, name_to_bins, bin_to_names)
    else:
        results = analyze_all_taxa_parallel(
            taxa, name_to_count, name_to_bins, bin_to_names,
            num_workers=args.workers,
            batch_size=args.batch_size
        )
    
    # Write results
    write_results(results, args.output, input_columns)
    
    # Print summary
    print_summary(results)
    
    logging.info("=" * 60)
    logging.info("Gap analysis complete!")
    logging.info("=" * 60)


if __name__ == '__main__':
    main()
