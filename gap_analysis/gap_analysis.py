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
        """All names (valid + synonyms) as normalized lowercase set."""
        names = {normalize_species_name(self.valid_name)}
        names.update(normalize_species_name(s) for s in self.synonyms)
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
    bin_uris: Set[str] = field(default_factory=set)  # Distinct bin_uri values
    otu_ids: Set[str] = field(default_factory=set)  # Distinct otu_id values


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# BOLD-specific filter settings
@dataclass
class BoldFilterSettings:
    """Settings for filtering raw BOLD data during index building."""
    enabled: bool = False           # Master switch: use BOLD parsing mode
    filter_species: bool = True     # Only include records with valid species names
    filter_kingdom: bool = True     # Filter by kingdom
    kingdom_list: List[str] = field(default_factory=lambda: ['Animalia'])
    marker: Optional[str] = 'COI-5P'  # Filter by marker_code column


# Values to treat as empty/missing for cluster IDs
EMPTY_CLUSTER_VALUES = frozenset(['', 'none', 'null', 'na', 'n/a', '-', '.'])

# Values to treat as empty/missing for species/subspecies names
# BOLD data commonly uses "None" for unidentified specimens
EMPTY_SPECIES_VALUES = frozenset(['', 'none', 'null', 'na', 'n/a', '-', '.'])


def sanitize_field(field_value: str) -> str:
    """
    Sanitize a field value to prevent parsing issues caused by unescaped quotes.
    
    BOLD raw TSV data contains unescaped double quotes (e.g. "Syn.) and embedded
    newlines/carriage returns that break standard CSV parsers.
    
    Only needed when reading raw BOLD data (--bold mode).
    
    Args:
        field_value: Raw field value from TSV
        
    Returns:
        Sanitized field value with quotes removed and whitespace cleaned
    """
    if not field_value:
        return ''
    
    # Remove embedded newlines and carriage returns that cause row merging
    field_value = field_value.replace('\r', ' ').replace('\n', ' ')
    
    # Remove all double quotes to prevent CSV parsing issues
    field_value = field_value.replace('"', '')
    
    return field_value.strip()


def normalize_species_name(name: str) -> str:
    """
    Normalize a species name for matching.
    
    Handles differences between naming conventions:
    - BOLD/UKSI use spaces: "Genus species"
    - UNITE uses underscores: "Genus_species"
    
    Converts to lowercase and replaces underscores with spaces.
    
    Args:
        name: Raw species name
    
    Returns:
        Normalized lowercase name with spaces instead of underscores
    """
    return name.lower().replace('_', ' ')


def is_valid_cluster_id(cluster_id: str) -> bool:
    """Check if a cluster ID is valid (not empty or a none-like string)."""
    return cluster_id.lower() not in EMPTY_CLUSTER_VALUES


def is_valid_species_name(species_name: str) -> bool:
    """
    Check if a species name is valid (not empty or a none-like placeholder).
    
    BOLD data commonly has "None" in species/subspecies columns for 
    specimens not identified to species level. These should be ignored
    as they don't represent real taxonomic identifications.
    """
    return species_name.lower() not in EMPTY_SPECIES_VALUES


def parse_cluster_ids(field_value: str) -> List[str]:
    """
    Parse cluster IDs from a field value.
    
    Handles pipe-separated values and filters out empty/none values.
    
    Args:
        field_value: Raw field value (may be pipe-separated)
    
    Returns:
        List of valid cluster IDs
    """
    if not field_value:
        return []
    
    ids = [c.strip() for c in field_value.split('|')]
    return [cid for cid in ids if cid and is_valid_cluster_id(cid)]


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

def detect_species_column(fieldnames: List[str]) -> str:
    """
    Detect which column to use for the valid species name.
    
    Prefers 'taxon_name' over 'species' if both are present.
    
    Returns column name to use.
    """
    if 'taxon_name' in fieldnames:
        return 'taxon_name'
    elif 'species' in fieldnames:
        return 'species'
    else:
        return None


def load_species_list(species_file: Path) -> Tuple[List[Taxon], List[str]]:
    """
    Load species list from TSV file.
    
    Expected columns:
    - 'taxon_name' or 'species': Valid name (required, prefers taxon_name)
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
                
                # Detect species column (prefer taxon_name over species)
                species_column = detect_species_column(input_columns)
                if species_column is None:
                    logging.error("Species list must have a 'taxon_name' or 'species' column")
                    sys.exit(1)
                
                logging.info(f"Using '{species_column}' column for valid species names")
                
                for row_idx, row in enumerate(reader):
                    valid_name = row.get(species_column, '').strip()
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


def detect_species_column_in_records(fieldnames: List[str]) -> str:
    """
    Detect which column to use for species name in records file.
    
    Prefers 'species' over 'organism' (NCBI files use 'organism').
    Note: 'ORGANISM_KEY' in UKSI files is NOT a species name column.
    
    Returns column name to use, or None if not found.
    """
    # Prefer 'species' (BOLD format)
    if 'species' in fieldnames:
        return 'species'
    # Fall back to 'organism' (NCBI format) - exact match only
    elif 'organism' in fieldnames:
        return 'organism'
    else:
        return None


def build_indices_from_records(
    records_file: Path,
    bold_filters: Optional[BoldFilterSettings] = None,
    chunk_size: int = 100000
) -> Tuple[Dict[str, int], Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]], str]:
    """
    Build indices from records file in a single pass.
    
    When bold_filters is provided and enabled, applies BOLD-specific parsing
    (quoting=QUOTE_NONE, field sanitisation) and row-level filters for
    species validity, kingdom, and marker code.
    
    Returns:
        - name_to_count: species_name (lowercase) -> record count
        - name_to_bins: species_name (lowercase) -> set of BIN/OTU IDs (for analysis)
        - bin_to_names: BIN/OTU ID -> set of species_names (lowercase)
        - name_to_bin_uris: species_name (lowercase) -> set of bin_uri values
        - name_to_otu_ids: species_name (lowercase) -> set of otu_id values
        - cluster_column: name of column used for primary clustering
    """
    logging.info(f"Building indices from {records_file}")
    start_time = time.time()
    
    name_to_count: Dict[str, int] = defaultdict(int)
    name_to_bins: Dict[str, Set[str]] = defaultdict(set)
    bin_to_names: Dict[str, Set[str]] = defaultdict(set)
    name_to_bin_uris: Dict[str, Set[str]] = defaultdict(set)
    name_to_otu_ids: Dict[str, Set[str]] = defaultdict(set)
    
    total_records = 0
    records_with_cluster = 0
    unique_names = set()
    unique_bins = set()
    
    cluster_column = None
    has_subspecies = False
    has_bin_uri = False
    has_otu_id = False
    
    # BOLD-specific filter tracking
    bold_mode = bold_filters is not None and bold_filters.enabled
    skipped_species = 0
    skipped_kingdom = 0
    skipped_marker = 0
    
    if bold_mode:
        logging.info("BOLD mode enabled - applying raw BOLD parsing and filters")
        if bold_filters.filter_species:
            logging.info("  Species filter: ON (skipping empty/None species)")
        if bold_filters.filter_kingdom:
            kingdom_set = frozenset(k.lower() for k in bold_filters.kingdom_list)
            logging.info(f"  Kingdom filter: ON (keeping: {', '.join(bold_filters.kingdom_list)})")
        if bold_filters.marker:
            logging.info(f"  Marker filter: ON (keeping: {bold_filters.marker})")
    
    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(records_file, 'r', encoding=encoding) as f:
                # Use QUOTE_NONE for raw BOLD data to handle unescaped quotes
                if bold_mode:
                    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
                else:
                    reader = csv.DictReader(f, delimiter='\t')
                fieldnames = reader.fieldnames or []
                
                # Detect cluster column (file-level decision)
                cluster_column = detect_cluster_column(fieldnames)
                logging.info(f"Using cluster column for analysis: {cluster_column}")
                
                # Check which cluster columns exist
                has_bin_uri = 'bin_uri' in fieldnames
                has_otu_id = 'otu_id' in fieldnames or 'OTU_ID' in fieldnames
                otu_col = 'otu_id' if 'otu_id' in fieldnames else 'OTU_ID' if 'OTU_ID' in fieldnames else None
                
                if has_bin_uri:
                    logging.info("bin_uri column detected")
                if has_otu_id:
                    logging.info(f"otu_id column detected ({otu_col})")
                
                # Detect species column (species or organism)
                species_column = detect_species_column_in_records(fieldnames)
                if species_column is None:
                    logging.error("Records file must have 'species' or 'organism' column")
                    sys.exit(1)
                logging.info(f"Using '{species_column}' column for species names")
                
                # Check for subspecies column
                has_subspecies = 'subspecies' in fieldnames
                if has_subspecies:
                    logging.info("Subspecies column detected - will match against both species and subspecies")
                
                # Detect BOLD-specific columns when in BOLD mode
                has_marker_col = 'marker_code' in fieldnames
                has_kingdom_col = 'kingdom' in fieldnames
                if bold_mode:
                    if bold_filters.marker and not has_marker_col:
                        logging.warning("Marker filter requested but 'marker_code' column not found in records file")
                    if bold_filters.filter_kingdom and not has_kingdom_col:
                        logging.warning("Kingdom filter requested but 'kingdom' column not found in records file")
                
                for row in reader:
                    total_records += 1
                    
                    # --- BOLD mode: sanitize fields to handle broken lines/quotes ---
                    if bold_mode:
                        row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                    
                    # --- BOLD filter: marker_code ---
                    if bold_mode and bold_filters.marker and has_marker_col:
                        marker_val = (row.get('marker_code', '') or '').strip()
                        if marker_val != bold_filters.marker:
                            skipped_marker += 1
                            continue
                    
                    # --- BOLD filter: kingdom ---
                    if bold_mode and bold_filters.filter_kingdom and has_kingdom_col:
                        kingdom_val = (row.get('kingdom', '') or '').strip().lower()
                        if not kingdom_val or kingdom_val not in kingdom_set:
                            skipped_kingdom += 1
                            continue
                    
                    # Get species name - skip if empty or placeholder like "None"
                    species = (row.get(species_column, '') or '').strip()
                    if not species or not is_valid_species_name(species):
                        # In BOLD mode with species filter, track this separately
                        if bold_mode and bold_filters.filter_species:
                            skipped_species += 1
                        continue
                    
                    # Get cluster IDs for primary analysis (may be pipe-separated)
                    cluster_field = (row.get(cluster_column, '') or '').strip()
                    cluster_ids = parse_cluster_ids(cluster_field)
                    
                    # Get bin_uri and otu_id separately for output columns
                    bin_uri_ids = []
                    otu_id_ids = []
                    if has_bin_uri:
                        bin_uri_ids = parse_cluster_ids((row.get('bin_uri', '') or '').strip())
                    if has_otu_id:
                        otu_id_ids = parse_cluster_ids((row.get(otu_col, '') or '').strip())
                    
                    # Process species name
                    species_lower = normalize_species_name(species)
                    name_to_count[species_lower] += 1
                    unique_names.add(species_lower)
                    
                    # Track bin_uri and otu_id separately
                    for bid in bin_uri_ids:
                        name_to_bin_uris[species_lower].add(bid)
                    for oid in otu_id_ids:
                        name_to_otu_ids[species_lower].add(oid)
                    
                    if cluster_ids:
                        records_with_cluster += 1
                        for cid in cluster_ids:
                            name_to_bins[species_lower].add(cid)
                            bin_to_names[cid].add(species_lower)
                            unique_bins.add(cid)
                    
                    # Also process subspecies if present
                    if has_subspecies:
                        subspecies = (row.get('subspecies', '') or '').strip()
                        if subspecies and is_valid_species_name(subspecies):
                            subspecies_lower = normalize_species_name(subspecies)
                            name_to_count[subspecies_lower] += 1
                            unique_names.add(subspecies_lower)
                            
                            # Track bin_uri and otu_id for subspecies too
                            for bid in bin_uri_ids:
                                name_to_bin_uris[subspecies_lower].add(bid)
                            for oid in otu_id_ids:
                                name_to_otu_ids[subspecies_lower].add(oid)
                            
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
            logging.info(f"  Records with valid cluster ID: {records_with_cluster:,}")
            logging.info(f"  Unique species names: {len(unique_names):,}")
            logging.info(f"  Unique BIN/OTU IDs: {len(unique_bins):,}")
            
            # BOLD filter summary
            if bold_mode:
                total_skipped = skipped_marker + skipped_kingdom + skipped_species
                records_kept = total_records - total_skipped
                logging.info(f"  BOLD filtering summary:")
                logging.info(f"    Records kept: {records_kept:,} of {total_records:,}")
                if bold_filters.marker:
                    logging.info(f"    Skipped (wrong marker): {skipped_marker:,}")
                if bold_filters.filter_kingdom:
                    logging.info(f"    Skipped (wrong kingdom): {skipped_kingdom:,}")
                if bold_filters.filter_species:
                    logging.info(f"    Skipped (no valid species): {skipped_species:,}")
            
            return (dict(name_to_count), dict(name_to_bins), dict(bin_to_names), 
                    dict(name_to_bin_uris), dict(name_to_otu_ids), cluster_column)
            
        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                # Reset
                name_to_count = defaultdict(int)
                name_to_bins = defaultdict(set)
                bin_to_names = defaultdict(set)
                name_to_bin_uris = defaultdict(set)
                name_to_otu_ids = defaultdict(set)
                total_records = 0
                records_with_cluster = 0
                unique_names = set()
                unique_bins = set()
                skipped_species = 0
                skipped_kingdom = 0
                skipped_marker = 0
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
    bin_to_names: Dict[str, Set[str]],
    name_to_bin_uris: Dict[str, Set[str]],
    name_to_otu_ids: Dict[str, Set[str]]
) -> TaxonResult:
    """
    Analyze a single taxon and determine its BAGS grade and status.
    
    Args:
        taxon: The taxon to analyze
        name_to_count: Mapping of species name -> record count
        name_to_bins: Mapping of species name -> set of BIN/OTU IDs
        bin_to_names: Mapping of BIN/OTU ID -> set of species names
        name_to_bin_uris: Mapping of species name -> set of bin_uri values
        name_to_otu_ids: Mapping of species name -> set of otu_id values
    
    Returns:
        TaxonResult with grade, status, and other metrics
    """
    result = TaxonResult(taxon=taxon)
    taxon_names = taxon.all_names  # normalized lowercase set
    valid_name_lower = normalize_species_name(taxon.valid_name)
    
    # Step 1: Count records and find which names are recorded
    for name in taxon_names:
        count = name_to_count.get(name, 0)
        if count > 0:
            result.number_records += count
            result.names_recorded.add(name)
        
        # Collect bin_uris and otu_ids for this taxon
        result.bin_uris.update(name_to_bin_uris.get(name, set()))
        result.otu_ids.update(name_to_otu_ids.get(name, set()))
    
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
    bin_to_names: Dict[str, Set[str]],
    name_to_bin_uris: Dict[str, Set[str]],
    name_to_otu_ids: Dict[str, Set[str]]
) -> List[TaxonResult]:
    """Analyze a batch of taxa (for parallel processing)."""
    return [analyze_taxon(t, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids) for t in taxa_batch]


def analyze_all_taxa_parallel(
    taxa: List[Taxon],
    name_to_count: Dict[str, int],
    name_to_bins: Dict[str, Set[str]],
    bin_to_names: Dict[str, Set[str]],
    name_to_bin_uris: Dict[str, Set[str]],
    name_to_otu_ids: Dict[str, Set[str]],
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
        name_to_bin_uris: Mapping of species name -> set of bin_uri values
        name_to_otu_ids: Mapping of species name -> set of otu_id values
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
        results = [analyze_taxon(t, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids) for t in taxa]
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
            executor.submit(analyze_taxa_batch, batch, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids): i
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
    bin_to_names: Dict[str, Set[str]],
    name_to_bin_uris: Dict[str, Set[str]],
    name_to_otu_ids: Dict[str, Set[str]]
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
        results.append(analyze_taxon(taxon, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids))
        
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
    analysis_columns = ['number_records', 'bags_grade', 'species_status', 'bin_uris', 'otu_ids', 'other_names']
    
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
                row['bin_uris'] = ';'.join(sorted(result.bin_uris))
                row['otu_ids'] = ';'.join(sorted(result.otu_ids))
                row['other_names'] = ';'.join(format_species_name(n) for n in result.other_names)
                
                writer.writerow(row)
        
        logging.info(f"Successfully wrote {len(results):,} results")
        
    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


def collect_relevant_names(
    taxa: List[Taxon],
    results: List[TaxonResult]
) -> Set[str]:
    """
    Collect all normalized species names relevant to the analysis.

    This includes:
    - All names (valid + synonyms) from the input species list
    - All "other names" discovered via BIN/OTU sharing during analysis

    Returns:
        Set of normalized (lowercase) species names
    """
    relevant = set()

    # All names from the input species list
    for taxon in taxa:
        relevant.update(taxon.all_names)

    # Other names found via BIN/OTU sharing
    for result in results:
        for name in result.other_names:
            relevant.add(normalize_species_name(name))

    return relevant


def write_filtered_records(
    records_file: Path,
    output_file: Path,
    relevant_names: Set[str],
    bold_filters: Optional[BoldFilterSettings] = None
) -> None:
    """
    Write a filtered copy of the records file containing only records
    whose species (or subspecies) name matches the set of relevant names.

    These are the records the gap analysis is based on: species list matches
    plus records from other taxa that share BIN/OTUs with listed species.

    All original columns are preserved.

    Args:
        records_file: Path to the input records TSV
        output_file: Path for the filtered output TSV
        relevant_names: Set of normalized species names to keep
        bold_filters: Optional BOLD filter settings (same as used for indexing)
    """
    logging.info(f"Writing filtered records to {output_file}")
    start_time = time.time()

    bold_mode = bold_filters is not None and bold_filters.enabled
    kingdom_set = None
    if bold_mode and bold_filters.filter_kingdom:
        kingdom_set = frozenset(k.lower() for k in bold_filters.kingdom_list)

    total_records = 0
    written_records = 0

    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(records_file, 'r', encoding=encoding) as fin, \
                 open(output_file, 'w', encoding='utf-8', newline='') as fout:

                if bold_mode:
                    reader = csv.DictReader(fin, delimiter='\t', quoting=csv.QUOTE_NONE)
                else:
                    reader = csv.DictReader(fin, delimiter='\t')

                fieldnames = reader.fieldnames or []

                species_column = detect_species_column_in_records(fieldnames)
                if species_column is None:
                    logging.error("Records file must have 'species' or 'organism' column")
                    sys.exit(1)

                has_subspecies = 'subspecies' in fieldnames
                has_marker_col = 'marker_code' in fieldnames
                has_kingdom_col = 'kingdom' in fieldnames

                writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter='\t',
                                        extrasaction='ignore')
                writer.writeheader()

                for row in reader:
                    total_records += 1

                    # Apply same BOLD sanitization and filters as index building
                    if bold_mode:
                        row = {k: sanitize_field(v) if v else '' for k, v in row.items()}

                    if bold_mode and bold_filters.marker and has_marker_col:
                        marker_val = (row.get('marker_code', '') or '').strip()
                        if marker_val != bold_filters.marker:
                            continue

                    if bold_mode and bold_filters.filter_kingdom and has_kingdom_col:
                        kingdom_val = (row.get('kingdom', '') or '').strip().lower()
                        if not kingdom_val or kingdom_val not in kingdom_set:
                            continue

                    species = (row.get(species_column, '') or '').strip()
                    if not species or not is_valid_species_name(species):
                        if bold_mode and bold_filters.filter_species:
                            continue
                        # Even outside BOLD mode, we skip empty species names
                        if not species:
                            continue

                    species_lower = normalize_species_name(species)

                    # Check if species or subspecies matches relevant names
                    match = species_lower in relevant_names
                    if not match and has_subspecies:
                        subspecies = (row.get('subspecies', '') or '').strip()
                        if subspecies and is_valid_species_name(subspecies):
                            match = normalize_species_name(subspecies) in relevant_names

                    if match:
                        writer.writerow(row)
                        written_records += 1

                    if total_records % 500000 == 0:
                        logging.info(f"  Scanned {total_records:,} records...")

            elapsed = time.time() - start_time
            logging.info(f"Filtered records complete in {elapsed:.1f} seconds")
            logging.info(f"  Total records scanned: {total_records:,}")
            logging.info(f"  Records written: {written_records:,}")
            return

        except UnicodeDecodeError:
            if encoding == 'utf-8':
                logging.warning("UTF-8 decoding failed, trying Latin-1")
                total_records = 0
                written_records = 0
                continue
            raise

    logging.error("Failed to read records file for filtering")
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

  # Raw BOLD data with default filters (Animalia, COI-5P, valid species only)
  python gap_analysis.py \\
      --species-list species.tsv \\
      --records BOLD_Public.tsv \\
      --output results/gap_analysis.tsv \\
      --bold

  # BOLD data with custom filters
  python gap_analysis.py \\
      --species-list species.tsv \\
      --records BOLD_Public.tsv \\
      --output results/gap_analysis.tsv \\
      --bold --marker matK --kingdom-list Plantae Fungi --no-filter-species
        """
    )
    
    parser.add_argument(
        '--species-list',
        type=Path,
        required=True,
        help='Path to species list TSV file (must have taxon_name or species column, plus optional synonyms column)'
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
        '--filtered-records',
        type=Path,
        default=None,
        help='Path to output filtered records TSV file (all records the analysis is based on). '
             'If not specified, derived from --output as <stem>_filtered_records.tsv'
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
    
    # BOLD-specific arguments
    bold_group = parser.add_argument_group(
        'BOLD mode',
        'Options for processing raw BOLD TSV data. Use --bold to enable. '
        'When --bold is set, default filters are applied (Animalia, COI-5P, valid species). '
        'Override individual filters as needed.'
    )
    
    bold_group.add_argument(
        '--bold',
        action='store_true',
        default=False,
        help='Enable BOLD mode: use QUOTE_NONE parsing, field sanitisation, and default filters'
    )
    
    bold_group.add_argument(
        '--filter-species',
        dest='filter_species',
        action='store_true',
        default=None,
        help='Only include records with valid species names (default: ON in BOLD mode, OFF otherwise)'
    )
    bold_group.add_argument(
        '--no-filter-species',
        dest='filter_species',
        action='store_false',
        help='Disable species name filtering'
    )
    
    bold_group.add_argument(
        '--filter-kingdom',
        dest='filter_kingdom',
        action='store_true',
        default=None,
        help='Filter records by kingdom (default: ON in BOLD mode, OFF otherwise)'
    )
    bold_group.add_argument(
        '--no-filter-kingdom',
        dest='filter_kingdom',
        action='store_false',
        help='Disable kingdom filtering'
    )
    
    bold_group.add_argument(
        '--kingdom-list',
        nargs='+',
        default=None,
        help='Kingdom(s) to include when --filter-kingdom is active (default: Animalia)'
    )
    
    bold_group.add_argument(
        '--marker',
        default=None,
        help='Marker code to filter by, e.g. "COI-5P" (default: COI-5P in BOLD mode, disabled otherwise)'
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
    
    # Build BOLD filter settings
    bold_filters = None
    if args.bold:
        bold_filters = BoldFilterSettings(
            enabled=True,
            filter_species=args.filter_species if args.filter_species is not None else True,
            filter_kingdom=args.filter_kingdom if args.filter_kingdom is not None else True,
            kingdom_list=args.kingdom_list if args.kingdom_list is not None else ['Animalia'],
            marker=args.marker if args.marker is not None else 'COI-5P'
        )
    else:
        # Non-BOLD mode: only apply explicit overrides
        if any(v is not None for v in [args.filter_species, args.filter_kingdom, args.kingdom_list, args.marker]):
            bold_filters = BoldFilterSettings(
                enabled=True,  # Enable parsing safeguards if any filter is set
                filter_species=args.filter_species if args.filter_species is not None else False,
                filter_kingdom=args.filter_kingdom if args.filter_kingdom is not None else False,
                kingdom_list=args.kingdom_list if args.kingdom_list is not None else ['Animalia'],
                marker=args.marker  # None = disabled
            )
    
    # Load species list
    taxa, input_columns = load_species_list(args.species_list)
    
    # Build indices from records file
    name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids, cluster_col = build_indices_from_records(
        args.records, bold_filters=bold_filters
    )
    
    # Analyze taxa
    if args.workers == 1:
        results = analyze_all_taxa_serial(taxa, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids)
    else:
        results = analyze_all_taxa_parallel(
            taxa, name_to_count, name_to_bins, bin_to_names, name_to_bin_uris, name_to_otu_ids,
            num_workers=args.workers,
            batch_size=args.batch_size
        )
    
    # Write results
    write_results(results, args.output, input_columns)

    # Write filtered records (second pass through records file)
    filtered_output = args.filtered_records
    if filtered_output is None:
        filtered_output = args.output.parent / f"{args.output.stem}_filtered_records.tsv"
    filtered_output.parent.mkdir(parents=True, exist_ok=True)
    relevant_names = collect_relevant_names(taxa, results)
    logging.info(f"Relevant species names for filtering: {len(relevant_names):,}")
    write_filtered_records(args.records, filtered_output, relevant_names, bold_filters=bold_filters)

    # Print summary
    print_summary(results)
    
    logging.info("=" * 60)
    logging.info("Gap analysis complete!")
    logging.info("=" * 60)


if __name__ == '__main__':
    main()
