#!/usr/bin/env python3
"""
BAGS Assessment Script

Calculates Barcode, Audit & Grade System (BAGS) grades for taxonomic data
based on OTU clustering results. Assesses the quality of species-level
barcode reference libraries.

BAGS Grades:
    A = 1 species name, 1 unshared OTU, >=11 specimens
    B = 1 species name, 1 unshared OTU, 3-10 specimens
    C = 1 species name split across multiple OTUs, each OTU unshared
    D = 1 species name, 1 unshared OTU, <3 specimens
    E = Multiple species names share the same OTU (BIN sharing)
    F = No OTU_ID assigned for this species

Grade precedence: E > C > A > B > D > F

Author: Ben Price / Claude
Date: 2025-01-23
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
from dataclasses import dataclass, field

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


@dataclass
class SpeciesData:
    """Data structure to hold aggregated data for a species/taxid."""
    taxid: int
    species_name: str
    otu_ids: Set[str] = field(default_factory=set)
    record_count: int = 0
    records_per_otu: Dict[str, int] = field(default_factory=dict)


@dataclass 
class OTUData:
    """Data structure to hold aggregated data for an OTU."""
    otu_id: str
    species_names: Set[str] = field(default_factory=set)  # lowercase for comparison
    species_names_original: Dict[str, str] = field(default_factory=dict)  # lowercase -> original case
    record_count: int = 0


def find_column_index(headers: List[str], candidates: List[str], required: bool = True) -> Optional[int]:
    """
    Find the index of a column from a list of candidate names.
    
    Args:
        headers: List of column headers
        candidates: List of possible column names to search for (case-insensitive)
        required: If True, raise error when not found
        
    Returns:
        Column index or None if not found and not required
    """
    headers_lower = [h.lower() for h in headers]
    for candidate in candidates:
        if candidate.lower() in headers_lower:
            return headers_lower.index(candidate.lower())
    
    if required:
        raise ValueError(f"Required column not found. Tried: {', '.join(candidates)}. "
                        f"Available columns: {', '.join(headers)}")
    return None


def read_input_file(input_path: Path, verbose: bool = False) -> Tuple[List[str], List[List[str]], int, int, int, int]:
    """
    Read input TSV file and identify relevant column indices.
    
    Args:
        input_path: Path to input TSV file
        verbose: Print progress information
        
    Returns:
        Tuple of (headers, rows, species_idx, otu_idx, accession_idx, taxid_idx or -1)
    """
    headers = []
    rows = []
    
    with open(input_path, 'r', encoding='utf-8') as f:
        header_line = f.readline().rstrip('\n\r')
        headers = header_line.split('\t')
        
        # Find required columns
        species_idx = find_column_index(headers, ['species', 'organism'])
        otu_idx = find_column_index(headers, ['OTU_ID', 'otu_id'])
        accession_idx = find_column_index(headers, ['accession', 'processid'])
        
        # Find optional taxid column
        taxid_idx = find_column_index(headers, ['taxid'], required=False)
        
        if verbose:
            logger.info(f"Species column: '{headers[species_idx]}' (index {species_idx})")
            logger.info(f"OTU_ID column: '{headers[otu_idx]}' (index {otu_idx})")
            logger.info(f"Accession column: '{headers[accession_idx]}' (index {accession_idx})")
            if taxid_idx is not None:
                logger.info(f"Taxid column: '{headers[taxid_idx]}' (index {taxid_idx})")
            else:
                logger.info("No taxid column found - will generate taxids")
        
        # Read all data rows
        line_num = 1
        for line in f:
            line_num += 1
            row = line.rstrip('\n\r').split('\t')
            rows.append(row)
            
            if verbose and line_num % 100000 == 0:
                logger.info(f"Read {line_num} rows...")
    
    if verbose:
        logger.info(f"Total rows read: {len(rows)}")
    
    return headers, rows, species_idx, otu_idx, accession_idx, taxid_idx if taxid_idx is not None else -1


def generate_taxids(rows: List[List[str]], species_idx: int, taxid_idx: int,
                    verbose: bool = False) -> Tuple[Dict[str, int], bool]:
    """
    Generate or validate taxids for species.
    
    Args:
        rows: Data rows
        species_idx: Index of species column
        taxid_idx: Index of taxid column (-1 if not present)
        verbose: Print progress information
        
    Returns:
        Tuple of (species_to_taxid dict, taxid_column_existed)
    """
    species_to_taxid: Dict[str, int] = {}
    taxid_column_existed = taxid_idx >= 0
    next_taxid = 1
    
    # First pass: collect existing taxids if column exists
    if taxid_column_existed:
        for row in rows:
            if len(row) <= species_idx:
                continue
            
            species = row[species_idx].strip() if len(row) > species_idx else ''
            if not species:
                continue
            
            species_lower = species.lower()
            
            # Get existing taxid if available
            taxid_str = row[taxid_idx].strip() if len(row) > taxid_idx else ''
            if taxid_str:
                try:
                    taxid = int(taxid_str)
                    if species_lower not in species_to_taxid:
                        species_to_taxid[species_lower] = taxid
                    next_taxid = max(next_taxid, taxid + 1)
                except ValueError:
                    pass
    
    # Second pass: assign taxids to species without them
    missing_taxid_species = []
    for row in rows:
        if len(row) <= species_idx:
            continue
        
        species = row[species_idx].strip() if len(row) > species_idx else ''
        if not species:
            continue
        
        species_lower = species.lower()
        
        if species_lower not in species_to_taxid:
            if taxid_column_existed:
                # Log error - this shouldn't happen
                logger.error(f"Species '{species}' has no taxid in existing taxid column - generating one")
                missing_taxid_species.append(species)
            species_to_taxid[species_lower] = next_taxid
            next_taxid += 1
    
    if verbose:
        logger.info(f"Total unique species: {len(species_to_taxid)}")
        if missing_taxid_species:
            logger.warning(f"Generated taxids for {len(missing_taxid_species)} species with missing values")
    
    return species_to_taxid, taxid_column_existed


def aggregate_data(rows: List[List[str]], species_idx: int, otu_idx: int,
                   species_to_taxid: Dict[str, int], verbose: bool = False
                   ) -> Tuple[Dict[int, SpeciesData], Dict[str, OTUData]]:
    """
    Aggregate data by species (taxid) and by OTU.
    
    Args:
        rows: Data rows
        species_idx: Index of species column
        otu_idx: Index of OTU_ID column
        species_to_taxid: Mapping of species name (lowercase) to taxid
        verbose: Print progress information
        
    Returns:
        Tuple of (species_data dict keyed by taxid, otu_data dict keyed by OTU_ID)
    """
    species_data: Dict[int, SpeciesData] = {}
    otu_data: Dict[str, OTUData] = {}
    
    # Track original case species names
    species_original_case: Dict[str, str] = {}  # lowercase -> original
    
    for row in rows:
        # Get species name
        species = row[species_idx].strip() if len(row) > species_idx else ''
        if not species:
            continue  # Skip rows without species
        
        species_lower = species.lower()
        
        # Track original case (keep first encountered)
        if species_lower not in species_original_case:
            species_original_case[species_lower] = species
        
        # Get taxid
        taxid = species_to_taxid.get(species_lower)
        if taxid is None:
            continue  # Shouldn't happen, but be safe
        
        # Get OTU_ID
        otu_id = row[otu_idx].strip() if len(row) > otu_idx else ''
        
        # Initialize species data if needed
        if taxid not in species_data:
            species_data[taxid] = SpeciesData(
                taxid=taxid,
                species_name=species_original_case[species_lower]
            )
        
        # Update species data
        sd = species_data[taxid]
        sd.record_count += 1
        
        if otu_id:
            sd.otu_ids.add(otu_id)
            sd.records_per_otu[otu_id] = sd.records_per_otu.get(otu_id, 0) + 1
            
            # Initialize OTU data if needed
            if otu_id not in otu_data:
                otu_data[otu_id] = OTUData(otu_id=otu_id)
            
            # Update OTU data
            od = otu_data[otu_id]
            od.species_names.add(species_lower)
            od.species_names_original[species_lower] = species_original_case[species_lower]
            od.record_count += 1
    
    if verbose:
        logger.info(f"Species with data: {len(species_data)}")
        logger.info(f"Unique OTUs: {len(otu_data)}")
    
    return species_data, otu_data


def calculate_bags_grade(species_data: SpeciesData, otu_data: Dict[str, OTUData]) -> Tuple[str, List[str], List[List[str]]]:
    """
    Calculate BAGS grade for a species.
    
    Grade precedence: E > C > A > B > D > F
    
    Args:
        species_data: Aggregated data for the species
        otu_data: Dictionary of all OTU data
        
    Returns:
        Tuple of (grade, list of OTU_IDs, list of sharers per OTU)
    """
    species_lower = species_data.species_name.lower()
    otu_ids = sorted(species_data.otu_ids)
    
    # Grade F: No OTU_ID assigned
    if not otu_ids:
        return 'F', [], []
    
    # Check for sharing (Grade E condition)
    # Also collect sharers for each OTU
    is_shared = False
    sharers_per_otu: List[List[str]] = []
    
    for otu_id in otu_ids:
        od = otu_data.get(otu_id)
        if od is None:
            sharers_per_otu.append([])
            continue
        
        # Get other species sharing this OTU (not the current species)
        other_species = [
            od.species_names_original[sp] 
            for sp in od.species_names 
            if sp != species_lower
        ]
        
        if other_species:
            is_shared = True
        
        sharers_per_otu.append(sorted(other_species))
    
    # Grade E: BIN sharing with other species (highest priority after F check)
    if is_shared:
        return 'E', otu_ids, sharers_per_otu
    
    # Grade C: Species split across multiple unshared OTUs
    if len(otu_ids) > 1:
        return 'C', otu_ids, sharers_per_otu
    
    # Single unshared OTU - grade based on specimen count
    record_count = species_data.record_count
    
    # Grade A: >=11 specimens
    if record_count >= 11:
        return 'A', otu_ids, sharers_per_otu
    
    # Grade B: 3-10 specimens
    if record_count >= 3:
        return 'B', otu_ids, sharers_per_otu
    
    # Grade D: <3 specimens
    return 'D', otu_ids, sharers_per_otu


def write_taxid_output(output_path: Path, headers: List[str], rows: List[List[str]],
                       species_idx: int, taxid_idx: int, species_to_taxid: Dict[str, int],
                       taxid_existed: bool, verbose: bool = False) -> None:
    """
    Write output TSV with taxid column appended (if it didn't exist).
    
    Args:
        output_path: Path to output file
        headers: Original headers
        rows: Data rows
        species_idx: Index of species column
        taxid_idx: Index of taxid column (-1 if didn't exist)
        species_to_taxid: Mapping of species name (lowercase) to taxid
        taxid_existed: Whether taxid column existed in input
        verbose: Print progress information
    """
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header
        if taxid_existed:
            f.write('\t'.join(headers) + '\n')
        else:
            f.write('\t'.join(headers) + '\ttaxid\n')
        
        # Write data rows
        for row in rows:
            if taxid_existed:
                # Check if taxid needs to be filled in
                species = row[species_idx].strip() if len(row) > species_idx else ''
                current_taxid = row[taxid_idx].strip() if len(row) > taxid_idx else ''
                
                if species and not current_taxid:
                    # Fill in missing taxid
                    taxid = species_to_taxid.get(species.lower(), '')
                    # Ensure row has enough columns
                    while len(row) <= taxid_idx:
                        row.append('')
                    row[taxid_idx] = str(taxid) if taxid else ''
                
                f.write('\t'.join(row) + '\n')
            else:
                # Append taxid column
                species = row[species_idx].strip() if len(row) > species_idx else ''
                taxid = species_to_taxid.get(species.lower(), '') if species else ''
                f.write('\t'.join(row) + '\t' + (str(taxid) if taxid else '') + '\n')
    
    if verbose:
        logger.info(f"Wrote taxid output to: {output_path}")


def write_bags_output(output_path: Path, species_data: Dict[int, SpeciesData],
                      otu_data: Dict[str, OTUData], verbose: bool = False) -> Dict[str, int]:
    """
    Write BAGS assessment output TSV.
    
    Args:
        output_path: Path to output file
        species_data: Dictionary of species data keyed by taxid
        otu_data: Dictionary of OTU data keyed by OTU_ID
        verbose: Print progress information
        
    Returns:
        Dictionary of grade counts for statistics
    """
    grade_counts: Dict[str, int] = defaultdict(int)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header
        f.write('taxid\tBAGS\tOTU_ID\tsharers\n')
        
        # Process each species in taxid order
        for taxid in sorted(species_data.keys()):
            sd = species_data[taxid]
            
            # Calculate BAGS grade
            grade, otu_ids, sharers_per_otu = calculate_bags_grade(sd, otu_data)
            grade_counts[grade] += 1
            
            # Format OTU_ID column (pipe-separated)
            otu_id_str = '|'.join(otu_ids)
            
            # Format sharers column (comma-separated within OTU, pipe-separated between OTUs)
            sharers_parts = []
            for sharers in sharers_per_otu:
                sharers_parts.append(','.join(sharers))
            sharers_str = '|'.join(sharers_parts)
            
            f.write(f'{taxid}\t{grade}\t{otu_id_str}\t{sharers_str}\n')
    
    if verbose:
        logger.info(f"Wrote BAGS assessment to: {output_path}")
        logger.info("Grade distribution:")
        for grade in ['A', 'B', 'C', 'D', 'E', 'F']:
            count = grade_counts.get(grade, 0)
            logger.info(f"  Grade {grade}: {count}")
    
    return dict(grade_counts)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate BAGS (Barcode, Audit & Grade System) grades for taxonomic data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
BAGS Grades:
  A = 1 species name, 1 unshared OTU, >=11 specimens (high quality)
  B = 1 species name, 1 unshared OTU, 3-10 specimens (good quality)
  C = 1 species name split across multiple OTUs, each unshared (split species)
  D = 1 species name, 1 unshared OTU, <3 specimens (low coverage)
  E = Multiple species names share the same OTU (BIN sharing/conflict)
  F = No OTU_ID assigned (no barcode data)

Grade precedence (highest to lowest): E > C > A > B > D > F

Examples:
  # Basic usage
  python bags_assessment.py input.tsv
  
  # With verbose output
  python bags_assessment.py -v input.tsv
  
  # Specify output directory
  python bags_assessment.py input.tsv -o /path/to/output/

Input Format:
  Tab-separated file with at minimum:
    - 'species' or 'organism' column (species names)
    - 'OTU_ID' column (OTU cluster assignments)
    - 'accession' or 'processid' column (record identifiers)
    - Optional: 'taxid' column (if not present, will be generated)

Output Files:
  - assessed_BAGS.tsv: BAGS grades with columns taxid, BAGS, OTU_ID, sharers
  - {input}_taxid.tsv: Input data with taxid column appended (if not present)
'''
    )
    
    parser.add_argument('input', type=str,
                        help='Input TSV file with OTU clustering results')
    parser.add_argument('-o', '--output-dir', type=str, default=None,
                        help='Output directory (default: same as input file)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress and statistics')
    
    args = parser.parse_args()
    
    # Set up paths
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    output_dir = Path(args.output_dir) if args.output_dir else input_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    bags_output = output_dir / 'assessed_BAGS.tsv'
    taxid_output = output_dir / f'{input_path.stem}_taxid.tsv'
    
    if args.verbose:
        logger.info(f"Input file: {input_path}")
        logger.info(f"BAGS output: {bags_output}")
        logger.info(f"Taxid output: {taxid_output}")
    
    try:
        # Step 1: Read input file
        if args.verbose:
            logger.info("Reading input file...")
        headers, rows, species_idx, otu_idx, accession_idx, taxid_idx = read_input_file(
            input_path, args.verbose
        )
        
        # Step 2: Generate or validate taxids
        if args.verbose:
            logger.info("Processing taxids...")
        species_to_taxid, taxid_existed = generate_taxids(
            rows, species_idx, taxid_idx, args.verbose
        )
        
        # Step 3: Aggregate data by species and OTU
        if args.verbose:
            logger.info("Aggregating data...")
        species_data, otu_data = aggregate_data(
            rows, species_idx, otu_idx, species_to_taxid, args.verbose
        )
        
        # Step 4: Write taxid output (only if taxid column was missing or had gaps)
        if args.verbose:
            logger.info("Writing taxid output...")
        write_taxid_output(
            taxid_output, headers, rows, species_idx, taxid_idx, 
            species_to_taxid, taxid_existed, args.verbose
        )
        
        # Step 5: Calculate and write BAGS assessment
        if args.verbose:
            logger.info("Calculating BAGS grades...")
        grade_counts = write_bags_output(
            bags_output, species_data, otu_data, args.verbose
        )
        
        if args.verbose:
            total_species = sum(grade_counts.values())
            logger.info(f"BAGS assessment completed for {total_species} species")
        
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == '__main__':
    main()
