#!/usr/bin/env python3
"""
uksi_jncc_annotation_v2.py

Annotates UKSI valid species list with JNCC conservation designations.
Matches UKSI taxon_version_key against the semicolon-separated values in
the JNCC 'included_tvk_list' column.

Usage:
    python uksi_jncc_annotation_v2.py --input uksi_species.tsv --jncc jncc_designations.tsv [--output-dir /path/to/output]

Outputs:
    - uksi_valid_species_jncc_annotated.tsv
    - uksi_jncc_annotation.log (includes errors for unmatched JNCC TVKs)
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional


# JNCC designation columns in order
JNCC_DESIGNATION_COLUMNS = [
    "A: Bern Convention",
    "C: Birds Directive",
    "C1: Convention on Migratory Species",
    "C2: OSPAR",
    "D: Habitats Directive",
    "E: EC Cites",
    "F: Global Red list status",
    "Fa: Red Listing based on pre 1994 IUCN guidelines",
    "Fb: Red Listing based on 1994 IUCN guidelines",
    "Fc: Red listing based on 2001 IUCN guidelines",
    "Fd: Red data categories  - birds (not based on IUCN criteria)",
    "Fe: Red data categories - Spiders (not based on IUCN criteria)",
    "Ga: Rare and scarce species",
    "Gb: Rare and scarce species (not based on IUCN criteria)",
    "Ha: Biodiversity Action Plan UK list of priority species",
    "Hb: Biodiversity Lists - England",
    "Hc: Biodiversity Lists - Scotland",
    "Hd: Biodiversity Lists - Wales",
    "He: Biodiversity Lists - Northern Ireland",
    "I: Wildlife and Countryside Act 1981",
    "J: The Wildlife (Northern Ireland) Order 1985",
    "K: The Conservation of Habitats and Species Regulations 2010",
    "L: The Conservation (Nature Habitats, etc_) Regulations (NI) 199",
    "M: Protection of Badgers Act",
]

# JNCC file column containing the semicolon-separated TVK list
JNCC_TVK_LIST_COLUMN = "included_tvk_list"

# Column names for matching info in output
JNCC_MATCHING_TVK_COL = "jncc_matching_tvk"
JNCC_ROW_INDEX_COL = "jncc_row_index"


def setup_logging(log_file: Path) -> logging.Logger:
    """Configure logging to both file and stderr."""
    logger = logging.getLogger("uksi_jncc_annotation")
    logger.setLevel(logging.DEBUG)
    
    # File handler - captures everything including DEBUG
    fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    
    # Stderr handler (for console output)
    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(sh)
    
    return logger



def load_jncc_data(jncc_file: Path, logger: logging.Logger) -> Tuple[Dict[str, List[dict]], Set[str]]:
    """
    Load JNCC conservation designations indexed by individual TVK.
    
    The JNCC file has an 'included_tvk_list' column containing semicolon-separated TVKs.
    We create a mapping from each individual TVK to the row's designation data.
    
    Returns:
        Tuple of:
        - Dict mapping individual TVK -> list of dicts containing designation data and row info
          (list because same TVK might appear in multiple rows)
        - Set of all unique TVKs found in the JNCC file
    """
    tvk_to_designations = defaultdict(list)
    all_jncc_tvks = set()
    
    with open(jncc_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # Verify required columns exist
        if JNCC_TVK_LIST_COLUMN not in reader.fieldnames:
            raise ValueError(f"JNCC file missing required column: {JNCC_TVK_LIST_COLUMN}")
        
        missing_cols = [col for col in JNCC_DESIGNATION_COLUMNS if col not in reader.fieldnames]
        if missing_cols:
            logger.warning(f"JNCC file missing designation columns: {missing_cols}")
        
        for row_idx, row in enumerate(reader, start=2):  # start=2 because row 1 is header
            tvk_list_str = row.get(JNCC_TVK_LIST_COLUMN, "").strip()
            if not tvk_list_str:
                continue
            
            # Parse semicolon-separated TVKs
            tvks = [t.strip() for t in tvk_list_str.split(';') if t.strip()]
            
            # Extract designation values for this row
            designations = {}
            for col in JNCC_DESIGNATION_COLUMNS:
                designations[col] = row.get(col, "").strip()
            
            # Store row metadata
            row_data = {
                'designations': designations,
                'row_index': row_idx,
                'recommended_name': row.get('Recommended_taxon_name', '').strip(),
                'source_tvk_list': tvk_list_str
            }
            
            # Map each individual TVK to this row's data
            for tvk in tvks:
                tvk_to_designations[tvk].append(row_data)
                all_jncc_tvks.add(tvk)
    
    logger.info(f"Loaded {len(all_jncc_tvks)} unique TVKs from JNCC file (from {row_idx - 1} rows)")
    return dict(tvk_to_designations), all_jncc_tvks


def merge_designations(designation_list: List[Dict[str, str]]) -> Dict[str, str]:
    """
    Merge multiple sets of designation values, combining with semicolons where different.
    """
    if not designation_list:
        return {col: "" for col in JNCC_DESIGNATION_COLUMNS}
    
    if len(designation_list) == 1:
        return designation_list[0].copy()
    
    merged = {}
    for col in JNCC_DESIGNATION_COLUMNS:
        values = set()
        for d in designation_list:
            val = d.get(col, "")
            if val:
                # Split existing semicolon-separated values
                for part in val.split(';'):
                    part = part.strip()
                    if part:
                        values.add(part)
        merged[col] = ';'.join(sorted(values)) if values else ""
    
    return merged



def process_species(
    input_file: Path,
    tvk_to_designations: Dict[str, List[dict]],
    output_file: Path,
    logger: logging.Logger
) -> Set[str]:
    """
    Process UKSI species file and annotate with JNCC designations.
    
    Matches UKSI taxon_version_key against JNCC included_tvk_list entries.
    
    Returns:
        Set of JNCC TVKs that were matched to UKSI species
    """
    matched_jncc_tvks = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        input_fieldnames = list(reader.fieldnames)
        
        # Build output fieldnames: input columns + matching info + JNCC designation columns
        output_fieldnames = input_fieldnames + [
            JNCC_MATCHING_TVK_COL,
            JNCC_ROW_INDEX_COL
        ] + JNCC_DESIGNATION_COLUMNS
        
        with open(output_file, 'w', encoding='utf-8', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=output_fieldnames, delimiter='\t')
            writer.writeheader()
            
            species_count = 0
            matched_count = 0
            
            for row in reader:
                species_count += 1
                
                # Get the UKSI taxon_version_key
                uksi_tvk = row.get('taxon_version_key', '').strip()
                
                # Build output row starting with input data
                output_row = dict(row)
                
                # Check if this TVK exists in JNCC data
                if uksi_tvk and uksi_tvk in tvk_to_designations:
                    matched_count += 1
                    matched_jncc_tvks.add(uksi_tvk)
                    
                    # Get all matching row data (could be multiple JNCC rows)
                    matching_rows = tvk_to_designations[uksi_tvk]
                    
                    # Collect row indices for reference
                    row_indices = [str(r['row_index']) for r in matching_rows]
                    
                    # Merge designations from all matching rows
                    designation_dicts = [r['designations'] for r in matching_rows]
                    merged = merge_designations(designation_dicts)
                    
                    output_row[JNCC_MATCHING_TVK_COL] = uksi_tvk
                    output_row[JNCC_ROW_INDEX_COL] = ';'.join(row_indices)
                    
                    for col in JNCC_DESIGNATION_COLUMNS:
                        output_row[col] = merged.get(col, "")
                    
                    if len(matching_rows) > 1:
                        logger.debug(f"Multiple JNCC rows matched for TVK {uksi_tvk}: rows {row_indices}")
                else:
                    # No match - empty strings for all JNCC columns
                    output_row[JNCC_MATCHING_TVK_COL] = ""
                    output_row[JNCC_ROW_INDEX_COL] = ""
                    for col in JNCC_DESIGNATION_COLUMNS:
                        output_row[col] = ""
                
                writer.writerow(output_row)
    
    logger.info(f"Processed {species_count} UKSI species, {matched_count} matched to JNCC designations")
    return matched_jncc_tvks



def log_unmatched_jncc_tvks(
    all_jncc_tvks: Set[str],
    matched_tvks: Set[str],
    tvk_to_designations: Dict[str, List[dict]],
    logger: logging.Logger
) -> None:
    """
    Log JNCC TVKs that were not found in the UKSI file.
    
    These are TVKs from the included_tvk_list that don't match any UKSI taxon_version_key.
    """
    unmatched = all_jncc_tvks - matched_tvks
    
    if unmatched:
        logger.warning(f"JNCC TVKs not found in UKSI file: {len(unmatched)} of {len(all_jncc_tvks)}")
        logger.info("Logging unmatched JNCC TVKs to error log...")
        
        for tvk in sorted(unmatched):
            # Get info about where this TVK came from
            row_data_list = tvk_to_designations.get(tvk, [])
            for row_data in row_data_list:
                logger.error(
                    f"JNCC TVK not found in UKSI: {tvk} "
                    f"(from JNCC row {row_data['row_index']}, "
                    f"name: '{row_data['recommended_name']}')"
                )
    else:
        logger.info("All JNCC TVKs matched to UKSI species")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate UKSI valid species with JNCC conservation designations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python uksi_jncc_annotation_v2.py --input uksi_species.tsv --jncc jncc_designations.tsv
    
    # With output directory
    python uksi_jncc_annotation_v2.py --input uksi_species.tsv --jncc jncc_designations.tsv --output-dir /results

Matching logic:
    - Looks up UKSI taxon_version_key in JNCC included_tvk_list
    - included_tvk_list contains semicolon-separated TVKs
    - If multiple JNCC rows contain the same TVK, designations are merged
    - Unmatched JNCC TVKs are logged as errors
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        type=Path,
        help='Path to UKSI valid species TSV file'
    )
    
    parser.add_argument(
        '--jncc', '-j',
        required=True,
        type=Path,
        help='Path to JNCC conservation designations TSV file'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: same as input file)'
    )
    
    parser.add_argument(
        '--log-file',
        type=Path,
        default=None,
        help='Log file path (default: uksi_jncc_annotation.log in output directory)'
    )
    
    return parser.parse_args()



def main():
    """Main entry point."""
    args = parse_args()
    
    # Validate input files exist
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not args.jncc.exists():
        print(f"Error: JNCC file not found: {args.jncc}", file=sys.stderr)
        sys.exit(1)
    
    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = args.input.parent
    
    # Setup output paths
    output_file = output_dir / "uksi_valid_species_jncc_annotated.tsv"
    log_file = args.log_file if args.log_file else output_dir / "uksi_jncc_annotation.log"
    
    # Setup logging
    logger = setup_logging(log_file)
    logger.info(f"Starting JNCC annotation (v2 - using included_tvk_list)")
    logger.info(f"UKSI input file: {args.input}")
    logger.info(f"JNCC input file: {args.jncc}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"Log file: {log_file}")
    
    try:
        # Load JNCC data - creates mapping from individual TVKs to designations
        logger.info("Loading JNCC conservation designations...")
        tvk_to_designations, all_jncc_tvks = load_jncc_data(args.jncc, logger)
        
        # Process UKSI species and annotate with JNCC designations
        logger.info("Processing UKSI species file...")
        matched_tvks = process_species(
            args.input,
            tvk_to_designations,
            output_file,
            logger
        )
        
        # Log unmatched JNCC TVKs as errors
        log_unmatched_jncc_tvks(all_jncc_tvks, matched_tvks, tvk_to_designations, logger)
        
        logger.info(f"JNCC annotation complete. Output written to: {output_file}")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
