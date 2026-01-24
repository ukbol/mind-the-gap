#!/usr/bin/env python3
"""
uksi_jncc_annotation.py

Annotates UKSI valid species list with JNCC conservation designations.
Matches on TVK (taxon version key) for both valid names and synonyms.

Usage:
    python uksi_jncc_annotation.py --input species.tsv --jncc jncc.tsv [--output-dir /path/to/output] [--stdout]

Outputs:
    - uksi_valid_species_jncc_annotated.tsv (or stdout if --stdout)
    - uksi_jncc_annotation.log
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
]
JNCC_DESIGNATION_COLUMNS_CONTINUED = [
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

# Combine all designation columns
JNCC_DESIGNATION_COLUMNS = JNCC_DESIGNATION_COLUMNS + JNCC_DESIGNATION_COLUMNS_CONTINUED

# Column names for JNCC matching info
JNCC_MATCHING_TVK_COL = "jncc_matching_tvk"
JNCC_MATCHING_STATUS_COL = "jncc_matching_tvk_status"

# JNCC file TVK column
JNCC_TVK_COLUMN = "Recommended_taxon_version"


def setup_logging(log_file: Path) -> logging.Logger:
    """Configure logging to both file and stderr."""
    logger = logging.getLogger("uksi_jncc_annotation")
    logger.setLevel(logging.DEBUG)
    
    # File handler
    fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    
    # Stderr handler (for console output that won't interfere with stdout)
    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(sh)
    
    return logger


def load_jncc_data(jncc_file: Path, logger: logging.Logger) -> Dict[str, Dict[str, str]]:
    """
    Load JNCC conservation designations indexed by TVK.
    
    Returns:
        Dict mapping TVK -> dict of designation column -> value
    """
    jncc_data = {}
    
    with open(jncc_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # Verify required columns exist
        if JNCC_TVK_COLUMN not in reader.fieldnames:
            raise ValueError(f"JNCC file missing required column: {JNCC_TVK_COLUMN}")
        
        missing_cols = [col for col in JNCC_DESIGNATION_COLUMNS if col not in reader.fieldnames]
        if missing_cols:
            logger.warning(f"JNCC file missing designation columns: {missing_cols}")
        
        for row in reader:
            tvk = row.get(JNCC_TVK_COLUMN, "").strip()
            if not tvk:
                continue
            
            # Extract designation values
            designations = {}
            for col in JNCC_DESIGNATION_COLUMNS:
                designations[col] = row.get(col, "").strip()
            
            if tvk in jncc_data:
                logger.warning(f"Duplicate TVK in JNCC file: {tvk}")
                # Merge designations (combine with semicolon if different)
                for col in JNCC_DESIGNATION_COLUMNS:
                    existing = jncc_data[tvk][col]
                    new_val = designations[col]
                    if new_val and existing and new_val != existing:
                        jncc_data[tvk][col] = f"{existing};{new_val}"
                    elif new_val and not existing:
                        jncc_data[tvk][col] = new_val
            else:
                jncc_data[tvk] = designations
    
    logger.info(f"Loaded {len(jncc_data)} TVKs from JNCC file")
    return jncc_data


def merge_designations(existing: Dict[str, str], new: Dict[str, str]) -> Dict[str, str]:
    """
    Merge two sets of designation values, combining with semicolons where different.
    """
    merged = {}
    for col in JNCC_DESIGNATION_COLUMNS:
        existing_val = existing.get(col, "")
        new_val = new.get(col, "")
        
        if not existing_val:
            merged[col] = new_val
        elif not new_val:
            merged[col] = existing_val
        elif existing_val == new_val:
            merged[col] = existing_val
        else:
            # Combine unique values with semicolon
            existing_parts = set(existing_val.split(';'))
            new_parts = set(new_val.split(';'))
            all_parts = existing_parts | new_parts
            # Remove empty strings
            all_parts.discard('')
            merged[col] = ';'.join(sorted(all_parts))
    
    return merged


def find_jncc_match(
    valid_tvk: str,
    synonym_tvks: List[str],
    jncc_data: Dict[str, Dict[str, str]],
    logger: logging.Logger
) -> Tuple[Optional[Dict[str, str]], str, str]:
    """
    Find JNCC designation match for a species.
    
    Checks valid TVK first, then synonym TVKs.
    Merges designations if multiple matches found.
    
    Returns:
        Tuple of (merged_designations, matching_tvks, status)
        - merged_designations: dict of designation column -> value, or None if no match
        - matching_tvks: semicolon-separated TVKs that matched
        - status: 'valid', 'synonym', or 'valid;synonym' if both matched
    """
    merged_designations = None
    matching_tvks = []
    statuses = []
    
    # Check valid TVK first
    if valid_tvk and valid_tvk in jncc_data:
        merged_designations = jncc_data[valid_tvk].copy()
        matching_tvks.append(valid_tvk)
        statuses.append("valid")
    
    # Check synonym TVKs
    matched_synonyms = []
    for syn_tvk in synonym_tvks:
        if syn_tvk and syn_tvk in jncc_data:
            matched_synonyms.append(syn_tvk)
            if merged_designations is None:
                merged_designations = jncc_data[syn_tvk].copy()
            else:
                merged_designations = merge_designations(merged_designations, jncc_data[syn_tvk])
            matching_tvks.append(syn_tvk)
    
    if matched_synonyms:
        statuses.append("synonym")
    
    # Log if multiple matches (valid + synonyms or multiple synonyms)
    if len(matching_tvks) > 1:
        logger.info(f"Multiple JNCC matches for valid TVK {valid_tvk}: matched TVKs = {matching_tvks}")
    
    matching_tvk_str = ";".join(matching_tvks) if matching_tvks else ""
    status_str = ";".join(statuses) if statuses else ""
    
    return merged_designations, matching_tvk_str, status_str


def process_species(
    input_file: Path,
    jncc_data: Dict[str, Dict[str, str]],
    output_file: Optional[Path],
    use_stdout: bool,
    logger: logging.Logger
) -> Set[str]:
    """
    Process valid species file and annotate with JNCC designations.
    
    Returns:
        Set of matched JNCC TVKs (for tracking unmatched)
    """
    matched_jncc_tvks = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        input_fieldnames = reader.fieldnames
        
        # Build output fieldnames: input columns + JNCC columns
        output_fieldnames = list(input_fieldnames) + [
            JNCC_MATCHING_TVK_COL,
            JNCC_MATCHING_STATUS_COL
        ] + JNCC_DESIGNATION_COLUMNS
        
        # Setup output
        if use_stdout:
            out_handle = sys.stdout
        else:
            out_handle = open(output_file, 'w', encoding='utf-8', newline='')
        
        try:
            writer = csv.DictWriter(out_handle, fieldnames=output_fieldnames, delimiter='\t')
            writer.writeheader()
            
            species_count = 0
            matched_count = 0
            
            for row in reader:
                species_count += 1
                
                # Get TVKs to check
                valid_tvk = row.get('taxon_version_key', '').strip()
                synonym_tvk_str = row.get('synonym_tvk_list', '').strip()
                
                # Parse synonym TVKs (semicolon-separated)
                synonym_tvks = []
                if synonym_tvk_str:
                    synonym_tvks = [t.strip() for t in synonym_tvk_str.split(';') if t.strip()]
                
                # Find JNCC match
                designations, matching_tvk, status = find_jncc_match(
                    valid_tvk, synonym_tvks, jncc_data, logger
                )
                
                # Build output row
                output_row = dict(row)
                output_row[JNCC_MATCHING_TVK_COL] = matching_tvk
                output_row[JNCC_MATCHING_STATUS_COL] = status
                
                if designations:
                    matched_count += 1
                    for tvk in matching_tvk.split(';'):
                        matched_jncc_tvks.add(tvk)
                    for col in JNCC_DESIGNATION_COLUMNS:
                        output_row[col] = designations.get(col, "")
                else:
                    # No match - empty strings for all designation columns
                    for col in JNCC_DESIGNATION_COLUMNS:
                        output_row[col] = ""
                
                writer.writerow(output_row)
            
            logger.info(f"Processed {species_count} species, {matched_count} matched to JNCC designations")
            
        finally:
            if not use_stdout:
                out_handle.close()
    
    return matched_jncc_tvks


def log_unmatched_jncc(
    jncc_data: Dict[str, Dict[str, str]],
    matched_tvks: Set[str],
    logger: logging.Logger
) -> None:
    """Log JNCC TVKs that were not matched to any species."""
    unmatched = set(jncc_data.keys()) - matched_tvks
    
    if unmatched:
        logger.info(f"JNCC TVKs not matched to any species: {len(unmatched)}")
        for tvk in sorted(unmatched):
            logger.debug(f"Unmatched JNCC TVK: {tvk}")
    else:
        logger.info("All JNCC TVKs matched to species")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate UKSI valid species with JNCC conservation designations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python uksi_jncc_annotation.py --input species.tsv --jncc jncc.tsv
    
    # With output directory
    python uksi_jncc_annotation.py --input species.tsv --jncc jncc.tsv --output-dir /results
    
    # Output to stdout (for piping in snakemake)
    python uksi_jncc_annotation.py --input species.tsv --jncc jncc.tsv --stdout
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
        '--stdout',
        action='store_true',
        help='Write annotated TSV to stdout instead of file (for pipeline integration)'
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
    logger.info(f"Starting JNCC annotation")
    logger.info(f"Input file: {args.input}")
    logger.info(f"JNCC file: {args.jncc}")
    logger.info(f"Output: {'stdout' if args.stdout else output_file}")
    
    try:
        # Load JNCC data
        jncc_data = load_jncc_data(args.jncc, logger)
        
        # Process species and annotate
        matched_tvks = process_species(
            args.input,
            jncc_data,
            output_file if not args.stdout else None,
            args.stdout,
            logger
        )
        
        # Log unmatched JNCC entries
        log_unmatched_jncc(jncc_data, matched_tvks, logger)
        
        logger.info("JNCC annotation complete")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
