#!/usr/bin/env python3
"""
Pantheon to UKSI Name Matcher

Matches taxon names from the Pantheon database against the UK Species Inventory (UKSI)
to retrieve taxonomic keys for downstream processing.

Author: Ben Price / Claude
Date: 2025-01-25
"""

import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime
import sys


def setup_logging(log_path: Path) -> list:
    """Initialize log buffer."""
    return []


def log_message(log_buffer: list, message: str, also_print: bool = True):
    """Add message to log buffer and optionally print to console."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"[{timestamp}] {message}"
    log_buffer.append(log_entry)
    if also_print:
        print(message)


def write_log(log_buffer: list, log_path: Path):
    """Write log buffer to file."""
    with open(log_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(log_buffer))


def load_uksi(uksi_path: Path, log_buffer: list) -> pd.DataFrame:
    """
    Load UKSI names file and prepare for matching.
    Prioritises recommended names (NAME_STATUS = 'R') over synonyms.
    """
    log_message(log_buffer, f"Loading UKSI file: {uksi_path}")
    
    df = pd.read_csv(uksi_path, sep='\t', dtype=str, keep_default_na=False)
    log_message(log_buffer, f"  Total UKSI records loaded: {len(df):,}")
    
    # Create lowercase name column for matching
    df['TAXON_NAME_LOWER'] = df['TAXON_NAME'].str.lower().str.strip()
    
    # Sort to prioritise recommended names (R before S, U, I, etc.)
    # R = Recommended, S = Synonym, U = Unverified, I = Invalid
    status_priority = {'R': 0, 'S': 1, 'U': 2, 'I': 3}
    df['STATUS_PRIORITY'] = df['NAME_STATUS'].map(status_priority).fillna(4)
    df = df.sort_values(['TAXON_NAME_LOWER', 'STATUS_PRIORITY'])
    
    # Keep first occurrence per lowercase name (will be the recommended one if available)
    df_dedup = df.drop_duplicates(subset=['TAXON_NAME_LOWER'], keep='first')
    
    log_message(log_buffer, f"  Unique UKSI names for matching: {len(df_dedup):,}")
    
    # Count by status
    status_counts = df_dedup['NAME_STATUS'].value_counts()
    for status, count in status_counts.items():
        log_message(log_buffer, f"    NAME_STATUS '{status}': {count:,}")
    
    return df_dedup


def load_pantheon(pantheon_path: Path, log_buffer: list) -> pd.DataFrame:
    """
    Load Pantheon file and deduplicate by Taxon name.
    Combines duplicate rows by joining values with commas.
    """
    log_message(log_buffer, f"Loading Pantheon file: {pantheon_path}")
    
    df = pd.read_csv(pantheon_path, sep='\t', dtype=str, keep_default_na=False)
    log_message(log_buffer, f"  Total Pantheon records loaded: {len(df):,}")
    
    # Check for Taxon column
    if 'Taxon' not in df.columns:
        raise ValueError("Pantheon file must contain a 'Taxon' column")
    
    # Check for original_taxon fallback column
    has_original_taxon = 'original_taxon' in df.columns
    if has_original_taxon:
        log_message(log_buffer, "  Found 'original_taxon' column for fallback matching")
    
    # Count duplicates before deduplication
    dup_counts = df['Taxon'].value_counts()
    n_duplicated = (dup_counts > 1).sum()
    log_message(log_buffer, f"  Taxa with duplicate rows: {n_duplicated:,}")
    
    # Deduplicate by combining rows with same Taxon
    # For each column, join unique non-empty values with comma
    def combine_values(series):
        unique_vals = series[series != ''].unique()
        return ', '.join(unique_vals)
    
    # Group by Taxon and aggregate
    df_dedup = df.groupby('Taxon', as_index=False).agg(
        lambda x: combine_values(x) if x.name != 'Taxon' else x.iloc[0]
    )
    
    log_message(log_buffer, f"  Unique Pantheon taxa after deduplication: {len(df_dedup):,}")
    
    # Create lowercase column for matching
    df_dedup['TAXON_LOWER'] = df_dedup['Taxon'].str.lower().str.strip()
    
    # Create lowercase fallback column if original_taxon exists
    if has_original_taxon:
        df_dedup['ORIGINAL_TAXON_LOWER'] = df_dedup['original_taxon'].str.lower().str.strip()
    
    return df_dedup


def match_names(pantheon_df: pd.DataFrame, uksi_df: pd.DataFrame, 
                log_buffer: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match Pantheon taxa against UKSI names.
    First tries 'Taxon' column, then falls back to 'original_taxon' if available.
    Returns (matched_df, unmatched_df).
    """
    log_message(log_buffer, "Matching Pantheon taxa against UKSI...")
    
    # Check if fallback column exists
    has_fallback = 'ORIGINAL_TAXON_LOWER' in pantheon_df.columns
    if has_fallback:
        log_message(log_buffer, "  Fallback column 'original_taxon' available")
    
    # Create lookup dict from UKSI: lowercase name -> (TVK, RTVK, NAME_STATUS)
    uksi_lookup = {}
    for _, row in uksi_df.iterrows():
        name_lower = row['TAXON_NAME_LOWER']
        uksi_lookup[name_lower] = (
            row['TAXON_VERSION_KEY'],
            row['RECOMMENDED_TAXON_VERSION_KEY'],
            row['NAME_STATUS']
        )
    
    # Match each Pantheon taxon
    matched_rows = []
    unmatched_rows = []
    matched_on_taxon = 0
    matched_on_fallback = 0
    
    for _, row in pantheon_df.iterrows():
        taxon_lower = row['TAXON_LOWER']
        row_dict = row.to_dict()
        
        # Try primary match on Taxon column
        if taxon_lower in uksi_lookup:
            tvk, rtvk, status = uksi_lookup[taxon_lower]
            row_dict['TAXON_VERSION_KEY'] = tvk
            row_dict['RECOMMENDED_TAXON_VERSION_KEY'] = rtvk
            row_dict['NAME_STATUS'] = status
            row_dict['MATCHED_ON'] = 'Taxon'
            matched_rows.append(row_dict)
            matched_on_taxon += 1
        # Try fallback on original_taxon column
        elif has_fallback and row['ORIGINAL_TAXON_LOWER'] in uksi_lookup:
            fallback_lower = row['ORIGINAL_TAXON_LOWER']
            tvk, rtvk, status = uksi_lookup[fallback_lower]
            row_dict['TAXON_VERSION_KEY'] = tvk
            row_dict['RECOMMENDED_TAXON_VERSION_KEY'] = rtvk
            row_dict['NAME_STATUS'] = status
            row_dict['MATCHED_ON'] = 'original_taxon'
            matched_rows.append(row_dict)
            matched_on_fallback += 1
        else:
            unmatched_rows.append(row_dict)
    
    matched_df = pd.DataFrame(matched_rows)
    unmatched_df = pd.DataFrame(unmatched_rows)
    
    # Remove helper columns
    helper_cols = ['TAXON_LOWER', 'ORIGINAL_TAXON_LOWER']
    for col in helper_cols:
        if len(matched_df) > 0 and col in matched_df.columns:
            matched_df = matched_df.drop(columns=[col])
        if len(unmatched_df) > 0 and col in unmatched_df.columns:
            unmatched_df = unmatched_df.drop(columns=[col])
    
    # Log results
    total = len(pantheon_df)
    n_matched = len(matched_df)
    n_unmatched = len(unmatched_df)
    match_pct = (n_matched / total * 100) if total > 0 else 0
    
    log_message(log_buffer, f"  Matched: {n_matched:,} ({match_pct:.1f}%)")
    log_message(log_buffer, f"    - on Taxon: {matched_on_taxon:,}")
    if has_fallback:
        log_message(log_buffer, f"    - on original_taxon (fallback): {matched_on_fallback:,}")
    log_message(log_buffer, f"  Unmatched: {n_unmatched:,} ({100-match_pct:.1f}%)")
    
    # Summary of matched name statuses
    if len(matched_df) > 0:
        status_counts = matched_df['NAME_STATUS'].value_counts()
        log_message(log_buffer, "  Matched by NAME_STATUS:")
        for status, count in status_counts.items():
            status_label = {'R': 'Recommended', 'S': 'Synonym', 
                          'U': 'Unverified', 'I': 'Invalid'}.get(status, status)
            log_message(log_buffer, f"    {status} ({status_label}): {count:,}")
    
    return matched_df, unmatched_df


def main():
    parser = argparse.ArgumentParser(
        description='Match Pantheon taxon names against UKSI to retrieve taxonomic keys.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python pantheon_uksi_matcher.py --uksi uksi_names.tsv --pantheon pantheon_input.tsv --output-dir ./output

Output files:
  - {input_name}_matched.tsv: Pantheon records with UKSI keys appended
  - {input_name}_unmatched.tsv: Pantheon records not found in UKSI
  - {input_name}_log.txt: Processing log with statistics
        """
    )
    
    parser.add_argument('--uksi', required=True, type=Path,
                        help='Path to UKSI names TSV file')
    parser.add_argument('--pantheon', required=True, type=Path,
                        help='Path to Pantheon input TSV file')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for results')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.uksi.exists():
        print(f"Error: UKSI file not found: {args.uksi}")
        sys.exit(1)
    if not args.pantheon.exists():
        print(f"Error: Pantheon file not found: {args.pantheon}")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate output filenames based on input pantheon name
    base_name = args.pantheon.stem
    matched_path = args.output_dir / f"{base_name}_matched.tsv"
    unmatched_path = args.output_dir / f"{base_name}_unmatched.tsv"
    log_path = args.output_dir / f"{base_name}_log.txt"
    
    # Initialize logging
    log_buffer = setup_logging(log_path)
    
    log_message(log_buffer, "=" * 60)
    log_message(log_buffer, "Pantheon to UKSI Name Matcher")
    log_message(log_buffer, "=" * 60)
    log_message(log_buffer, f"UKSI input: {args.uksi}")
    log_message(log_buffer, f"Pantheon input: {args.pantheon}")
    log_message(log_buffer, f"Output directory: {args.output_dir}")
    log_message(log_buffer, "")
    
    try:
        # Load data
        uksi_df = load_uksi(args.uksi, log_buffer)
        log_message(log_buffer, "")
        
        pantheon_df = load_pantheon(args.pantheon, log_buffer)
        log_message(log_buffer, "")
        
        # Match names
        matched_df, unmatched_df = match_names(pantheon_df, uksi_df, log_buffer)
        log_message(log_buffer, "")
        
        # Write outputs
        log_message(log_buffer, "Writing output files...")
        
        if len(matched_df) > 0:
            matched_df.to_csv(matched_path, sep='\t', index=False)
            log_message(log_buffer, f"  Matched: {matched_path}")
        else:
            log_message(log_buffer, "  No matched records to write")
        
        if len(unmatched_df) > 0:
            unmatched_df.to_csv(unmatched_path, sep='\t', index=False)
            log_message(log_buffer, f"  Unmatched: {unmatched_path}")
        else:
            log_message(log_buffer, "  No unmatched records to write")
        
        log_message(log_buffer, "")
        log_message(log_buffer, "Processing complete.")
        
    except Exception as e:
        log_message(log_buffer, f"ERROR: {e}")
        write_log(log_buffer, log_path)
        raise
    
    # Write final log
    write_log(log_buffer, log_path)
    print(f"Log written to: {log_path}")


if __name__ == '__main__':
    main()
