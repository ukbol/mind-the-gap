#!/usr/bin/env python3
"""
MIDORI FASTA to TSV Converter

Parses MIDORI reference library FASTA files and extracts taxonomic information
into a structured TSV format suitable for downstream analysis.

Author: Ben Price / Claude
Date: 2025-01-22
"""

import argparse
import sys
import re
from pathlib import Path


def parse_taxonomy(header: str) -> dict:
    """
    Parse MIDORI FASTA header to extract taxonomic ranks and taxids.
    
    Args:
        header: FASTA header line (without leading '>')
        
    Returns:
        Dictionary with taxonomic ranks and metadata
    """
    # Split header into accession part and taxonomy part
    parts = header.split(' ', 1)
    
    # Extract accession (e.g., "KF807066" from "KF807066.1.<1.>569")
    accession_full = parts[0]
    accession = accession_full.split('.')[0]
    
    # Initialize result with empty values
    result = {
        'accession': accession,
        'kingdom': '',
        'phylum': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': '',
        'species': '',
        'taxid': ''
    }
    
    if len(parts) < 2:
        return result
    
    taxonomy_str = parts[1]
    
    # Split by semicolon to get individual taxonomic entries
    taxa = taxonomy_str.split(';')
    
    # Track the lowest taxid found
    lowest_taxid = ''
    
    # Target ranks we want to extract
    target_ranks = {'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}
    
    for taxon in taxa:
        taxon = taxon.strip()
        if not taxon:
            continue
            
        # Pattern: rank_name_taxid (e.g., "kingdom_Metazoa_33208" or "species_Silverstoneia nubicola_384849")
        # The taxid is always the last underscore-separated element
        # The rank is the first underscore-separated element
        # The name is everything in between
        
        # Find the last underscore to get taxid
        last_underscore = taxon.rfind('_')
        if last_underscore == -1:
            continue
            
        taxid_part = taxon[last_underscore + 1:]
        
        # Check if it's a valid taxid (numeric)
        if not taxid_part.isdigit():
            continue
            
        # Get the rank_name part
        rank_name_part = taxon[:last_underscore]
        
        # Find the first underscore to get rank
        first_underscore = rank_name_part.find('_')
        if first_underscore == -1:
            continue
            
        rank = rank_name_part[:first_underscore].lower()
        name = rank_name_part[first_underscore + 1:]
        
        # Store if it's a target rank
        if rank in target_ranks:
            result[rank] = name
            
        # Update lowest taxid (last valid one in the chain)
        lowest_taxid = taxid_part
    
    result['taxid'] = lowest_taxid
    
    return result


def process_fasta(input_file: Path, output_file: Path, verbose: bool = False) -> int:
    """
    Process MIDORI FASTA file and write TSV output.
    
    Args:
        input_file: Path to input FASTA file
        output_file: Path to output TSV file
        verbose: Print progress information
        
    Returns:
        Number of sequences processed
    """
    columns = ['accession', 'kingdom', 'phylum', 'class', 'order', 
               'family', 'genus', 'species', 'taxid', 'sequence']
    
    count = 0
    current_header = None
    current_sequence = []
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Write header
        fout.write('\t'.join(columns) + '\n')
        
        for line in fin:
            line = line.rstrip('\n\r')
            
            if line.startswith('>'):
                # Process previous record if exists
                if current_header is not None:
                    record = parse_taxonomy(current_header)
                    record['sequence'] = ''.join(current_sequence)
                    fout.write('\t'.join(record[col] for col in columns) + '\n')
                    count += 1
                    
                    if verbose and count % 10000 == 0:
                        print(f"Processed {count} sequences...", file=sys.stderr)
                
                # Start new record
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                # Accumulate sequence
                current_sequence.append(line)
        
        # Process final record
        if current_header is not None:
            record = parse_taxonomy(current_header)
            record['sequence'] = ''.join(current_sequence)
            fout.write('\t'.join(record[col] for col in columns) + '\n')
            count += 1
    
    return count


def main():
    parser = argparse.ArgumentParser(
        description='Convert MIDORI FASTA files to TSV format with structured taxonomy.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage
  python process_midori.py input.fasta output.tsv
  
  # With verbose output
  python process_midori.py -v MIDORI2_LONGEST_NUC_GB264_lrRNA_QIIME.fasta lrRNA_processed.tsv
  
  # Using stdin/stdout for pipeline integration
  zcat input.fasta.gz | python process_midori.py - - > output.tsv
'''
    )
    
    parser.add_argument('input', type=str,
                        help='Input FASTA file (use "-" for stdin)')
    parser.add_argument('output', type=str,
                        help='Output TSV file (use "-" for stdout)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress information to stderr')
    
    args = parser.parse_args()
    
    # Handle stdin/stdout
    if args.input == '-':
        input_path = Path('/dev/stdin')
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)
    
    if args.output == '-':
        output_path = Path('/dev/stdout')
    else:
        output_path = Path(args.output)
        # Create parent directories if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if args.verbose:
        print(f"Processing: {args.input}", file=sys.stderr)
    
    count = process_fasta(input_path, output_path, args.verbose)
    
    if args.verbose:
        print(f"Complete. Processed {count} sequences.", file=sys.stderr)


if __name__ == '__main__':
    main()
