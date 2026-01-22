#!/usr/bin/env python3
"""
UNITE FASTA to TSV Converter

Parses UNITE reference library FASTA files and extracts taxonomic information
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
    Parse UNITE FASTA header to extract taxonomic ranks and cluster info.
    
    Header format:
    >Name|ACCESSION|CLUSTER|type|k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    
    Args:
        header: FASTA header line (without leading '>')
        
    Returns:
        Dictionary with taxonomic ranks and metadata
    """
    # Initialize result with empty values
    result = {
        'accession': '',
        'cluster': '',
        'kingdom': '',
        'phylum': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': '',
        'species': ''
    }
    
    # Split by pipe to get sections
    parts = header.split('|')
    
    if len(parts) < 5:
        return result
    
    # Accession is after the first pipe (index 1)
    result['accession'] = parts[1]
    
    # Cluster is the third section (index 2)
    result['cluster'] = parts[2]
    
    # Taxonomy string is the last section (index 4)
    taxonomy_str = parts[4]
    
    # Mapping of UNITE prefixes to our column names
    prefix_map = {
        'k__': 'kingdom',
        'p__': 'phylum',
        'c__': 'class',
        'o__': 'order',
        'f__': 'family',
        'g__': 'genus',
        's__': 'species'
    }
    
    # Split taxonomy by semicolon
    taxa = taxonomy_str.split(';')
    
    for taxon in taxa:
        taxon = taxon.strip()
        if not taxon:
            continue
            
        # Check each prefix
        for prefix, column in prefix_map.items():
            if taxon.startswith(prefix):
                # Extract the name after the prefix
                name = taxon[len(prefix):]
                result[column] = name
                break
    
    return result


def process_fasta(input_file: Path, output_file: Path, verbose: bool = False) -> int:
    """
    Process UNITE FASTA file and write TSV output.
    
    Args:
        input_file: Path to input FASTA file
        output_file: Path to output TSV file
        verbose: Print progress information
        
    Returns:
        Number of sequences processed
    """
    columns = ['accession', 'cluster', 'kingdom', 'phylum', 'class', 'order', 
               'family', 'genus', 'species', 'sequence']
    
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
        description='Convert UNITE FASTA files to TSV format with structured taxonomy.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage
  python process_unite.py input.fasta output.tsv
  
  # With verbose output
  python process_unite.py -v sh_general_release_all_10.0.fasta unite_processed.tsv
  
  # Using stdin/stdout for pipeline integration
  zcat input.fasta.gz | python process_unite.py - - > output.tsv
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
