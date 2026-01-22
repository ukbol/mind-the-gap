#!/usr/bin/env python3
"""
BOLD Gene Extractor

Filters BOLD TSV data files to extract rows matching specified marker genes.
Supports case-insensitive matching and multiple gene names.

Author: Ben Price / Claude
Date: 2025-01-22
"""

import argparse
import sys
from pathlib import Path
from typing import List, Set, TextIO


def normalize_gene_name(gene: str) -> str:
    """
    Normalize gene name for case-insensitive comparison.
    
    Args:
        gene: Gene name to normalize
        
    Returns:
        Lowercase gene name
    """
    return gene.lower().strip()


def parse_gene_list(genes_arg: List[str]) -> Set[str]:
    """
    Parse gene arguments into a set of normalized gene names.
    
    Handles both comma-separated values and multiple arguments.
    
    Args:
        genes_arg: List of gene arguments (may contain comma-separated values)
        
    Returns:
        Set of normalized gene names
    """
    genes = set()
    for item in genes_arg:
        # Handle comma-separated values
        for gene in item.split(','):
            normalized = normalize_gene_name(gene)
            if normalized:
                genes.add(normalized)
    return genes


def find_marker_column(header: str, delimiter: str = '\t') -> int:
    """
    Find the index of the marker_code column in the header.
    
    Args:
        header: Header line from TSV file
        delimiter: Field delimiter
        
    Returns:
        Column index of marker_code
        
    Raises:
        ValueError: If marker_code column not found
    """
    columns = header.rstrip('\n\r').split(delimiter)
    try:
        return columns.index('marker_code')
    except ValueError:
        raise ValueError(
            f"Column 'marker_code' not found in header. "
            f"Available columns: {', '.join(columns[:10])}..."
        )


def process_file(
    input_file: TextIO,
    output_file: TextIO,
    target_genes: Set[str],
    delimiter: str = '\t',
    verbose: bool = False
) -> tuple:
    """
    Process a BOLD TSV file and extract rows matching target genes.
    
    Args:
        input_file: Input file handle
        output_file: Output file handle
        target_genes: Set of normalized target gene names
        delimiter: Field delimiter
        verbose: Print progress information
        
    Returns:
        Tuple of (processed_count, matched_count, skipped_count)
    """
    # Read and write header
    header = input_file.readline()
    if not header:
        raise ValueError("Input file is empty")
    
    marker_col = find_marker_column(header, delimiter)
    output_file.write(header)
    
    if verbose:
        print(f"Found marker_code at column index {marker_col}", file=sys.stderr)
    
    processed = 0
    matched = 0
    skipped = 0
    
    for line in input_file:
        processed += 1
        
        # Skip empty lines
        if not line.strip():
            skipped += 1
            continue
        
        fields = line.rstrip('\n\r').split(delimiter)
        
        # Check if we have enough columns
        if len(fields) <= marker_col:
            skipped += 1
            if verbose:
                print(f"WARNING: Line {processed + 1} has insufficient columns, skipping", 
                      file=sys.stderr)
            continue
        
        # Get marker code and check for match
        marker_code = normalize_gene_name(fields[marker_col])
        
        if marker_code in target_genes:
            output_file.write(line)
            matched += 1
        else:
            skipped += 1
        
        # Progress reporting
        if verbose and processed % 100000 == 0:
            print(f"Processed {processed:,} rows, matched {matched:,}...", 
                  file=sys.stderr)
    
    return processed, matched, skipped


def main():
    parser = argparse.ArgumentParser(
        description='Extract rows from BOLD TSV files matching specified marker genes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Extract rbcL sequences
  python bold_gene_extract.py -g rbcL input.tsv output.tsv
  
  # Extract multiple related genes (comma-separated)
  python bold_gene_extract.py -g rbcL,rbcLa input.tsv output.tsv
  
  # Extract multiple genes (multiple -g flags)
  python bold_gene_extract.py -g rbcL -g matK -g ITS2 input.tsv output.tsv
  
  # Case-insensitive matching (these are equivalent)
  python bold_gene_extract.py -g rbcL input.tsv output.tsv
  python bold_gene_extract.py -g RBCL input.tsv output.tsv
  python bold_gene_extract.py -g RbcL input.tsv output.tsv
  
  # Process from stdin/stdout for pipeline integration
  cat input.tsv | python bold_gene_extract.py -g COI-5P - - > output.tsv
  
  # With verbose output
  python bold_gene_extract.py -v -g rbcL,rbcLa large_dataset.tsv rbcL_only.tsv

Output Format:
  Tab-separated file with the same columns as input.
  Header row is always copied to output.
  Only rows where marker_code matches one of the target genes are included.

Notes:
  - Gene matching is case-insensitive
  - Multiple genes can be specified with comma separation or multiple -g flags
  - Both BOLD web downloads and data package formats are supported
'''
    )
    
    parser.add_argument('input', type=str,
                        help='Input BOLD TSV file (use "-" for stdin)')
    parser.add_argument('output', type=str,
                        help='Output TSV file (use "-" for stdout)')
    parser.add_argument('-g', '--gene', type=str, action='append', required=True,
                        dest='genes', metavar='GENE',
                        help='Target gene name(s) to extract. Can be specified multiple times '
                             'or as comma-separated values (e.g., -g rbcL,rbcLa or -g rbcL -g rbcLa). '
                             'Case-insensitive matching.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress and statistics to stderr')
    parser.add_argument('-d', '--delimiter', type=str, default='\t',
                        help='Field delimiter (default: tab)')
    
    args = parser.parse_args()
    
    # Parse target genes
    target_genes = parse_gene_list(args.genes)
    
    if not target_genes:
        print("Error: No valid gene names provided", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print(f"Target genes: {', '.join(sorted(target_genes))}", file=sys.stderr)
    
    # Set up input
    if args.input == '-':
        input_file = sys.stdin
        input_name = 'stdin'
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)
        input_file = open(input_path, 'r', encoding='utf-8')
        input_name = str(input_path)
    
    # Set up output
    if args.output == '-':
        output_file = sys.stdout
        output_name = 'stdout'
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_file = open(output_path, 'w', encoding='utf-8')
        output_name = str(output_path)
    
    try:
        if args.verbose:
            print(f"Processing: {input_name}", file=sys.stderr)
            print(f"Output: {output_name}", file=sys.stderr)
        
        processed, matched, skipped = process_file(
            input_file, output_file, target_genes,
            delimiter=args.delimiter, verbose=args.verbose
        )
        
        if args.verbose:
            print(f"Complete. Processed: {processed:,}, Matched: {matched:,}, "
                  f"Skipped: {skipped:,}", file=sys.stderr)
    
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    finally:
        # Close files if not stdin/stdout
        if args.input != '-':
            input_file.close()
        if args.output != '-':
            output_file.close()


if __name__ == '__main__':
    main()
