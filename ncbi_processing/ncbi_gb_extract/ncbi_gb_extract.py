#!/usr/bin/env python3
"""
NCBI GenBank to TSV Converter

Parses NCBI GenBank flat files and extracts structured data for a specified
target gene into TSV format suitable for downstream analysis.

Author: Ben Price / Claude
Date: 2025-01-22
"""

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation


def extract_locus_info(record) -> Dict[str, str]:
    """
    Extract LOCUS line components from BioPython record annotations.
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        Dictionary with locus components
    """
    return {
        'locus_name': record.name,
        'locus_length': str(len(record.seq)),
        'locus_mol_type': record.annotations.get('molecule_type', ''),
        'locus_topology': record.annotations.get('topology', ''),
        'locus_division': record.annotations.get('data_file_division', ''),
        'locus_date': record.annotations.get('date', '')
    }


def extract_sequence_region(full_sequence: str, location) -> str:
    """
    Extract nucleotide sequence based on feature location.
    
    Handles simple locations, complement, join, and partial indicators.
    
    Args:
        full_sequence: Complete sequence from ORIGIN
        location: BioPython location object
        
    Returns:
        Extracted nucleotide sequence (reverse complemented if needed)
    """
    try:
        return str(location.extract(full_sequence))
    except Exception:
        return ''


def format_location(location) -> str:
    """
    Format a BioPython location object as a string.
    
    Args:
        location: BioPython location object
        
    Returns:
        String representation of the location
    """
    try:
        # Handle compound locations (join)
        if isinstance(location, CompoundLocation):
            parts = []
            for part in location.parts:
                start = int(part.start) + 1  # Convert to 1-based
                end = int(part.end)
                if part.strand == -1:
                    parts.append(f"complement({start}..{end})")
                else:
                    parts.append(f"{start}..{end}")
            if location.strand == -1:
                return f"complement(join({','.join(parts)}))"
            return f"join({','.join(parts)})"
        else:
            start = int(location.start) + 1  # Convert to 1-based
            end = int(location.end)
            if location.strand == -1:
                return f"complement({start}..{end})"
            return f"{start}..{end}"
    except Exception:
        return str(location)


def get_organism_lineage(record) -> str:
    """
    Extract full organism lineage from record annotations.
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        Semicolon-separated lineage string
    """
    taxonomy = record.annotations.get('taxonomy', [])
    if taxonomy:
        return '; '.join(taxonomy)
    return ''


def get_reference_info(record) -> Dict[str, str]:
    """
    Extract reference information (first reference only).
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        Dictionary with authors, title, journal
    """
    result = {
        'ref_authors': '',
        'ref_title': '',
        'ref_journal': ''
    }
    
    references = record.annotations.get('references', [])
    if references:
        ref = references[0]
        if hasattr(ref, 'authors') and ref.authors:
            result['ref_authors'] = ref.authors
        if hasattr(ref, 'title') and ref.title:
            result['ref_title'] = ref.title
        if hasattr(ref, 'journal') and ref.journal:
            result['ref_journal'] = ref.journal
    
    return result


def find_target_cds_features(record, target_gene: str) -> List[SeqFeature]:
    """
    Find all CDS features matching the target gene name (case-insensitive).
    
    Args:
        record: BioPython SeqRecord
        target_gene: Gene name to search for
        
    Returns:
        List of matching CDS features
    """
    target_lower = target_gene.lower()
    matching_cds = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            gene_name = feature.qualifiers.get('gene', [''])[0]
            if gene_name.lower() == target_lower:
                matching_cds.append(feature)
    
    return matching_cds


def find_source_feature(record) -> Optional[SeqFeature]:
    """
    Find the source feature in a record.
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        Source feature or None
    """
    for feature in record.features:
        if feature.type == 'source':
            return feature
    return None


def extract_qualifiers(feature: SeqFeature, exclude_keys: Set[str] = None) -> Dict[str, str]:
    """
    Extract all qualifiers from a feature as flat key-value pairs.
    
    Args:
        feature: BioPython SeqFeature
        exclude_keys: Keys to exclude from extraction
        
    Returns:
        Dictionary of qualifier key-value pairs
    """
    if exclude_keys is None:
        exclude_keys = set()
    
    result = {}
    for key, values in feature.qualifiers.items():
        if key in exclude_keys:
            continue
        # Join multiple values with semicolon
        if isinstance(values, list):
            result[key] = '; '.join(str(v) for v in values)
        else:
            result[key] = str(values)
    
    return result


def process_record(record, target_gene: str, all_columns: Set[str], 
                   verbose: bool = False) -> Tuple[Optional[Dict[str, str]], str]:
    """
    Process a single GenBank record and extract data for target gene.
    
    Args:
        record: BioPython SeqRecord
        target_gene: Gene name to extract
        all_columns: Set of all columns seen so far (will be updated)
        verbose: Print warnings
        
    Returns:
        Tuple of (data dict or None if skipped, status message)
    """
    accession = record.id
    
    # Find matching CDS features
    cds_features = find_target_cds_features(record, target_gene)
    
    if not cds_features:
        return None, f"SKIPPED: {accession} - no CDS feature for gene '{target_gene}'"
    
    # Warn if multiple CDS features found
    if len(cds_features) > 1 and verbose:
        msg = f"WARNING: {accession} - found {len(cds_features)} CDS features for gene '{target_gene}', consolidating"
        print(msg, file=sys.stderr)
    
    # Initialize result with header information
    result = {}
    
    # LOCUS information
    locus_data = extract_locus_info(record)
    result.update(locus_data)
    
    # Basic header fields
    result['definition'] = record.description
    result['accession'] = record.annotations.get('accessions', [accession])[0] if record.annotations.get('accessions') else accession
    result['version'] = record.annotations.get('sequence_version', '')
    if result['version']:
        result['version'] = f"{result['accession']}.{result['version']}"
    else:
        result['version'] = record.id
    result['keywords'] = '; '.join(record.annotations.get('keywords', []))
    result['organism'] = record.annotations.get('organism', '')
    result['taxonomy'] = get_organism_lineage(record)
    
    # Reference information
    ref_info = get_reference_info(record)
    result.update(ref_info)
    
    # Comment/structured comment
    comment = record.annotations.get('comment', '')
    result['comment'] = comment
    
    # Source feature qualifiers
    source_feature = find_source_feature(record)
    if source_feature:
        source_quals = extract_qualifiers(source_feature)
        for key, value in source_quals.items():
            col_name = f"source_{key}"
            result[col_name] = value
            all_columns.add(col_name)
    
    # Full sequence for extraction
    full_sequence = record.seq
    
    # Process CDS features
    for i, cds in enumerate(cds_features):
        suffix = '' if i == 0 else f'_{i+1}'
        
        # CDS location
        col_name = f'cds_location{suffix}'
        result[col_name] = format_location(cds.location)
        all_columns.add(col_name)
        
        # CDS qualifiers (exclude translation as we're extracting nucleotides)
        cds_quals = extract_qualifiers(cds, exclude_keys={'translation'})
        for key, value in cds_quals.items():
            col_name = f'cds_{key}{suffix}'
            result[col_name] = value
            all_columns.add(col_name)
        
        # Extract nucleotide sequence
        col_name = f'nucleotide_sequence{suffix}'
        result[col_name] = extract_sequence_region(full_sequence, cds.location)
        all_columns.add(col_name)
    
    # Update all_columns with any new keys
    all_columns.update(result.keys())
    
    return result, f"PROCESSED: {accession}"


def get_ordered_columns(all_columns: Set[str]) -> List[str]:
    """
    Return columns in a logical order.
    
    Args:
        all_columns: Set of all column names
        
    Returns:
        Ordered list of column names
    """
    # Define preferred order for known columns
    priority_order = [
        'locus_name', 'locus_length', 'locus_mol_type', 'locus_topology', 
        'locus_division', 'locus_date',
        'definition', 'accession', 'version', 'keywords',
        'organism', 'taxonomy',
        'ref_authors', 'ref_title', 'ref_journal',
        'comment'
    ]
    
    # Source columns
    source_cols = sorted([c for c in all_columns if c.startswith('source_')])
    
    # CDS columns (primary)
    cds_primary = sorted([c for c in all_columns if c.startswith('cds_') and not re.search(r'_\d+$', c)])
    
    # Nucleotide sequence (primary)
    nuc_primary = ['nucleotide_sequence'] if 'nucleotide_sequence' in all_columns else []
    
    # Additional CDS columns (from multiple features)
    cds_additional = sorted([c for c in all_columns if c.startswith('cds_') and re.search(r'_\d+$', c)])
    nuc_additional = sorted([c for c in all_columns if c.startswith('nucleotide_sequence_')])
    
    # Build ordered list
    ordered = []
    for col in priority_order:
        if col in all_columns:
            ordered.append(col)
    
    ordered.extend(source_cols)
    ordered.extend(cds_primary)
    ordered.extend(nuc_primary)
    ordered.extend(cds_additional)
    ordered.extend(nuc_additional)
    
    # Add any remaining columns not yet included
    remaining = sorted(all_columns - set(ordered))
    ordered.extend(remaining)
    
    return ordered


def process_genbank_file(input_path: Path, target_gene: str, 
                         all_columns: Set[str], all_records: List[Dict],
                         verbose: bool = False) -> Tuple[int, int]:
    """
    Process a single GenBank file.
    
    Args:
        input_path: Path to GenBank file
        target_gene: Gene name to extract
        all_columns: Set of all columns (updated in place)
        all_records: List of all records (updated in place)
        verbose: Print progress information
        
    Returns:
        Tuple of (processed count, skipped count)
    """
    processed = 0
    skipped = 0
    
    try:
        for record in SeqIO.parse(input_path, "genbank"):
            data, status = process_record(record, target_gene, all_columns, verbose)
            
            if data is None:
                skipped += 1
                if verbose:
                    print(status, file=sys.stderr)
            else:
                all_records.append(data)
                processed += 1
                
                if verbose and processed % 1000 == 0:
                    print(f"Processed {processed} records from {input_path.name}...", 
                          file=sys.stderr)
    except Exception as e:
        print(f"ERROR processing {input_path}: {e}", file=sys.stderr)
    
    return processed, skipped


def write_output(output_path: Path, all_records: List[Dict], 
                 all_columns: Set[str]) -> None:
    """
    Write all records to TSV output file.
    
    Args:
        output_path: Path to output TSV file
        all_records: List of all record dictionaries
        all_columns: Set of all column names
    """
    columns = get_ordered_columns(all_columns)
    
    with open(output_path, 'w', encoding='utf-8') as fout:
        # Write header
        fout.write('\t'.join(columns) + '\n')
        
        # Write records
        for record in all_records:
            row = []
            for col in columns:
                value = record.get(col, '')
                # Escape tabs and newlines in values
                if value:
                    value = str(value).replace('\t', ' ').replace('\n', ' ').replace('\r', '')
                row.append(value)
            fout.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Convert NCBI GenBank files to TSV format, extracting data for a specific gene.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Extract rbcL gene data from a single file
  python ncbi_gb_extract.py -g rbcL input.gb output.tsv
  
  # Process all GenBank files in a folder
  python ncbi_gb_extract.py -g COI -i /path/to/gb_files/ -o results.tsv
  
  # With verbose output
  python ncbi_gb_extract.py -v -g matK sequences.gb matK_data.tsv
  
  # Case-insensitive gene matching (these are equivalent)
  python ncbi_gb_extract.py -g rbcL input.gb output.tsv
  python ncbi_gb_extract.py -g RBCL input.gb output.tsv
  python ncbi_gb_extract.py -g RbcL input.gb output.tsv

Output Format:
  Tab-separated file with columns dynamically determined from input data.
  Core columns include locus info, definition, accession, version, organism,
  taxonomy, source qualifiers, CDS qualifiers, and nucleotide sequence.
  
  Records without the target gene CDS are skipped and logged to stderr.
  Multiple CDS features for the same gene are consolidated as additional
  columns (e.g., cds_protein_id, cds_protein_id_2).
'''
    )
    
    # Input options (mutually exclusive: single file or folder)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('input', type=str, nargs='?', default=None,
                             help='Input GenBank file (use "-" for stdin)')
    input_group.add_argument('-i', '--input-dir', type=str,
                             help='Input directory containing GenBank files (.gb, .gbk, .genbank)')
    
    parser.add_argument('output', type=str, nargs='?', default=None,
                        help='Output TSV file (use "-" for stdout)')
    parser.add_argument('-o', '--output-file', type=str,
                        help='Output TSV file (alternative to positional argument)')
    
    parser.add_argument('-g', '--gene', type=str, required=True,
                        help='Target gene name to extract (case-insensitive)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress and skipped records to stderr')
    
    args = parser.parse_args()
    
    # Resolve output path
    output_path_str = args.output or args.output_file
    if not output_path_str:
        print("Error: Output file must be specified", file=sys.stderr)
        sys.exit(1)
    
    # Collect input files
    input_files = []
    
    if args.input_dir:
        input_dir = Path(args.input_dir)
        if not input_dir.exists():
            print(f"Error: Input directory not found: {input_dir}", file=sys.stderr)
            sys.exit(1)
        
        # Find all GenBank files
        for ext in ['*.gb', '*.gbk', '*.genbank']:
            input_files.extend(input_dir.glob(ext))
        
        if not input_files:
            print(f"Error: No GenBank files found in {input_dir}", file=sys.stderr)
            sys.exit(1)
        
        input_files = sorted(input_files)
        
    elif args.input:
        if args.input == '-':
            # Handle stdin - write to temp file first (BioPython needs seekable file)
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as tmp:
                tmp.write(sys.stdin.read())
                input_files = [Path(tmp.name)]
        else:
            input_path = Path(args.input)
            if not input_path.exists():
                print(f"Error: Input file not found: {input_path}", file=sys.stderr)
                sys.exit(1)
            input_files = [input_path]
    
    # Process output path
    if output_path_str == '-':
        output_path = Path('/dev/stdout')
    else:
        output_path = Path(output_path_str)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if args.verbose:
        print(f"Target gene: {args.gene}", file=sys.stderr)
        print(f"Input files: {len(input_files)}", file=sys.stderr)
    
    # Process all files
    all_columns: Set[str] = set()
    all_records: List[Dict] = []
    total_processed = 0
    total_skipped = 0
    
    for input_file in input_files:
        if args.verbose:
            print(f"Processing: {input_file}", file=sys.stderr)
        
        processed, skipped = process_genbank_file(
            input_file, args.gene, all_columns, all_records, args.verbose
        )
        total_processed += processed
        total_skipped += skipped
    
    # Write output
    if args.verbose:
        print(f"Writing output to: {output_path_str}", file=sys.stderr)
    
    write_output(output_path, all_records, all_columns)
    
    if args.verbose:
        print(f"Complete. Processed: {total_processed}, Skipped: {total_skipped}", 
              file=sys.stderr)
    
    # Clean up temp file if used
    if args.input == '-':
        import os
        os.unlink(input_files[0])


if __name__ == '__main__':
    main()
