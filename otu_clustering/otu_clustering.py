#!/usr/bin/env python3
"""
OTU Clustering Script

Clusters DNA sequences into Operational Taxonomic Units (OTUs) using VSEARCH.
Accepts TSV input with sequence data and outputs the same TSV with an additional
OTU_ID column appended.

Sequences can be in either orientation - the script handles both forward and
reverse complement sequences through VSEARCH's strand parameter.

Author: Ben Price / Claude
Date: 2025-01-22
"""

import argparse
import sys
import subprocess
import tempfile
import shutil
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def log_message(level: str, message: str, verbose: bool = True) -> None:
    """
    Print a log message to stderr.
    
    Args:
        level: Log level (INFO, WARNING, ERROR)
        message: Message to print
        verbose: Only print if True
    """
    if verbose or level == 'ERROR':
        from datetime import datetime
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"[{timestamp}] {level}: {message}", file=sys.stderr)


def find_vsearch_binary() -> Optional[str]:
    """
    Find the VSEARCH binary in various possible locations.
    
    Returns:
        Path to vsearch binary or None if not found
    """
    import os
    
    # List of possible locations
    possible_paths = [
        'vsearch',                      # In PATH
        '/usr/bin/vsearch',             # System install
        '/usr/local/bin/vsearch',       # Local install
    ]
    
    # Check conda environment
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        possible_paths.insert(0, os.path.join(conda_prefix, 'bin', 'vsearch'))
    
    # Try which command first (most reliable for finding executables)
    try:
        result = subprocess.run(['which', 'vsearch'], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except Exception:
        pass
    
    # Windows: try where command
    try:
        result = subprocess.run(['where', 'vsearch'], capture_output=True, text=True, shell=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip().split('\n')[0]
    except Exception:
        pass
    
    # Check each possible path
    for path in possible_paths:
        expanded = os.path.expandvars(path)
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
    
    return None


def validate_sequence(sequence: str, min_length: int = 100) -> Tuple[bool, str]:
    """
    Validate and clean a nucleotide sequence.
    
    Args:
        sequence: Raw sequence string
        min_length: Minimum acceptable length
        
    Returns:
        Tuple of (is_valid, cleaned_sequence)
    """
    if not sequence:
        return False, ''
    
    # Clean sequence - remove whitespace and convert to uppercase
    cleaned = re.sub(r'\s+', '', sequence).upper()
    
    # Trim leading and trailing gaps (alignment artifacts)
    cleaned = cleaned.strip('-')
    
    # Replace internal gaps with N (unknown nucleotide)
    cleaned = cleaned.replace('-', 'N')
    
    # Check for valid IUPAC nucleotide codes only
    # VSEARCH accepts: ACGTURYSWKMDBHVN
    if re.search(r'[^ACGTURYSWKMDBHVN]', cleaned):
        return False, cleaned
    
    # Check minimum length
    if len(cleaned) < min_length:
        return False, cleaned
    
    return True, cleaned


def process_input_file(input_path: Path, accession_col: str,
                       sequence_cols: List[str], min_length: int,
                       verbose: bool = False
                       ) -> Tuple[List[str], List[List[str]], Dict[str, str], int, int, int]:
    """
    Read and process input TSV file, extracting valid sequences for clustering.
    
    Args:
        input_path: Path to input TSV file
        accession_col: Name of accession column
        sequence_cols: List of possible sequence column names
        min_length: Minimum sequence length
        verbose: Print progress information
        
    Returns:
        Tuple of (headers, all_rows, valid_sequences dict, acc_idx, seq_idx, valid_count)
    """
    headers = []
    all_rows = []
    valid_sequences = {}
    valid_count = 0
    no_sequence_count = 0
    
    with open(input_path, 'r', encoding='utf-8') as f:
        # Read header
        header_line = f.readline().rstrip('\n\r')
        headers = header_line.split('\t')
        
        # Find accession column index
        try:
            acc_idx = headers.index(accession_col)
        except ValueError:
            raise ValueError(f"Accession column '{accession_col}' not found in header. "
                           f"Available columns: {', '.join(headers)}")
        
        # Find sequence column index
        seq_idx = None
        seq_col_found = None
        for col in sequence_cols:
            if col in headers:
                seq_idx = headers.index(col)
                seq_col_found = col
                break
        
        if seq_idx is None:
            raise ValueError(f"No sequence column found. Tried: {', '.join(sequence_cols)}. "
                           f"Available columns: {', '.join(headers)}")
        
        if verbose:
            log_message('INFO', f"Using accession column: '{accession_col}' (index {acc_idx})")
            log_message('INFO', f"Using sequence column: '{seq_col_found}' (index {seq_idx})")
        
        # Read data rows
        line_num = 1
        for line in f:
            line_num += 1
            line = line.rstrip('\n\r')
            
            # Store the raw row (even if empty - preserve all rows)
            fields = line.split('\t')
            all_rows.append(fields)
            
            # Skip empty lines for sequence extraction
            if not line.strip():
                no_sequence_count += 1
                continue
            
            # Get accession - handle rows with fewer columns
            if len(fields) <= acc_idx:
                no_sequence_count += 1
                continue
            
            accession = fields[acc_idx].strip()
            if not accession:
                no_sequence_count += 1
                continue
            
            # Get sequence - handle rows with fewer columns
            if len(fields) <= seq_idx:
                no_sequence_count += 1
                continue
            
            sequence = fields[seq_idx].strip()
            
            # Validate sequence
            is_valid, cleaned_seq = validate_sequence(sequence, min_length)
            
            if is_valid:
                valid_sequences[accession] = cleaned_seq
                valid_count += 1
            else:
                no_sequence_count += 1
            
            if verbose and line_num % 100000 == 0:
                log_message('INFO', f"Read {line_num} rows...")
    
    if verbose:
        log_message('INFO', f"Total rows: {len(all_rows)}")
        log_message('INFO', f"Valid sequences for clustering: {valid_count}")
        log_message('INFO', f"Rows without valid sequence: {no_sequence_count}")
    
    return headers, all_rows, valid_sequences, acc_idx, seq_idx, valid_count


def write_fasta(sequences: Dict[str, str], output_path: Path) -> None:
    """
    Write sequences to FASTA format.
    
    Args:
        sequences: Dictionary of accession -> sequence
        output_path: Path to output FASTA file
    """
    with open(output_path, 'w') as f:
        for accession, sequence in sequences.items():
            f.write(f">{accession}\n{sequence}\n")


def run_vsearch_clustering(fasta_path: Path, output_path: Path,
                          threshold: float, threads: int,
                          strand: str, verbose: bool = False) -> bool:
    """
    Run VSEARCH clustering.
    
    Args:
        fasta_path: Path to input FASTA file
        output_path: Path for VSEARCH UC output
        threshold: Similarity threshold (0.0-1.0)
        threads: Number of threads to use
        strand: Strand option ('plus' or 'both')
        verbose: Print command details
        
    Returns:
        True if successful, False otherwise
    """
    vsearch_binary = find_vsearch_binary()
    
    if not vsearch_binary:
        raise RuntimeError("VSEARCH binary not found. Please ensure VSEARCH is installed and in PATH.")
    
    cmd = [
        vsearch_binary,
        '--cluster_fast', str(fasta_path),
        '--id', str(threshold),
        '--uc', str(output_path),
        '--threads', str(threads),
        '--strand', strand,
        '--quiet'
    ]
    
    if verbose:
        log_message('INFO', f"VSEARCH binary: {vsearch_binary}")
        log_message('INFO', f"VSEARCH command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        log_message('ERROR', f"VSEARCH failed: {result.stderr}")
        return False
    
    return True


def parse_vsearch_uc(uc_path: Path) -> Dict[str, str]:
    """
    Parse VSEARCH UC format output to extract OTU assignments.
    
    Args:
        uc_path: Path to UC file
        
    Returns:
        Dictionary of accession -> OTU_ID
    """
    otu_assignments = {}
    otu_centroids = {}  # Maps centroid accession to OTU_ID
    otu_counter = 1
    
    with open(uc_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            
            # VSEARCH UC format:
            # Type H = hit (member of cluster)
            # Type C = centroid (cluster representative)
            # Type S = singleton (single member cluster)
            
            record_type = fields[0]
            query_id = fields[8]
            
            if record_type in ('C', 'S'):
                # This is a centroid/singleton - create new OTU
                otu_id = f"OTU_{otu_counter:06d}"
                otu_assignments[query_id] = otu_id
                otu_centroids[query_id] = otu_id
                otu_counter += 1
            elif record_type == 'H':
                # This is a hit - assign to existing OTU
                target_id = fields[9]
                if target_id in otu_centroids:
                    otu_assignments[query_id] = otu_centroids[target_id]
    
    return otu_assignments


def calculate_statistics(otu_assignments: Dict[str, str], total_rows: int) -> Dict[str, any]:
    """
    Calculate summary statistics for OTU clustering.
    
    Args:
        otu_assignments: Dictionary of accession -> OTU_ID
        total_rows: Total number of data rows in input
        
    Returns:
        Dictionary of statistics
    """
    # Count sequences per OTU
    otu_counts = {}
    for otu_id in otu_assignments.values():
        otu_counts[otu_id] = otu_counts.get(otu_id, 0) + 1
    
    total_otus = len(otu_counts)
    total_seqs = len(otu_assignments)
    
    # Count singletons
    singletons = sum(1 for count in otu_counts.values() if count == 1)
    
    # Find largest OTU
    largest_otu = max(otu_counts.items(), key=lambda x: x[1]) if otu_counts else ('', 0)
    
    # Calculate mean and median OTU size
    sizes = list(otu_counts.values())
    mean_size = sum(sizes) / len(sizes) if sizes else 0
    
    sorted_sizes = sorted(sizes)
    n = len(sorted_sizes)
    if n == 0:
        median_size = 0
    elif n % 2 == 1:
        median_size = sorted_sizes[n // 2]
    else:
        median_size = (sorted_sizes[n // 2 - 1] + sorted_sizes[n // 2]) / 2
    
    return {
        'total_rows': total_rows,
        'clustered_sequences': total_seqs,
        'rows_without_otu': total_rows - total_seqs,
        'total_otus': total_otus,
        'singletons': singletons,
        'largest_otu_id': largest_otu[0],
        'largest_otu_size': largest_otu[1],
        'mean_otu_size': mean_size,
        'median_otu_size': median_size
    }


def write_annotated_output(output_path: Path, headers: List[str], 
                           all_rows: List[List[str]], acc_idx: int,
                           otu_assignments: Dict[str, str]) -> None:
    """
    Write annotated TSV file with OTU_ID column appended.
    
    Args:
        output_path: Path to output file
        headers: Original header columns
        all_rows: All data rows from input
        acc_idx: Index of accession column
        otu_assignments: Dictionary of accession -> OTU_ID
    """
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header with OTU_ID column appended
        f.write('\t'.join(headers) + '\tOTU_ID\n')
        
        # Write data rows with OTU_ID appended
        for fields in all_rows:
            # Get accession if available
            accession = ''
            if len(fields) > acc_idx:
                accession = fields[acc_idx].strip()
            
            # Look up OTU assignment (empty string if no valid sequence)
            otu_id = otu_assignments.get(accession, '')
            
            # Write original row with OTU_ID appended
            f.write('\t'.join(fields) + '\t' + otu_id + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Cluster DNA sequences into OTUs using VSEARCH.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic clustering at 99% identity
  python otu_clustering.py -t 0.99 input.tsv output.tsv
  
  # Clustering at 97% identity with more threads
  python otu_clustering.py -t 0.97 --threads 16 sequences.tsv otus.tsv
  
  # Specify custom column names
  python otu_clustering.py -t 0.99 --accession-col recordid --sequence-col nuc input.tsv output.tsv
  
  # With verbose output
  python otu_clustering.py -v -t 0.99 input.tsv output.tsv

Input Format:
  Tab-separated file with at minimum an accession column and a sequence column.
  Default accession column: 'accession'
  Default sequence columns (tried in order): 'sequence', 'nucleotide_sequence', 'nuc'

Output Format:
  The input TSV file with an additional OTU_ID column appended.
  Rows without valid sequences will have an empty OTU_ID field.
'''
    )
    
    # Input/output arguments
    parser.add_argument('input', type=str,
                        help='Input TSV file with sequences (use "-" for stdin)')
    parser.add_argument('output', type=str,
                        help='Output TSV file with OTU_ID column appended (use "-" for stdout)')
    
    # Clustering parameters
    parser.add_argument('-t', '--threshold', type=float, required=True,
                        help='Similarity threshold for clustering (0.0-1.0, e.g., 0.99 for 99%%)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads for VSEARCH (default: 8)')
    parser.add_argument('--strand', choices=['plus', 'both'], default='both',
                        help='Strand orientation: "plus" (forward only) or "both" (default: both)')
    parser.add_argument('--min-length', type=int, default=100,
                        help='Minimum sequence length to include (default: 100)')
    
    # Column specification
    parser.add_argument('--accession-col', type=str, default='accession',
                        help='Name of accession column (default: accession)')
    parser.add_argument('--sequence-col', type=str, action='append',
                        help='Name of sequence column (can specify multiple, tried in order)')
    
    # Temp directory
    parser.add_argument('--temp-dir', type=str, default=None,
                        help='Temporary directory for intermediate files (default: system temp)')
    
    # Output options
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress and statistics to stderr')
    
    args = parser.parse_args()
    
    # Validate threshold
    if not 0.0 < args.threshold <= 1.0:
        print("Error: Threshold must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)
    
    # Set default sequence columns if not specified
    if args.sequence_col is None:
        sequence_cols = ['sequence', 'nucleotide_sequence', 'nuc']
    else:
        sequence_cols = args.sequence_col
    
    # Handle stdin
    if args.input == '-':
        # Write stdin to temp file
        temp_input = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        temp_input.write(sys.stdin.read())
        temp_input.close()
        input_path = Path(temp_input.name)
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)
    
    # Handle stdout
    if args.output == '-':
        output_path = Path('/dev/stdout')
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create temp directory
    temp_dir = Path(args.temp_dir) if args.temp_dir else Path(tempfile.mkdtemp(prefix='otu_'))
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        if args.verbose:
            log_message('INFO', "Starting OTU clustering")
            log_message('INFO', f"Input: {args.input}")
            log_message('INFO', f"Output: {args.output}")
            log_message('INFO', f"Threshold: {args.threshold}")
            log_message('INFO', f"Threads: {args.threads}")
            log_message('INFO', f"Strand: {args.strand}")
            log_message('INFO', f"Min length: {args.min_length}")
        
        # Read and process input file
        if args.verbose:
            log_message('INFO', "Reading input TSV file...")
        
        headers, all_rows, valid_sequences, acc_idx, seq_idx, valid_count = process_input_file(
            input_path, args.accession_col, sequence_cols, args.min_length, args.verbose
        )
        
        # Handle case with no valid sequences
        if valid_count == 0:
            if args.verbose:
                log_message('WARNING', "No valid sequences found for clustering")
            # Write output with empty OTU_ID column
            write_annotated_output(output_path, headers, all_rows, acc_idx, {})
            if args.verbose:
                log_message('INFO', "Output written with empty OTU_ID column")
            sys.exit(0)
        
        # Write FASTA for VSEARCH
        fasta_path = temp_dir / 'sequences.fasta'
        write_fasta(valid_sequences, fasta_path)
        
        if args.verbose:
            log_message('INFO', f"Wrote {len(valid_sequences)} sequences to FASTA for clustering")
        
        # Run VSEARCH clustering
        if args.verbose:
            log_message('INFO', "Running VSEARCH clustering...")
        
        uc_path = temp_dir / 'clusters.uc'
        success = run_vsearch_clustering(
            fasta_path, uc_path, args.threshold, args.threads, args.strand, args.verbose
        )
        
        if not success:
            print("Error: VSEARCH clustering failed", file=sys.stderr)
            sys.exit(1)
        
        if args.verbose:
            log_message('INFO', "VSEARCH clustering completed")
        
        # Parse results
        if args.verbose:
            log_message('INFO', "Parsing clustering results...")
        
        otu_assignments = parse_vsearch_uc(uc_path)
        
        # Write annotated output
        if args.verbose:
            log_message('INFO', f"Writing annotated output to {args.output}")
        
        write_annotated_output(output_path, headers, all_rows, acc_idx, otu_assignments)
        
        # Calculate and print statistics
        if args.verbose:
            stats = calculate_statistics(otu_assignments, len(all_rows))
            log_message('INFO', "Summary statistics:")
            log_message('INFO', f"  Total rows in input: {stats['total_rows']}")
            log_message('INFO', f"  Sequences clustered: {stats['clustered_sequences']}")
            log_message('INFO', f"  Rows without OTU assignment: {stats['rows_without_otu']}")
            log_message('INFO', f"  Total OTUs created: {stats['total_otus']}")
            log_message('INFO', f"  Singleton OTUs: {stats['singletons']}")
            if stats['largest_otu_id']:
                log_message('INFO', f"  Largest OTU: {stats['largest_otu_id']} ({stats['largest_otu_size']} sequences)")
            log_message('INFO', f"  Mean OTU size: {stats['mean_otu_size']:.2f}")
            log_message('INFO', f"  Median OTU size: {stats['median_otu_size']:.1f}")
            log_message('INFO', "OTU clustering completed successfully")
    
    finally:
        # Clean up temp files
        if args.input == '-':
            input_path.unlink(missing_ok=True)
        
        if not args.temp_dir:  # Only clean up if we created the temp dir
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == '__main__':
    main()
