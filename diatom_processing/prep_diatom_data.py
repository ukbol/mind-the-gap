#!/usr/bin/env python3
"""
Diatom Dataset Preparation Script

Reformats the diatom TSV dataset for compatibility with the mind-the-gap
pipeline (otu_clustering.py, bags_assessment.py, gap_analysis.py).

Column transformations applied:
  - 'Sequence ID'   -> 'accession'   (otu_clustering.py uses exact match)
  - 'Sequence'      -> 'sequence'    (otu_clustering.py uses exact match)
  - 'Species'       -> 'species'     (gap_analysis.py uses exact lowercase match)
  - 'Sequence aligned - rbcl only - aligned 04/09/2025' -> dropped
    (aligned column is only populated for rbcl records; we want the unaligned sequence)

All other columns are preserved unchanged, including 'Original taxon name'
for reference.

Author: Ben Price / Claude
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path


# Column name in the raw diatom TSV -> column name expected by the pipeline
RENAMES = {
    'Sequence ID': 'accession',
    'Sequence': 'sequence',
    'Species': 'species',
}

# Columns to drop from the output
DROP_COLUMNS = {
    'Sequence aligned - rbcl only - aligned 04/09/2025',
}


def prep_diatom_file(input_path: Path, output_path: Path, verbose: bool = False) -> None:
    """
    Reformat the diatom TSV for pipeline compatibility.

    Args:
        input_path: Path to the raw diatom TSV file
        output_path: Path to write the reformatted TSV
        verbose: Print progress information
    """
    with open(input_path, 'r', encoding='utf-8') as fin:
        # --- Header ---
        header_line = fin.readline().rstrip('\n\r')
        original_headers = header_line.split('\t')

        # Validate expected columns exist
        missing_renames = [col for col in RENAMES if col not in original_headers]
        if missing_renames:
            print(f"ERROR: Expected column(s) not found in input: {missing_renames}",
                  file=sys.stderr)
            print(f"Available columns: {original_headers}", file=sys.stderr)
            sys.exit(1)

        missing_drops = [col for col in DROP_COLUMNS if col not in original_headers]
        if missing_drops and verbose:
            print(f"WARNING: Column(s) to drop were not found (may have been renamed): "
                  f"{missing_drops}", file=sys.stderr)

        # Build output header: apply renames, skip dropped columns
        output_headers = []
        col_indices = []  # indices of columns to include, in order
        for idx, col in enumerate(original_headers):
            if col in DROP_COLUMNS:
                if verbose:
                    print(f"  Dropping column: '{col}'", file=sys.stderr)
                continue
            new_name = RENAMES.get(col, col)
            output_headers.append(new_name)
            col_indices.append(idx)
            if col != new_name and verbose:
                print(f"  Renaming column: '{col}' -> '{new_name}'", file=sys.stderr)

        if verbose:
            print(f"Input columns:  {len(original_headers)}", file=sys.stderr)
            print(f"Output columns: {len(output_headers)}", file=sys.stderr)

        with open(output_path, 'w', encoding='utf-8') as fout:
            fout.write('\t'.join(output_headers) + '\n')

            row_count = 0
            for line in fin:
                row_count += 1
                line = line.rstrip('\n\r')
                fields = line.split('\t')

                # Select only the kept columns (pad if row is short)
                out_fields = []
                for idx in col_indices:
                    out_fields.append(fields[idx] if idx < len(fields) else '')

                fout.write('\t'.join(out_fields) + '\n')

                if verbose and row_count % 10000 == 0:
                    print(f"  Processed {row_count} rows...", file=sys.stderr)

    if verbose:
        print(f"Done. {row_count} data rows written to: {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Prepare diatom TSV dataset for the mind-the-gap pipeline.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Column transformations:
  'Sequence ID'                                    -> 'accession'
  'Sequence'                                       -> 'sequence'
  'Species'                                        -> 'species'
  'Sequence aligned - rbcl only - ...'             -> (dropped)

The output is compatible with:
  otu_clustering.py   (uses 'accession' and 'sequence' columns)
  bags_assessment.py  (uses 'accession', 'species', 'OTU_ID' columns)
  gap_analysis.py     (uses 'species' and 'OTU_ID' columns)

Example:
  python prep_diatom_data.py 2026-02-20_diatoms.tsv 2026-02-20_diatoms_prepped.tsv
  python prep_diatom_data.py -v 2026-02-20_diatoms.tsv 2026-02-20_diatoms_prepped.tsv
'''
    )
    parser.add_argument('input', type=str,
                        help='Input diatom TSV file')
    parser.add_argument('output', type=str,
                        help='Output reformatted TSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress information to stderr')

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        print(f"Input:  {input_path}", file=sys.stderr)
        print(f"Output: {output_path}", file=sys.stderr)

    prep_diatom_file(input_path, output_path, args.verbose)


if __name__ == '__main__':
    main()
