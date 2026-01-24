#!/usr/bin/env python3
"""
JNCC Conservation Designations to UKSI TVK Mapper

Maps JNCC conservation designation taxa to UKSI taxon version keys (TVKs),
expanding to include related taxa based on taxonomic rank:
- Above species: includes all descendant TVKs
- Below species: includes parent species TVK
- Species level: includes self and recommended TVK

Author: Generated for NHM London biodiversity genomics work
"""

import argparse
import csv
import logging
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Set, Optional, Tuple

# Rank classification lookup
RANK_CLASSIFICATION = {
    # Above species
    'Genus': 'above', 'Species aggregate': 'above', 'Order': 'above',
    'Infraclass': 'above', 'Subclass': 'above', 'Family': 'above',
    'Suborder': 'above', 'Superorder': 'above', 'Superfamily': 'above',
    'Subgenus': 'above', 'Cohort': 'above', 'Supercohort': 'above',
    'Subcohort': 'above', 'Subfamily': 'above', 'Class': 'above',
    'Phylum': 'above', 'Genus aggregate': 'above', 'Division': 'above',
    'Tribe': 'above', 'Species group': 'above', 'Ordinal aggregate': 'above',
    'Family aggregate': 'above', 'Subphylum': 'above', 'Subdivision': 'above',
    'Functional group': 'above', 'Facies': 'above', 'Superclass': 'above',
    'Species sensu lato': 'above', 'Infraorder': 'above', 'Subkingdom': 'above',
    'Pathovar': 'above', 'Generic hybrid': 'above', 'Section': 'above',
    'Infrakingdom': 'above', 'ecad': 'above', 'Infraphylum': 'above',
    'Parvorder': 'above', 'Subsection': 'above', 'Subtribe': 'above',
    'Sublusus': 'above', 'Series': 'above', 'Subseries': 'above',
    'Supertribe': 'above', 'Subterclass': 'above', 'Kingdom': 'above',
    'Domain': 'above',
    
    # Species level
    'Species': 'species', 'Unranked': 'species', 'Species hybrid': 'species',
    
    # Below species
    'Subspecies': 'below', 'Variety': 'below', 'Form': 'below',
    'Species pro parte': 'below', 'Breed': 'below', 'Subspecies aggregate': 'below',
    'Race': 'below', 'Subspecies hybrid': 'below', 'Cultivar': 'below',
    'Morphotype': 'below', 'Nothosubspecies': 'below', 'Varietal hybrid': 'below',
    'Praespecies': 'below', 'Subvariety': 'below', 'Microspecies': 'below',
    'Convariety': 'below', 'Microgene': 'below', 'Nothovariety': 'below',
    'Forma specialis': 'below', 'Abberation': 'below', 'Subform': 'below',
    
    # Ambiguous
    'Unknown': 'unknown'
}

SPECIES_RANKS = {'Species', 'Unranked', 'Species hybrid'}


def setup_logging(log_file: Path) -> logging.Logger:
    """Configure logging to both file and console."""
    logger = logging.getLogger('jncc_mapper')
    logger.setLevel(logging.DEBUG)
    
    # File handler - detailed logging
    fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    
    # Console handler - info and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger


def read_tsv_with_fallback(filepath: Path, logger: logging.Logger) -> list:
    """Read TSV file with encoding fallback."""
    encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'windows-1252']
    
    for encoding in encodings:
        try:
            with open(filepath, 'r', encoding=encoding, newline='') as f:
                reader = csv.DictReader(f, delimiter='\t')
                data = list(reader)
                logger.info(f"Successfully read {filepath.name} with {encoding} encoding ({len(data)} rows)")
                return data
        except UnicodeDecodeError:
            continue
        except Exception as e:
            logger.error(f"Error reading {filepath}: {e}")
            raise
    
    raise ValueError(f"Could not read {filepath} with any supported encoding")


class TaxonomyTree:
    """Builds and queries taxonomic hierarchy from UKSI data."""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        # TVK -> parent TVK mapping
        self.parent_map: Dict[str, Optional[str]] = {}
        # TVK -> set of child TVKs
        self.children_map: Dict[str, Set[str]] = defaultdict(set)
        # TVK -> rank
        self.rank_map: Dict[str, str] = {}
        # TVK -> recommended TVK (from names file)
        self.recommended_map: Dict[str, str] = {}
        # All known TVKs
        self.all_tvks: Set[str] = set()
    
    def load_taxa(self, taxa_data: list):
        """Load taxonomic hierarchy from taxa file."""
        for row in taxa_data:
            tvk = row.get('TAXON_VERSION_KEY', '').strip()
            parent_tvk = row.get('PARENT_TVK', '').strip()
            rank = row.get('RANK', '').strip()
            
            if not tvk:
                continue
            
            self.all_tvks.add(tvk)
            self.rank_map[tvk] = rank
            
            if parent_tvk:
                self.parent_map[tvk] = parent_tvk
                self.children_map[parent_tvk].add(tvk)
            else:
                self.parent_map[tvk] = None
        
        self.logger.info(f"Loaded {len(self.all_tvks)} taxa into hierarchy")

    
    def load_names(self, names_data: list):
        """Load TVK to recommended TVK mappings from names file."""
        for row in names_data:
            tvk = row.get('TAXON_VERSION_KEY', '').strip()
            rec_tvk = row.get('RECOMMENDED_TAXON_VERSION_KEY', '').strip()
            
            if tvk:
                self.all_tvks.add(tvk)
                if rec_tvk:
                    self.recommended_map[tvk] = rec_tvk
                    self.all_tvks.add(rec_tvk)
        
        self.logger.info(f"Loaded {len(self.recommended_map)} TVK->recommended mappings")
    
    def get_rank_category(self, tvk: str) -> str:
        """Get rank category (above/species/below/unknown) for a TVK."""
        rank = self.rank_map.get(tvk, 'Unknown')
        return RANK_CLASSIFICATION.get(rank, 'unknown')
    
    def infer_rank_category(self, tvk: str) -> str:
        """
        Infer rank category from hierarchy position when rank is ambiguous.
        If parent is species-level, treat as below species.
        If any child is species-level, treat as above species.
        """
        category = self.get_rank_category(tvk)
        
        if category != 'unknown':
            return category
        
        # Check parent
        parent = self.parent_map.get(tvk)
        if parent:
            parent_rank = self.rank_map.get(parent, '')
            if parent_rank in SPECIES_RANKS:
                self.logger.debug(f"Inferred {tvk} as 'below' species (parent {parent} is species)")
                return 'below'

        
        # Check children
        children = self.children_map.get(tvk, set())
        for child in children:
            child_rank = self.rank_map.get(child, '')
            if child_rank in SPECIES_RANKS:
                self.logger.debug(f"Inferred {tvk} as 'above' species (child {child} is species)")
                return 'above'
        
        # Check if any descendant is species level
        if self._has_species_descendant(tvk):
            self.logger.debug(f"Inferred {tvk} as 'above' species (has species descendant)")
            return 'above'
        
        # Check if any ancestor is species level
        if self._has_species_ancestor(tvk):
            self.logger.debug(f"Inferred {tvk} as 'below' species (has species ancestor)")
            return 'below'
        
        self.logger.warning(f"Could not infer rank category for {tvk}, treating as species")
        return 'species'
    
    def _has_species_descendant(self, tvk: str, visited: Set[str] = None) -> bool:
        """Check if any descendant is at species level."""
        if visited is None:
            visited = set()
        
        if tvk in visited:
            return False
        visited.add(tvk)
        
        for child in self.children_map.get(tvk, set()):
            child_rank = self.rank_map.get(child, '')
            if child_rank in SPECIES_RANKS:
                return True
            if self._has_species_descendant(child, visited):
                return True
        return False

    
    def _has_species_ancestor(self, tvk: str, visited: Set[str] = None) -> bool:
        """Check if any ancestor is at species level."""
        if visited is None:
            visited = set()
        
        if tvk in visited:
            return False
        visited.add(tvk)
        
        parent = self.parent_map.get(tvk)
        if not parent:
            return False
        
        parent_rank = self.rank_map.get(parent, '')
        if parent_rank in SPECIES_RANKS:
            return True
        
        return self._has_species_ancestor(parent, visited)
    
    def get_all_descendants(self, tvk: str, visited: Set[str] = None) -> Set[str]:
        """Get all descendant TVKs (recursive)."""
        if visited is None:
            visited = set()
        
        if tvk in visited:
            return set()
        visited.add(tvk)
        
        descendants = set()
        for child in self.children_map.get(tvk, set()):
            descendants.add(child)
            descendants.update(self.get_all_descendants(child, visited))
        
        return descendants

    
    def get_parent_species(self, tvk: str, visited: Set[str] = None) -> Optional[str]:
        """Walk up hierarchy to find parent at species level."""
        if visited is None:
            visited = set()
        
        if tvk in visited:
            return None
        visited.add(tvk)
        
        parent = self.parent_map.get(tvk)
        if not parent:
            return None
        
        parent_rank = self.rank_map.get(parent, '')
        if parent_rank in SPECIES_RANKS:
            return parent
        
        # Continue walking up
        return self.get_parent_species(parent, visited)
    
    def get_included_tvks(self, tvk: str) -> Tuple[Set[str], str]:
        """
        Get all TVKs that should be included for a given taxon.
        Returns (set of TVKs, rank_category used).
        """
        included = set()
        
        # Always include the original TVK
        included.add(tvk)
        
        # Add recommended TVK if different
        rec_tvk = self.recommended_map.get(tvk)
        if rec_tvk and rec_tvk != tvk:
            included.add(rec_tvk)

        
        # Determine rank category (with inference for ambiguous cases)
        rank_category = self.infer_rank_category(tvk)
        
        if rank_category == 'above':
            # Include all descendants
            descendants = self.get_all_descendants(tvk)
            included.update(descendants)
            # Also add recommended TVKs for all descendants
            for desc in descendants:
                desc_rec = self.recommended_map.get(desc)
                if desc_rec:
                    included.add(desc_rec)
        
        elif rank_category == 'below':
            # Include parent species
            parent_species = self.get_parent_species(tvk)
            if parent_species:
                included.add(parent_species)
                # Add recommended for parent species
                parent_rec = self.recommended_map.get(parent_species)
                if parent_rec:
                    included.add(parent_rec)
        
        # For 'species' category, we just have self and recommended (already added)
        
        return included, rank_category


def process_jncc_file(
    jncc_data: list,
    tree: TaxonomyTree,
    logger: logging.Logger
) -> Tuple[list, Dict]:
    """Process JNCC data and add included_tvk_list column."""
    
    results = []
    stats = {
        'total': 0,
        'matched': 0,
        'unmatched': 0,
        'above_species': 0,
        'species': 0,
        'below_species': 0
    }

    
    for row in jncc_data:
        stats['total'] += 1
        
        # Get the JNCC taxon version key
        jncc_tvk = row.get('Recommended_taxon_version', '').strip()
        
        # Create output row with all original columns
        out_row = dict(row)
        
        if not jncc_tvk:
            logger.warning(f"Row {stats['total']}: Empty Recommended_taxon_version")
            out_row['included_tvk_list'] = ''
            stats['unmatched'] += 1
            results.append(out_row)
            continue
        
        # Check if TVK exists in our data
        if jncc_tvk not in tree.all_tvks:
            logger.warning(f"Row {stats['total']}: TVK {jncc_tvk} not found in UKSI data")
            out_row['included_tvk_list'] = ''
            stats['unmatched'] += 1
            results.append(out_row)
            continue
        
        # Get included TVKs
        included_tvks, rank_category = tree.get_included_tvks(jncc_tvk)
        
        stats['matched'] += 1
        if rank_category == 'above':
            stats['above_species'] += 1
        elif rank_category == 'below':
            stats['below_species'] += 1
        else:
            stats['species'] += 1
        
        # Format as semicolon-separated list
        out_row['included_tvk_list'] = ';'.join(sorted(included_tvks))
        
        logger.debug(
            f"Row {stats['total']}: {jncc_tvk} ({rank_category}) -> "
            f"{len(included_tvks)} TVKs"
        )
        
        results.append(out_row)
    
    return results, stats



def write_output(
    results: list,
    output_path: Path,
    original_fieldnames: list,
    logger: logging.Logger
):
    """Write results to TSV file."""
    fieldnames = original_fieldnames + ['included_tvk_list']
    
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    logger.info(f"Wrote {len(results)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Map JNCC conservation designations to UKSI TVKs'
    )
    parser.add_argument(
        '--jncc', '-j',
        default=r'C:\GitHub\mind-the-gap\uksi_processing\20231206_jncc_conservation_designations_taxon.tsv',
        help='Path to JNCC designations TSV file'
    )
    parser.add_argument(
        '--names', '-n',
        default=r'C:\GitHub\mind-the-gap\uksi_processing\uksi_20251203a_input_names.tsv',
        help='Path to UKSI names TSV file'
    )
    parser.add_argument(
        '--taxa', '-t',
        default=r'C:\GitHub\mind-the-gap\uksi_processing\uksi_20251203a_input_taxa.tsv',
        help='Path to UKSI taxa hierarchy TSV file'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default=r'C:\GitHub\mind-the-gap\uksi_processing\jncc_mapping',
        help='Output directory for results'
    )
    
    args = parser.parse_args()

    
    # Setup paths
    jncc_path = Path(args.jncc)
    names_path = Path(args.names)
    taxa_path = Path(args.taxa)
    output_dir = Path(args.output_dir)
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate output filenames
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_filename = f"{jncc_path.stem}_uksi_mapped.tsv"
    output_path = output_dir / output_filename
    log_path = output_dir / f"jncc_mapping_{timestamp}.log"
    
    # Setup logging
    logger = setup_logging(log_path)
    logger.info("=" * 60)
    logger.info("JNCC Conservation Designations to UKSI TVK Mapper")
    logger.info("=" * 60)
    logger.info(f"JNCC file: {jncc_path}")
    logger.info(f"Names file: {names_path}")
    logger.info(f"Taxa file: {taxa_path}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"Log file: {log_path}")
    
    try:
        # Load data files
        logger.info("Loading input files...")
        jncc_data = read_tsv_with_fallback(jncc_path, logger)
        names_data = read_tsv_with_fallback(names_path, logger)
        taxa_data = read_tsv_with_fallback(taxa_path, logger)
        
        # Build taxonomy tree
        logger.info("Building taxonomy tree...")
        tree = TaxonomyTree(logger)
        tree.load_taxa(taxa_data)
        tree.load_names(names_data)

        
        # Process JNCC file
        logger.info("Processing JNCC designations...")
        original_fieldnames = list(jncc_data[0].keys()) if jncc_data else []
        results, stats = process_jncc_file(jncc_data, tree, logger)
        
        # Write output
        logger.info("Writing output file...")
        write_output(results, output_path, original_fieldnames, logger)
        
        # Report statistics
        logger.info("=" * 60)
        logger.info("Processing complete!")
        logger.info(f"  Total rows processed: {stats['total']}")
        logger.info(f"  Matched: {stats['matched']}")
        logger.info(f"  Unmatched: {stats['unmatched']}")
        logger.info(f"  Above species (expanded to descendants): {stats['above_species']}")
        logger.info(f"  At species level: {stats['species']}")
        logger.info(f"  Below species (expanded to parent): {stats['below_species']}")
        logger.info("=" * 60)
        
        print(f"\nOutput written to: {output_path}")
        print(f"Log written to: {log_path}")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
