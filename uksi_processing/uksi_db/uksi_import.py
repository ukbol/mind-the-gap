#!/usr/bin/env python3
"""
UKSI Database Import Script
===========================
Creates a SQLite database linking UKSI names and taxa tables,
with additional Pantheon, JNCC, and freshwater species list data.

Author: Generated for Ben Price, NHM London
Date: 2025-01-25

Database Schema:
- names: All name variants from UKSI nameserver
- taxa: Backbone taxonomy with hierarchical structure
- pantheon: Invertebrate ecological traits (linked via RECOMMENDED_TAXON_VERSION_KEY)
- jncc: Conservation designations (linked via Recommended_taxon_version)
- freshbase: FreshBase freshwater species list (linked via TAXON_VERSION_KEY)
- ukceh_freshwater: UKCEH freshwater species list (linked via TAXON_VERSION_KEY)
"""

import sqlite3
import csv
import os
import sys
from datetime import datetime
from pathlib import Path

# Configuration
BASE_DIR = Path(r"C:\GitHub\mind-the-gap\uksi_processing")
DB_PATH = BASE_DIR / "uksi_db" / "uksi.db"

INPUT_FILES = {
    "names": BASE_DIR / "uksi_20251203a_input_names.tsv",
    "taxa": BASE_DIR / "uksi_20251203a_input_taxa.tsv",
    "pantheon": BASE_DIR / "pantheon_mapping" / "output" / "pantheon_input_cleaned_matched.tsv",
    "jncc": BASE_DIR / "jncc_mapping" / "20231206_jncc_conservation_designations_taxon.tsv",
    "freshbase": BASE_DIR / "freshwater" / "2026-02-10_freshbase.tsv",
    "ukceh_freshwater": BASE_DIR / "freshwater" / "UKCEH_freshwater_list.tsv",
}

# Logging setup
LOG_PATH = BASE_DIR / "uksi_db" / "uksi_import.log"

def log(message: str):
    """Log message to both console and file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] {message}"
    print(log_line)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(log_line + "\n")


def create_database_schema(conn: sqlite3.Connection):
    """Create the database schema with all tables and indexes."""
    log("Creating database schema...")
    
    cursor = conn.cursor()
    
    # Names table - all name variants from UKSI nameserver
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS names (
            INFORMAL_GROUP TEXT,
            TAXON_VERSION_KEY TEXT PRIMARY KEY,
            TAXON_NAME TEXT,
            TAXON_AUTHORITY TEXT,
            TAXON_QUALIFIER TEXT,
            LANGUAGE TEXT,
            RANK TEXT,
            NAME_FORM TEXT,
            NAME_STATUS TEXT,
            NAME_TYPE TEXT,
            DEPRECATED_DATE TEXT,
            RECOMMENDED_SCIENTIFIC_NAME TEXT,
            RECOMMENDED_NAME_AUTHORITY TEXT,
            RECOMMENDED_NAME_QUALIFIER TEXT,
            RECOMMENDED_NAME_RANK TEXT,
            RECOMMENDED_TAXON_VERSION_KEY TEXT,
            DATE_RECORD_ADDED TEXT,
            DATE_RECORD_LAST_CHANGED TEXT
        )
    """)
    
    # Taxa table - backbone taxonomy with hierarchy
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS taxa (
            LINEAGE TEXT,
            SORT_ORDER TEXT,
            ORGANISM_KEY TEXT,
            PARENT_KEY TEXT,
            TAXON_VERSION_KEY TEXT PRIMARY KEY,
            PARENT_TVK TEXT,
            INFORMAL_GROUP TEXT,
            RANK TEXT,
            TAXON_NAME TEXT,
            TAXON_AUTHORITY TEXT,
            TAXON_QUALIFIER TEXT,
            REDUNDANT_FLAG TEXT,
            NON_NATIVE_FLAG TEXT,
            TERRESTRIAL_FRESHWATER_FLAG TEXT,
            FRESHWATER TEXT,
            MARINE_FLAG TEXT
        )
    """)
    
    # Pantheon table - invertebrate ecological traits
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pantheon (
            Taxon TEXT,
            "Row" INTEGER,
            original_taxon TEXT,
            Vernacular TEXT,
            Family TEXT,
            "Order" TEXT,
            SQS TEXT,
            "Conservation status" TEXT,
            "Current conservation status" TEXT,
            "Designation summary" TEXT,
            "Larval feeding guild" TEXT,
            "Adult feeding guild" TEXT,
            "Broad biotope" TEXT,
            Habitat TEXT,
            Resources TEXT,
            "Specific assemblage type" TEXT,
            "Habitat score" TEXT,
            Associations TEXT,
            Genus TEXT,
            subgenus TEXT,
            species TEXT,
            subspecies TEXT,
            TAXON_VERSION_KEY TEXT,
            RECOMMENDED_TAXON_VERSION_KEY TEXT,
            NAME_STATUS TEXT,
            MATCHED_ON TEXT
        )
    """)
    
    # JNCC table - conservation designations
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS jncc (
            informal_group TEXT,
            Recommended_taxon_name TEXT,
            Recommended_authority TEXT,
            Recommended_qualifier TEXT,
            Recommended_taxon_version TEXT,
            "A: Bern Convention" TEXT,
            "C: Birds Directive" TEXT,
            "C1: Convention on Migratory Species" TEXT,
            "C2: OSPAR" TEXT,
            "D: Habitats Directive" TEXT,
            "E: EC Cites" TEXT,
            "F: Global Red list status" TEXT,
            "Fa: Red Listing based on pre 1994 IUCN guidelines" TEXT,
            "Fb: Red Listing based on 1994 IUCN guidelines" TEXT,
            "Fc: Red listing based on 2001 IUCN guidelines" TEXT,
            "Fd: Red data categories  - birds (not based on IUCN criteria)" TEXT,
            "Fe: Red data categories - Spiders (not based on IUCN criteria)" TEXT,
            "Ga: Rare and scarce species" TEXT,
            "Gb: Rare and scarce species (not based on IUCN criteria)" TEXT,
            "Ha: Biodiversity Action Plan UK list of priority species" TEXT,
            "Hb: Biodiversity Lists - England" TEXT,
            "Hc: Biodiversity Lists - Scotland" TEXT,
            "Hd: Biodiversity Lists - Wales" TEXT,
            "He: Biodiversity Lists - Northern Ireland" TEXT,
            "I: Wildlife and Countryside Act 1981" TEXT,
            "J: The Wildlife (Northern Ireland) Order 1985" TEXT,
            "K: The Conservation of Habitats and Species Regulations 2010" TEXT,
            "L: The Conservation (Nature Habitats, etc_) Regulations (NI) 199" TEXT,
            "M: Protection of Badgers Act" TEXT
        )
    """)
    
    # FreshBase table - freshwater species list
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS freshbase (
            ORGANISM_KEY TEXT,
            TAXON_VERSION_KEY TEXT,
            Phylum TEXT,
            Subphylum TEXT,
            "Class" TEXT,
            "Order" TEXT,
            Family TEXT,
            Genus TEXT,
            species TEXT,
            Taxon TEXT
        )
    """)

    # UKCEH freshwater table - UKCEH freshwater species list
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ukceh_freshwater (
            TAXON_GROUP_NAME TEXT,
            TAXON_LIST_KEY TEXT,
            TAXON_LIST_VERSION_KEY TEXT,
            TAXON_LIST_ITEM_KEY TEXT,
            PARENT TEXT,
            TAXON_VERSION_KEY TEXT,
            PREFERRED_NAME TEXT,
            SORT_CODE TEXT,
            LST_ITM_CODE TEXT,
            TAXON_RANK_KEY TEXT,
            LONG_NAME TEXT,
            ITEM_NAME TEXT,
            AUTHORITY TEXT,
            ATTRIBUTE TEXT,
            LANGUAGE TEXT,
            NOTE TEXT,
            CHANGED_BY TEXT,
            CHANGED_DATE TEXT,
            DELETED_DATE TEXT
        )
    """)

    # Create indexes for efficient lookups
    log("Creating indexes...")
    
    # Names table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_names_recommended_tvk ON names(RECOMMENDED_TAXON_VERSION_KEY)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_names_language ON names(LANGUAGE)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_names_rank ON names(RANK)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_names_status ON names(NAME_STATUS)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_names_form ON names(NAME_FORM)")
    
    # Taxa table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxa_parent_tvk ON taxa(PARENT_TVK)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxa_organism_key ON taxa(ORGANISM_KEY)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxa_parent_key ON taxa(PARENT_KEY)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxa_rank ON taxa(RANK)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxa_lineage ON taxa(LINEAGE)")
    
    # Pantheon table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_pantheon_rec_tvk ON pantheon(RECOMMENDED_TAXON_VERSION_KEY)")
    
    # JNCC table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_jncc_rec_tvk ON jncc(Recommended_taxon_version)")

    # FreshBase table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_freshbase_tvk ON freshbase(TAXON_VERSION_KEY)")

    # UKCEH freshwater table indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ukceh_freshwater_tvk ON ukceh_freshwater(TAXON_VERSION_KEY)")

    conn.commit()
    log("Schema created successfully")


def import_tsv_to_table(conn: sqlite3.Connection, file_path: Path, table_name: str, 
                         skip_header_check: bool = False):
    """
    Import a TSV file into a SQLite table.
    
    Args:
        conn: SQLite connection
        file_path: Path to the TSV file
        table_name: Name of the target table
        skip_header_check: If True, skip the first column name check (for JNCC file)
    """
    log(f"Importing {file_path.name} into {table_name}...")
    
    if not file_path.exists():
        log(f"ERROR: File not found: {file_path}")
        return 0
    
    cursor = conn.cursor()
    row_count = 0
    error_count = 0
    
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        
        # Clean headers - remove BOM if present
        if headers[0].startswith('\ufeff'):
            headers[0] = headers[0][1:]
        
        # For JNCC file, first column is oddly named - rename it
        if table_name == "jncc":
            headers[0] = "informal_group"
        
        # Build the INSERT statement
        placeholders = ','.join(['?' for _ in headers])
        quoted_headers = [f'"{h}"' for h in headers]
        insert_sql = f"INSERT OR REPLACE INTO {table_name} ({','.join(quoted_headers)}) VALUES ({placeholders})"
        
        # Batch insert for performance
        batch = []
        batch_size = 10000
        
        for row in reader:
            # Pad row with empty strings if shorter than headers
            while len(row) < len(headers):
                row.append('')
            # Truncate if longer
            row = row[:len(headers)]
            
            batch.append(row)
            row_count += 1
            
            if len(batch) >= batch_size:
                try:
                    cursor.executemany(insert_sql, batch)
                    conn.commit()
                except sqlite3.Error as e:
                    error_count += len(batch)
                    log(f"  Batch error: {e}")
                batch = []
                
                if row_count % 50000 == 0:
                    log(f"  Processed {row_count:,} rows...")
        
        # Insert remaining rows
        if batch:
            try:
                cursor.executemany(insert_sql, batch)
                conn.commit()
            except sqlite3.Error as e:
                error_count += len(batch)
                log(f"  Final batch error: {e}")
    
    log(f"  Imported {row_count:,} rows ({error_count} errors)")
    return row_count


def validate_import(conn: sqlite3.Connection):
    """Validate the imported data with summary statistics."""
    log("\n=== Import Validation ===")
    cursor = conn.cursor()
    
    # Names table stats
    cursor.execute("SELECT COUNT(*) FROM names")
    names_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT RECOMMENDED_TAXON_VERSION_KEY) FROM names")
    unique_recommended = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM names WHERE NAME_STATUS = 'R'")
    recommended_names = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM names WHERE LANGUAGE = 'la'")
    latin_names = cursor.fetchone()[0]
    log(f"Names table: {names_count:,} total, {unique_recommended:,} unique recommended TVKs")
    log(f"  - Recommended names (R): {recommended_names:,}")
    log(f"  - Latin names: {latin_names:,}")
    
    # Taxa table stats
    cursor.execute("SELECT COUNT(*) FROM taxa")
    taxa_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM taxa WHERE RANK = 'Species'")
    species_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT RANK) FROM taxa")
    rank_count = cursor.fetchone()[0]
    log(f"Taxa table: {taxa_count:,} total, {species_count:,} species, {rank_count} unique ranks")
    
    # Show rank distribution for valid species identification
    cursor.execute("""
        SELECT RANK, COUNT(*) as cnt 
        FROM taxa 
        GROUP BY RANK 
        ORDER BY cnt DESC 
        LIMIT 20
    """)
    log("  Top ranks in taxa:")
    for row in cursor.fetchall():
        log(f"    {row[0]}: {row[1]:,}")
    
    # Pantheon table stats
    cursor.execute("SELECT COUNT(*) FROM pantheon")
    pantheon_count = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM pantheon p
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = p.RECOMMENDED_TAXON_VERSION_KEY)
    """)
    pantheon_matched = cursor.fetchone()[0]
    log(f"Pantheon table: {pantheon_count:,} total, {pantheon_matched:,} matched to taxa")
    
    # JNCC table stats
    cursor.execute("SELECT COUNT(*) FROM jncc")
    jncc_count = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM jncc j
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = j.Recommended_taxon_version)
    """)
    jncc_direct_matched = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM jncc_resolved jr
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = jr.resolved_tvk)
    """)
    jncc_resolved_matched = cursor.fetchone()[0]
    log(f"JNCC table: {jncc_count:,} total")
    log(f"  - Direct match to taxa: {jncc_direct_matched:,}")
    log(f"  - Resolved via names table: {jncc_resolved_matched:,} ({jncc_resolved_matched - jncc_direct_matched:,} additional)")
    
    # FreshBase table stats
    cursor.execute("SELECT COUNT(*) FROM freshbase")
    freshbase_count = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM freshbase f
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = f.TAXON_VERSION_KEY)
    """)
    freshbase_direct_matched = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM freshbase_resolved fr
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = fr.resolved_tvk)
    """)
    freshbase_resolved_matched = cursor.fetchone()[0]
    log(f"FreshBase table: {freshbase_count:,} total")
    log(f"  - Direct match to taxa: {freshbase_direct_matched:,}")
    log(f"  - Resolved via names table: {freshbase_resolved_matched:,} ({freshbase_resolved_matched - freshbase_direct_matched:,} additional)")

    # UKCEH freshwater table stats
    cursor.execute("SELECT COUNT(*) FROM ukceh_freshwater")
    ukceh_count = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM ukceh_freshwater u
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = u.TAXON_VERSION_KEY)
    """)
    ukceh_direct_matched = cursor.fetchone()[0]
    cursor.execute("""
        SELECT COUNT(*) FROM ukceh_freshwater_resolved ur
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = ur.resolved_tvk)
    """)
    ukceh_resolved_matched = cursor.fetchone()[0]
    log(f"UKCEH freshwater table: {ukceh_count:,} total")
    log(f"  - Direct match to taxa: {ukceh_direct_matched:,}")
    log(f"  - Resolved via names table: {ukceh_resolved_matched:,} ({ukceh_resolved_matched - ukceh_direct_matched:,} additional)")

    # Check linkage between names and taxa via recommended TVK
    cursor.execute("""
        SELECT COUNT(DISTINCT n.RECOMMENDED_TAXON_VERSION_KEY)
        FROM names n
        WHERE EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = n.RECOMMENDED_TAXON_VERSION_KEY)
    """)
    linked_count = cursor.fetchone()[0]
    log(f"\nNames->Taxa linkage: {linked_count:,} recommended TVKs link to taxa table")


def create_utility_views(conn: sqlite3.Connection):
    """Create useful views for querying the database."""
    log("\nCreating utility views...")
    cursor = conn.cursor()
    
    # View: Valid species (species-level taxa that aren't subgenera)
    # This includes Species, Microspecies, and similar species-level ranks
    cursor.execute("""
        CREATE VIEW IF NOT EXISTS v_valid_species AS
        SELECT t.*
        FROM taxa t
        WHERE t.RANK IN (
            'Species', 'Microspecies', 'Species hybrid', 'Species aggregate',
            'Intergeneric hybrid', 'Species sensu lato', 'Species sensu stricto'
        )
        AND t.REDUNDANT_FLAG IS NULL OR t.REDUNDANT_FLAG = ''
    """)
    
    # View: All Latin synonyms for each recommended TVK
    cursor.execute("""
        CREATE VIEW IF NOT EXISTS v_latin_synonyms AS
        SELECT 
            n.RECOMMENDED_TAXON_VERSION_KEY,
            GROUP_CONCAT(
                CASE 
                    WHEN n.TAXON_AUTHORITY IS NOT NULL AND n.TAXON_AUTHORITY != '' 
                    THEN n.TAXON_NAME || ' ' || n.TAXON_AUTHORITY
                    ELSE n.TAXON_NAME
                END, 
                ' | '
            ) as all_synonyms,
            COUNT(*) as synonym_count
        FROM names n
        WHERE n.LANGUAGE = 'la'
        AND n.TAXON_VERSION_KEY != n.RECOMMENDED_TAXON_VERSION_KEY
        GROUP BY n.RECOMMENDED_TAXON_VERSION_KEY
    """)
    
    # View: Species with higher taxonomy (via lineage traversal)
    # This is a helper view that will be expanded in the export script
    cursor.execute("""
        CREATE VIEW IF NOT EXISTS v_species_with_taxonomy AS
        SELECT 
            t.TAXON_VERSION_KEY,
            t.TAXON_NAME,
            t.TAXON_AUTHORITY,
            t.RANK,
            t.LINEAGE,
            t.INFORMAL_GROUP,
            t.NON_NATIVE_FLAG,
            t.TERRESTRIAL_FRESHWATER_FLAG,
            t.FRESHWATER,
            t.MARINE_FLAG
        FROM v_valid_species t
    """)
    
    # View: JNCC resolved via names table
    # The JNCC file uses some older/synonym TVKs that need resolution through names table
    # This view resolves those to their current recommended TVK in taxa
    cursor.execute("""
        CREATE VIEW IF NOT EXISTS v_jncc_resolved AS
        SELECT 
            j.*,
            COALESCE(
                -- First try direct match to taxa
                CASE WHEN EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = j.Recommended_taxon_version)
                     THEN j.Recommended_taxon_version END,
                -- Otherwise resolve via names table
                n.RECOMMENDED_TAXON_VERSION_KEY
            ) as resolved_tvk
        FROM jncc j
        LEFT JOIN names n ON n.TAXON_VERSION_KEY = j.Recommended_taxon_version
    """)
    
    # Create index on the resolved view (as a materialized approach - create a table)
    # This is faster than computing on-the-fly for large joins
    cursor.execute("DROP TABLE IF EXISTS jncc_resolved")
    cursor.execute("""
        CREATE TABLE jncc_resolved AS
        SELECT 
            j.*,
            COALESCE(
                CASE WHEN EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = j.Recommended_taxon_version)
                     THEN j.Recommended_taxon_version END,
                n.RECOMMENDED_TAXON_VERSION_KEY
            ) as resolved_tvk
        FROM jncc j
        LEFT JOIN names n ON n.TAXON_VERSION_KEY = j.Recommended_taxon_version
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_jncc_resolved_tvk ON jncc_resolved(resolved_tvk)")

    # FreshBase resolved via names table
    # Resolves TVKs to their current recommended TVK in taxa
    cursor.execute("DROP TABLE IF EXISTS freshbase_resolved")
    cursor.execute("""
        CREATE TABLE freshbase_resolved AS
        SELECT
            f.*,
            COALESCE(
                CASE WHEN EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = f.TAXON_VERSION_KEY)
                     THEN f.TAXON_VERSION_KEY END,
                n.RECOMMENDED_TAXON_VERSION_KEY
            ) as resolved_tvk
        FROM freshbase f
        LEFT JOIN names n ON n.TAXON_VERSION_KEY = f.TAXON_VERSION_KEY
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_freshbase_resolved_tvk ON freshbase_resolved(resolved_tvk)")

    # UKCEH freshwater resolved via names table
    # Resolves TVKs to their current recommended TVK in taxa
    cursor.execute("DROP TABLE IF EXISTS ukceh_freshwater_resolved")
    cursor.execute("""
        CREATE TABLE ukceh_freshwater_resolved AS
        SELECT
            u.*,
            COALESCE(
                CASE WHEN EXISTS (SELECT 1 FROM taxa t WHERE t.TAXON_VERSION_KEY = u.TAXON_VERSION_KEY)
                     THEN u.TAXON_VERSION_KEY END,
                n.RECOMMENDED_TAXON_VERSION_KEY
            ) as resolved_tvk
        FROM ukceh_freshwater u
        LEFT JOIN names n ON n.TAXON_VERSION_KEY = u.TAXON_VERSION_KEY
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ukceh_freshwater_resolved_tvk ON ukceh_freshwater_resolved(resolved_tvk)")

    conn.commit()
    log("Utility views created")


def main():
    """Main entry point for the import script."""
    log("=" * 60)
    log("UKSI Database Import Script")
    log("=" * 60)
    
    # Ensure output directory exists
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    
    # Clear log file
    if LOG_PATH.exists():
        LOG_PATH.unlink()
    
    # Check all input files exist
    log("\nChecking input files...")
    all_exist = True
    for name, path in INPUT_FILES.items():
        exists = path.exists()
        status = "OK" if exists else "MISSING"
        log(f"  {name}: {status} ({path})")
        if not exists:
            all_exist = False
    
    if not all_exist:
        log("\nERROR: Some input files are missing. Aborting.")
        sys.exit(1)
    
    # Remove existing database to start fresh
    if DB_PATH.exists():
        log(f"\nRemoving existing database: {DB_PATH}")
        DB_PATH.unlink()
    
    # Create database connection
    log(f"\nCreating database: {DB_PATH}")
    conn = sqlite3.connect(str(DB_PATH))
    
    # Enable optimizations
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=-64000")  # 64MB cache
    
    try:
        # Create schema
        create_database_schema(conn)
        
        # Import data
        log("\n=== Importing Data ===")
        import_tsv_to_table(conn, INPUT_FILES["names"], "names")
        import_tsv_to_table(conn, INPUT_FILES["taxa"], "taxa")
        import_tsv_to_table(conn, INPUT_FILES["pantheon"], "pantheon")
        import_tsv_to_table(conn, INPUT_FILES["jncc"], "jncc")
        import_tsv_to_table(conn, INPUT_FILES["freshbase"], "freshbase")
        import_tsv_to_table(conn, INPUT_FILES["ukceh_freshwater"], "ukceh_freshwater")
        
        # Create utility views
        create_utility_views(conn)
        
        # Validate import
        validate_import(conn)
        
        # Final database optimization
        log("\nOptimizing database...")
        conn.execute("ANALYZE")
        conn.execute("VACUUM")
        
        log("\n" + "=" * 60)
        log("Import completed successfully!")
        log(f"Database saved to: {DB_PATH}")
        log("=" * 60)
        
    except Exception as e:
        log(f"\nERROR: {e}")
        import traceback
        log(traceback.format_exc())
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
