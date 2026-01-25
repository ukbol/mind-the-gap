<div align="center">
  <img src="docs/boromir-gap-analysis.jpeg" alt="Mind the Gap Logo">
</div>

# [mind-the-gap](https://youtu.be/QExoX4ls9OM?si=EmPShkAIIfieSnqX)

UK species gap filling analysis and associated tools.

## Overview

This repository contains a collection of Python scripts for processing DNA barcode reference libraries, taxonomic data, and conservation designations for UK species gap analysis. The tools support workflows for:

- Processing reference sequences from BOLD, NCBI GenBank, MIDORI, and UNITE databases
- OTU (Operational Taxonomic Unit) clustering
- BAGS (Barcode, Audit & Grade System) quality assessment
- UKSI (UK Species Inventory) taxonomy processing
- JNCC conservation designation mapping

## Directory Structure

```
mind-the-gap/
├── bold_processing/          # BOLD database processing
│   └── bold_gene_extract/    # Extract sequences by marker gene
├── ncbi_processing/          # NCBI GenBank processing
│   ├── ncbi_fetch/           # Fetch sequences from NCBI
│   └── ncbi_gb_extract/      # Parse GenBank flat files
├── midori_processing/        # MIDORI reference library processing
├── unite_processing/         # UNITE fungal database processing
├── otu_clustering/           # VSEARCH-based OTU clustering
├── bags_assessment/          # BAGS quality grading
└── uksi_processing/          # UK Species Inventory processing
    ├── jncc_mapping/         # Map JNCC taxa to UKSI TVKs
    ├── jncc_annotation/      # Annotate species with JNCC designations
    ├── pantheon_mapping/     # Match Pantheon data to UKSI
    └── uksi_db/              # SQLite database pipeline
```

## Scripts

### Reference Library Processing

| Script | Location | Description |
|--------|----------|-------------|
| `bold_gene_extract.py` | `bold_processing/bold_gene_extract/` | Filters BOLD TSV data to extract rows matching specified marker genes (e.g., COI-5P, rbcL). Supports case-insensitive matching and pipeline integration. |
| `ncbi_gb_extract.py` | `ncbi_processing/ncbi_gb_extract/` | Extracts gene data from NCBI GenBank flat files using BioPython. Outputs structured TSV with taxonomy, sequences, and feature qualifiers. |
| `process_midori.py` | `midori_processing/` | Converts MIDORI FASTA reference library files to TSV format with parsed taxonomic information. |
| `process_unite.py` | `unite_processing/` | Converts UNITE fungal FASTA files to TSV format with parsed taxonomy using prefix notation (k__, p__, etc.). |

### Sequence Analysis

| Script | Location | Description |
|--------|----------|-------------|
| `otu_clustering.py` | `otu_clustering/` | Clusters DNA sequences into OTUs using VSEARCH. Accepts TSV input and appends OTU_ID column to output. Supports configurable identity thresholds and strand handling. |
| `bags_assessment.py` | `bags_assessment/` | Calculates BAGS grades (A-F) for taxonomic data based on OTU clustering results. Assesses quality of species-level barcode reference libraries. |

### UKSI Processing

| Script | Location | Description |
|--------|----------|-------------|
| `uksi_species_pipeline.py` | `uksi_processing/` | Extracts valid species from UKSI TAXA and NAMES tables with higher taxonomy, synonyms, and comprehensive validation logging. |
| `jncc_uksi_mapper.py` | `uksi_processing/jncc_mapping/` | Maps JNCC conservation designation taxa to UKSI TVKs with hierarchical expansion (descendants for higher ranks, parents for subspecies). |
| `pantheon_uksi_matcher.py` | `uksi_processing/pantheon_mapping/` | Matches Pantheon invertebrate database taxa against UKSI to retrieve TVKs for downstream processing. |
| `uksi_jncc_annotation.py` | `uksi_processing/jncc_annotation/` | Annotates UKSI species with JNCC conservation designations by matching TVKs. |
| `uksi_jncc_annotation_v2.py` | `uksi_processing/jncc_annotation/` | Enhanced annotation using `included_tvk_list` column for flexible TVK matching. |
| `uksi_import.py` | `uksi_processing/uksi_db/` | Creates SQLite database linking UKSI names/taxa with Pantheon and JNCC data. |
| `uksi_export.py` | `uksi_processing/uksi_db/` | Exports valid species from SQLite with full taxonomy, synonyms, and conservation designations with hierarchical propagation. |

## Requirements

- Python 3.7+
- BioPython (for NCBI GenBank parsing)
- VSEARCH (for OTU clustering)
- pandas (for data manipulation)

## Usage

Each script includes detailed help via `--help` flag:

```bash
python bold_processing/bold_gene_extract/bold_gene_extract.py --help
python otu_clustering/otu_clustering.py --help
python uksi_processing/uksi_db/uksi_export.py --help
```

Most scripts support stdin/stdout for pipeline integration:

```bash
cat sequences.tsv | python bold_gene_extract.py -g COI-5P | python otu_clustering.py > clustered.tsv
```

See individual README files in each directory for detailed documentation.

## BAGS Grading System

The BAGS assessment assigns quality grades based on OTU clustering results:

| Grade | Criteria |
|-------|----------|
| A | Single species in single OTU with 11+ sequences |
| B | Single species in single OTU with 3-10 sequences |
| C | Species split across multiple OTUs |
| D | Single species in single OTU with <3 sequences |
| E | Multiple species sharing same OTU (BIN sharing) |
| F | No OTU assignment |

## References

- [UKSI Nameserver](https://www.nhm.ac.uk/ukspecies/) - UK Species Inventory
- [JNCC Conservation Designations](https://jncc.gov.uk/our-work/conservation-designations-for-uk-taxa/)
- [BOLD Systems](https://www.boldsystems.org/) - Barcode of Life Data System
- [MIDORI](http://www.reference-midori.info/) - Metazoan mitochondrial reference database
- [UNITE](https://unite.ut.ee/) - Fungal ITS reference database
- [Pantheon](https://www.brc.ac.uk/pantheon/) - Biological Records Centre invertebrate database


