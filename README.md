<div align="center">
  <img src="docs/boromir-gap-analysis.jpeg" alt="Mind the Gap Logo">
</div>

[![DOI](https://zenodo.org/badge/1139918985.svg)](https://doi.org/10.5281/zenodo.18378731)

# [mind-the-gap](https://youtu.be/QExoX4ls9OM?si=EmPShkAIIfieSnqX)

DNA reference library gap analysis and associated tools.

## Overview

This repository contains a collection of Python scripts for processing DNA barcode reference libraries, taxonomic data, and conservation designations for UK species gap analysis. The tools support workflows for:

- Processing reference sequences from BOLD, NCBI GenBank, MIDORI, and UNITE databases
- OTU (Operational Taxonomic Unit) clustering
- BAGS (Barcode, Audit & Grade System) quality assessment
- Gap analysis with parallel processing for HPC environments
- Darwin Tree of Life (DToL) genome sequencing status matching
- UKSI (UK Species Inventory) taxonomy processing
- JNCC conservation designation mapping and Pantheon ecological trait integration

## Directory Structure

```
mind-the-gap/
├── bold_processing/          # BOLD database processing
│   └── bold_gene_extract/    # Extract sequences by marker gene
├── ncbi_processing/          # NCBI GenBank processing
│   ├── ncbi_fetch/           # Batch download GenBank records via EDirect
│   └── ncbi_gb_extract/      # Parse GenBank flat files to TSV
├── midori_processing/        # MIDORI reference library processing
├── unite_processing/         # UNITE fungal database processing
├── otu_clustering/           # VSEARCH-based OTU clustering
├── bags_assessment/          # BAGS quality grading
├── gap_analysis/             # Taxon-centric gap analysis (HPC-optimised)
├── dtol_processing/          # Darwin Tree of Life status matching
├── uksi_processing/          # UK Species Inventory processing
│   ├── uksi_db/              # SQLite database pipeline (import/export)
│   ├── jncc_mapping/         # Map JNCC taxa to UKSI TVKs
│   ├── jncc_annotation/      # Annotate species with JNCC designations
│   ├── pantheon_mapping/     # Match Pantheon data to UKSI
│   └── freshwater/           # Freshwater species reference lists
├── final_result/             # Output gap analysis results
├── archive/                  # Archived result sets
└── docs/                     # Documentation and images
```

## Scripts

### Reference Library Processing

| Script | Location | Description |
|--------|----------|-------------|
| `bold_gene_extract.py` | `bold_processing/bold_gene_extract/` | Filters BOLD TSV data to extract rows matching specified marker genes (e.g., COI-5P, rbcL). Supports case-insensitive matching and pipeline integration via stdin/stdout. |
| `edirect_genbank_fetch.sh` | `ncbi_processing/ncbi_fetch/` | Batch downloads GenBank flat files from NCBI using EDirect. Configurable query, resumable, rate-limit compliant. |
| `ncbi_gb_extract.py` | `ncbi_processing/ncbi_gb_extract/` | Extracts gene data from NCBI GenBank flat files using BioPython. Outputs structured TSV with taxonomy, sequences, and feature qualifiers. |
| `process_midori.py` | `midori_processing/` | Converts MIDORI FASTA reference library files to TSV format with parsed taxonomic hierarchy and NCBI taxids. |
| `process_unite.py` | `unite_processing/` | Converts UNITE fungal FASTA files to TSV format with parsed taxonomy using prefix notation (k__, p__, etc.) and Species Hypothesis cluster codes. |

### Sequence Analysis

| Script | Location | Description |
|--------|----------|-------------|
| `otu_clustering.py` | `otu_clustering/` | Clusters DNA sequences into OTUs using VSEARCH. Accepts TSV input and appends OTU_ID column to output. Supports configurable identity thresholds, strand handling, and parallel threading. |
| `bags_assessment.py` | `bags_assessment/` | Calculates BAGS grades (A-F) for taxonomic data based on OTU clustering results. Assesses quality of species-level barcode reference libraries. |

### Gap Analysis & Assessment

| Script | Location | Description |
|--------|----------|-------------|
| `gap_analysis.py` | `gap_analysis/` | HPC-optimised taxon-centric gap analysis. Loads a species list with synonyms, scans records for matches, detects BIN/OTU sharing, and assigns BAGS grades and traffic-light status. Supports parallel processing. |
| `sh-run_gap_analysis.sh` | `gap_analysis/` | SLURM batch job submission script for running gap_analysis.py on HPC clusters. |
| `dtol_status.py` | `dtol_processing/` | Matches Darwin Tree of Life genome sequencing metadata against a UKSI species list. Maps DToL pipeline stages to traffic-light status colours. |

### UKSI Processing

| Script | Location | Description |
|--------|----------|-------------|
| `uksi_import.py` | `uksi_processing/uksi_db/` | Creates SQLite database linking UKSI names/taxa with Pantheon and JNCC data. Resolves JNCC synonym TVKs through the names table. |
| `uksi_export.py` | `uksi_processing/uksi_db/` | Exports valid species from SQLite with full taxonomy, synonyms, and conservation designations with hierarchical propagation. |
| `jncc_uksi_mapper.py` | `uksi_processing/jncc_mapping/` | Maps JNCC conservation designation taxa to UKSI TVKs with hierarchical expansion (descendants for higher ranks, parents for subspecies). |
| `uksi_jncc_annotation.py` | `uksi_processing/jncc_annotation/` | Annotates UKSI species with JNCC conservation designations by matching TVKs. |
| `uksi_jncc_annotation_v2.py` | `uksi_processing/jncc_annotation/` | Enhanced annotation using `included_tvk_list` column for flexible TVK matching. |
| `pantheon_uksi_matcher.py` | `uksi_processing/pantheon_mapping/` | Matches Pantheon invertebrate database taxa against UKSI to retrieve TVKs for downstream processing. |

## Requirements

- Python 3.7+
- BioPython (for NCBI GenBank parsing)
- VSEARCH (for OTU clustering)
- pandas (for Pantheon matching)
- EDirect (for NCBI batch downloads)

## Usage

Each script includes detailed help via `--help` flag:

```bash
python bold_processing/bold_gene_extract/bold_gene_extract.py --help
python otu_clustering/otu_clustering.py --help
python gap_analysis/gap_analysis.py --help
python dtol_processing/dtol_status.py --help
python uksi_processing/uksi_db/uksi_export.py --help
```

Most scripts support stdin/stdout for pipeline integration:

```bash
cat sequences.tsv | python bold_gene_extract.py -g COI-5P | python otu_clustering.py -t 0.99 > clustered.tsv
```

See individual README files in each directory for detailed documentation.

## Automated runs (GitHub Actions + Zenodo)

The full pipeline can be run end-to-end on GitHub Actions, with large
inputs hosted on Zenodo and results published back to a Zenodo record.

- Edit `pipeline/inputs.yaml` to point at your Zenodo concept DOIs.
- The workflow polls Zenodo weekly, or can be fired manually via the
  *Run gap-analysis pipeline* workflow on the Actions tab.
- See [`docs/automation.md`](docs/automation.md) for the one-time
  setup (Larger Runner label, `ZENODO_TOKEN` secret, container build)
  and operator runbook.

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

## Traffic Light Status

The gap analysis assigns a traffic-light status to each taxon:

| Status | Meaning |
|--------|---------|
| GREEN | Only valid name recorded, no conflicts |
| AMBER | Both valid name and synonym(s) recorded (nomenclatural cleanup needed) |
| BLUE | Only synonym(s) recorded, valid name missing from database |
| RED | Names from outside the taxon share BIN/OTU (taxonomic conflict) |
| BLACK | No records for this taxon |

## References

- [UKSI Nameserver](https://www.nhm.ac.uk/ukspecies/) - UK Species Inventory
- [JNCC Conservation Designations](https://jncc.gov.uk/our-work/conservation-designations-for-uk-taxa/)
- [BOLD Systems](https://www.boldsystems.org/) - Barcode of Life Data System
- [MIDORI](http://www.reference-midori.info/) - Metazoan mitochondrial reference database
- [UNITE](https://unite.ut.ee/) - Fungal ITS reference database
- [Pantheon](https://www.brc.ac.uk/pantheon/) - Biological Records Centre invertebrate database
- [Darwin Tree of Life](https://portal.darwintreeoflife.org/) - Genome sequencing initiative
