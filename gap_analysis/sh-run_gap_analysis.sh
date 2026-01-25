#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=%j_gap_analysis_out.out
#SBATCH --error=%j_gap_analysis_err.err
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

module load python/3.7

# Define paths
SCRIPT_DIR="/hpc/groups/genomics-collections/ukbol_accelerated/gap_fill/mind-the-gap/gap_analysis"
SPECIES_LIST="/hpc/groups/genomics-collections/ukbol_accelerated/gap_fill/uksi_valid_species_output.tsv"
RECORDS_FILE="/hpc/groups/genomics-collections/ukbol_accelerated/gap_fill/midori_lrRNA/MIDORI2_TOTAL_NUC_GB269_lrRNA_otus_taxid.tsv"
OUTPUT_FILE="/hpc/groups/genomics-collections/ukbol_accelerated/gap_fill/midori_lrRNA/lrRNA_gap_analysis.tsv"

# Number of workers (match cpus-per-task)
WORKERS=32

cd $SCRIPT_DIR

python gap_analysis.py \
    --species-list "$SPECIES_LIST" \
    --records "$RECORDS_FILE" \
    --output "$OUTPUT_FILE" \
    --workers $WORKERS \
    --batch-size 2000

echo "Complete!"
