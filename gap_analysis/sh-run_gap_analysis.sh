#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=%j_gap_analysis_out.out
#SBATCH --error=%j_gap_analysis_err.err
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=ALL

module load python/3.7

# Define paths
SCRIPT_DIR="/hpc/path"
SPECIES_LIST="/hpc/path/species_list.tsv"
RECORDS_FILE="/hpc/path/records.tsv"
OUTPUT_FILE="/hpc/path/output_results.tsv"

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
