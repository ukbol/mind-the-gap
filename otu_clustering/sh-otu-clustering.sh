#!/bin/bash
#SBATCH --partition=medium
#SBATCH --output=%j_otu_out.out
#SBATCH --error=%j_otu_err.err
#SBATCH --mem=40G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=ALL


## Conda environment
source /mnt/apps/users/bprice/conda/etc/profile.d/conda.sh
conda activate vsearch_env

# adjust threshold through t parameter

python otu_clustering.py \
    -v \
    -t 0.995 \
    --threads 32 \
    "/input/file.tsv" \
    "/output/file.tsv"

echo "Complete!"