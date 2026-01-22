#!/bin/bash
#SBATCH --job-name=16s_gb_download
#SBATCH --output=16s_download_%j.log
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --partition=himem

OUTDIR="/mnt/shared/scratch/museomix/rbcl_genbank"
mkdir -p $OUTDIR
cd $OUTDIR

QUERY='(rbcl[Gene Name]) OR (rbc-l[Gene Name]) OR (rubisco[Gene Name])'

# Get count first
COUNT=$(esearch -db nucleotide -query "$QUERY" | xtract -pattern ENTREZ_DIRECT -element Count)
echo "Total records: $COUNT"

# Download in batches
BATCH=5000
for ((i=0; i<COUNT; i+=BATCH)); do
    OUTFILE="rbcl_batch_$(printf '%06d' $i).gb"
    
    if [[ -f "$OUTFILE" && -s "$OUTFILE" ]]; then
        echo "Skipping $OUTFILE (already exists)"
        continue
    fi
    
    # Note: -start is 1-indexed (not 0) in efetch
    START=$((i + 1))
    STOP=$((i + BATCH))
    
    echo "$(date): Fetching records $START to $STOP..."
    esearch -db nucleotide -query "$QUERY" | \
        efetch -format gb -start $START -stop $STOP \
        > "$OUTFILE"
    
    if [[ ! -s "$OUTFILE" ]]; then
        echo "WARNING: $OUTFILE is empty, retrying after 30s..."
        sleep 30
        esearch -db nucleotide -query "$QUERY" | \
            efetch -format gb -start $START -stop $STOP \
            > "$OUTFILE"
    fi
    
    sleep 3
done

echo "$(date): Download complete"
