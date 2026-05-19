#!/usr/bin/env bash
# Pipeline driver for the automated mind-the-gap run.
#
# Reads:
#   $DATA_DIR/manifest.json  - produced by pipeline/fetch_zenodo.py
#   pipeline/inputs.yaml     - parameters block
#
# Writes:
#   $OUT_DIR/                - all intermediate + final TSVs
#
# Stages: bold_gene_extract -> otu_clustering -> gap_analysis
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"
OUT_DIR="${OUT_DIR:-$REPO_ROOT/final_result/run}"
mkdir -p "$OUT_DIR"

read_yaml() {
    # $1: dotted path into inputs.yaml, e.g. parameters.marker
    python - "$1" <<'PY'
import sys, yaml, pathlib
key = sys.argv[1].split(".")
data = yaml.safe_load(pathlib.Path("pipeline/inputs.yaml").read_text())
for k in key:
    data = data[k]
print(data)
PY
}

read_path() {
    # $1: input name in manifest.json (e.g. bold, species_list)
    python - "$1" <<'PY'
import json, sys, os, pathlib
name = sys.argv[1]
manifest = json.loads(pathlib.Path(os.environ["DATA_DIR"], "manifest.json").read_text())
print(manifest["paths"][name])
PY
}

export DATA_DIR

cd "$REPO_ROOT"

MARKER="$(read_yaml parameters.marker)"
OTU_THRESHOLD="$(read_yaml parameters.otu_threshold)"
WORKERS="$(read_yaml parameters.workers)"
BATCH_SIZE="$(read_yaml parameters.batch_size)"

BOLD_TSV="$(read_path bold)"
SPECIES_LIST="$(read_path species_list)"

FILTERED="$OUT_DIR/bold_filtered.tsv"
CLUSTERED="$OUT_DIR/bold_clustered.tsv"
GAP_OUT="$OUT_DIR/gap_analysis.tsv"

echo "[run] Stage 1/3: bold_gene_extract  (marker=$MARKER)"
python bold_processing/bold_gene_extract/bold_gene_extract.py \
    -g "$MARKER" \
    "$BOLD_TSV" "$FILTERED"

echo "[run] Stage 2/3: otu_clustering    (threshold=$OTU_THRESHOLD, threads=$WORKERS)"
python otu_clustering/otu_clustering.py \
    -t "$OTU_THRESHOLD" \
    --threads "$WORKERS" \
    "$FILTERED" "$CLUSTERED"

echo "[run] Stage 3/3: gap_analysis      (workers=$WORKERS, batch_size=$BATCH_SIZE)"
python gap_analysis/gap_analysis.py \
    --species-list "$SPECIES_LIST" \
    --records "$CLUSTERED" \
    --output "$GAP_OUT" \
    --workers "$WORKERS" \
    --batch-size "$BATCH_SIZE"

# Mark state.json with completion so the bot commit reflects a successful run.
python - <<'PY'
import json, pathlib, time
p = pathlib.Path("pipeline/state.json")
state = json.loads(p.read_text())
state.setdefault("last_run", {})
state["last_run"]["completed"] = True
state["last_run"]["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
p.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
PY

echo "[run] complete. Results in $OUT_DIR"
ls -lh "$OUT_DIR"
