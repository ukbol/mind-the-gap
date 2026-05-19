# Automated pipelines (GitHub Actions + Zenodo)

This document describes how the gap-analysis suite runs end-to-end on
GitHub Actions, how to add or update inputs, and what to configure on
the GitHub side before the workflow can fire.

## Pipelines

The workflow drives nine pipelines, in two layers:

```
+----------------------------------------------------------+
| Layer 1: producer                                        |
|                                                          |
|   uksi_species_list  (UKSI/Pantheon/JNCC/freshwater)     |
|     uksi_import -> uksi_export -> uksi_species_export.tsv|
+--------------------+-------------------------------------+
                     |
                     v   (every consumer below depends_on this)
+----------------------------------------------------------+
| Layer 2: gap-analysis pipelines                          |
|                                                          |
|   bold_coi           BOLD_Public.tsv.gz  -> bold_gene_extract -> gap_analysis --bold
|   bold_rbcl          BOLD_Public.tsv.gz  -> bold_gene_extract -> otu_clustering -> gap_analysis
|   midori_12s         midori_12s.fasta.gz -> process_midori    -> otu_clustering -> gap_analysis
|   midori_16s         midori_16s.fasta.gz -> process_midori    -> otu_clustering -> gap_analysis
|   unite_its          unite.fasta.gz      -> process_unite     -> gap_analysis
|   ncbi_rbcl          ncbi_rbcl.gb.tar.gz -> ncbi_gb_extract   -> otu_clustering -> gap_analysis
|   dtol_genomes       dtol.csv            -> dtol_status
|   ena_mitogenomes    mito_species_ena.tsv -> mito_status
+----------------------------------------------------------+
```

Each pipeline:
- has its own raw Zenodo input(s) tracked by **concept DOI**;
- writes a deliverable that's published as a **new version** of its
  own results record;
- runs only when one of its inputs (or, for consumers, its producer's
  published results) has advanced since the last successful run.

`bold_coi` skips OTU clustering because BOLD already provides BIN
URIs. `unite_its` skips clustering because UNITE's source FASTA
already carries OTU IDs.

## Files in this layout

| Path | Purpose |
|---|---|
| `pipeline/inputs.yaml` | The **only** file you edit to add or update a pipeline. Maps logical pipeline names to inputs, stages, and a results record. |
| `pipeline/state.json` | Bot-managed. Per-pipeline `input_versions` + `last_run` timestamps. The `check` job diffs against this to decide which matrix legs to run. |
| `pipeline/stages.py` | Stage dispatch table. Maps `name:` in inputs.yaml to a script + invocation convention. Add new stages here. |
| `pipeline/fetch_zenodo.py` | Resolves concept DOIs, downloads inputs, decompresses `.gz`/`.tar.gz`, resolves `depends_on` producer outputs. `--list-pipelines-needing-run` emits the matrix list (topologically ordered, with dependency closure). |
| `pipeline/run.py` | Per-pipeline orchestrator. `python pipeline/run.py --pipeline <name>` reads the manifest produced by `fetch_zenodo.py`, walks the stage list, threads outputs between stages, and updates `state.json`. |
| `pipeline/publish_zenodo.py` | `--pipeline <name>` publishes that pipeline's `final_result/<name>/` to a new version of its results concept DOI. |
| `.github/workflows/run-pipeline.yml` | Matrix workflow: `check` job emits `pipelines_to_run` JSON; `run` job runs each pipeline in a Larger Runner container, serially. |
| `.github/workflows/build-image.yml` | Builds and pushes the runtime container. Unchanged from the single-pipeline version. |
| `docker/Dockerfile`, `requirements.txt` | Runtime image (Python 3.11 + VSEARCH + EDirect + deps). |
| `uksi_processing/uksi_db/uksi_import.py`, `uksi_export.py` | UKSI species-list builder. Hardcoded paths now have env-var overrides (`UKSI_NAMES_TSV`, `UKSI_DB`, `UKSI_VALID_OUT`, etc.) so the pipeline can drive them. Interactive use with the historical defaults still works. |

## One-time setup

1. **Mint Zenodo records** for every raw input you want to track. Each
   needs at least a v1 publish so a concept DOI exists. Then mint a
   **results** record per pipeline (eight of them — `uksi_species_list`
   plus the seven gap-analysis pipelines). Publish a v1 with a
   placeholder file in each.
2. **Edit `pipeline/inputs.yaml`** and replace every `PLACEHOLDER_*`
   concept DOI with the real one from step 1. Leave `pinned_version`
   absent unless you want to lock to a specific historical version.
3. **Add a Zenodo personal access token** as the GitHub Actions secret
   `ZENODO_TOKEN`. Required scopes: `deposit:write`.
4. **Enable a Larger Runner.** The workflow targets
   `ubuntu-latest-32-cores`; change this label in
   `.github/workflows/run-pipeline.yml` if your org uses a different
   SKU. See
   <https://docs.github.com/en/actions/using-github-hosted-runners/about-larger-runners>.
5. **Allow workflow writes.** Repo Settings -> Actions -> General ->
   Workflow permissions -> *Read and write* (for the bot commit of
   `state.json`).
6. **Build the container at least once** by triggering
   `build-image.yml`. It auto-fires on `docker/` and
   `requirements.txt` changes.

## Updating inputs

The user-facing surface for "I want to re-run with new data" is
`pipeline/inputs.yaml`. Two common patterns:

- **You control the Zenodo deposit for a raw input.** Publish a new
  version of the existing concept DOI on Zenodo. The next scheduled
  poll (or a manual *Run workflow* click) will pick it up. Every
  pipeline that consumes that concept DOI will re-run; pipelines
  that don't are skipped.
- **You want to pin a single pipeline to a specific historical
  version.** Set `pinned_version` on the relevant input row in
  `inputs.yaml`. Clear it (delete the line) to resume tracking
  latest.

## Adding a new pipeline

Edit `pipeline/inputs.yaml`. Most new pipelines are a YAML-only edit:

```yaml
pipelines:
  bold_plants_matk:
    inputs:
      raw:
        concept_doi: 10.5281/zenodo.PLACEHOLDER_BOLD_RAW
        filename: BOLD_Public.tsv.gz
        gzipped: true
    depends_on:
      species_list: uksi_species_list
    stages:
      - name: bold_gene_extract
        args: { gene: matK }
      - name: otu_clustering
        args: { threshold: 0.97, threads: 32 }
      - name: gap_analysis
        args:
          workers: 32
          batch_size: 2000
          cluster_column: otu_id
    results:
      concept_doi: 10.5281/zenodo.PLACEHOLDER_RESULTS_BOLD_MATK
      primary_file: 03_gap_analysis.tsv
```

A truly new stage (a new analysis script) needs a new entry in
`pipeline/stages.py` mapping the stage name to its script and
invocation convention. Three conventions exist:

- `positional_inout` — `script INPUT OUTPUT [flags...]`.
- `flagged` — `script --primary-flag IN --output OUT [aux flags]`.
- `script_main` — for scripts driven by env vars / no CLI today
  (currently only the UKSI import/export pair).

## Triggers

- **Scheduled poll** — daily at 06:00 UTC. The `check` job runs on the
  free standard runner; if no pipeline's inputs advanced, the matrix
  doesn't fan out and the run costs nothing.
- **Manual** — *Actions -> Run gap-analysis pipelines -> Run workflow*.
  - `pipelines` input: leave as `all` to consider every pipeline, or
    pass a comma-separated subset like `bold_coi,midori_12s`.
  - `force`: re-run the selected pipelines even if their inputs
    haven't advanced.
- **External webhook** — anything that can POST to
  `https://api.github.com/repos/<owner>/mind-the-gap/dispatches`
  with body `{"event_type":"zenodo-new-version"}` will fire the
  workflow immediately.

## Dependency closure and ordering

If `uksi_species_list` has new inputs, every consumer of its results
is added to the matrix automatically (closure). The matrix runs
strictly serially (`max-parallel: 1`) in topological order, so
`uksi_species_list` finishes and publishes its results record before
any consumer's `fetch_zenodo.py` resolves the producer's latest
version.

The producer-then-consumer dance happens via Zenodo: consumers
dereference the producer's `results.concept_doi` -> latest version ->
`primary_file`. They don't share workspace with the producer.

## Costs

GitHub-hosted Larger Runners are billed per minute. Approximate per-run
wall-clock per pipeline:

| Pipeline | Heavy stage | Typical runtime |
|---|---|---|
| `bold_coi` | bold_gene_extract + gap_analysis (BOLD mode on ~20 GB) | 20-40 min |
| `bold_rbcl` | same + otu_clustering | 30-60 min |
| `midori_12s` / `midori_16s` | process_midori + otu_clustering | ~10-20 min each |
| `unite_its` | process_unite + gap_analysis | ~10 min |
| `ncbi_rbcl` | ncbi_gb_extract + otu_clustering | ~10-20 min |
| `uksi_species_list` | uksi_import + uksi_export | <5 min |
| `dtol_genomes` / `ena_mitogenomes` | single annotator step | <5 min each |

Add a safety margin for Zenodo download time on the large compressed
inputs (BOLD ~4-5 GB gzipped; MIDORI/UNITE/NCBI smaller).

With `max-parallel: 1`, the bill is the sum of per-pipeline times — at
~$0.064/min for `ubuntu-latest-32-cores` a "everything reran" day is
in the ~$5-8 range. Set a spending limit in GitHub billing to cap
exposure.

The `check` job runs on the free standard runner, so polling daily
with no new versions is effectively free.

## Troubleshooting

- **`check` emits an empty `pipelines_to_run`**: either nothing changed,
  or every candidate pipeline still has placeholder DOIs. Look at the
  `Resolve pipelines needing a run` step's log: each pipeline is logged
  as either `up to date`, `NEW versions detected`, or
  `skipping (placeholder DOIs)`.
- **`fetch_zenodo.py` reports a checksum mismatch**: the file on
  Zenodo changed mid-download. Re-run the workflow.
- **A consumer pipeline fails with "producer has no published results
  record"**: the producer's `results.concept_doi` still contains
  `PLACEHOLDER_*`. Publish a v1 of the producer's results record and
  update inputs.yaml.
- **`publish_zenodo.py` says "skipping"**: either the pipeline's
  `results.concept_doi` is a placeholder, or `ZENODO_TOKEN` is not
  set. Decide which you want; if you want publishing, fix both.
- **Workflow keeps re-triggering on schedule**: confirm the bot commit
  to `state.json` is actually landing (Repo Settings -> Actions ->
  General -> Workflow permissions -> Read and write).

## Alternative compute backends

This layout targets GitHub-hosted Larger Runners because no HPC node
is available today. If that changes, you can swap the `runs-on:` label
to a `[self-hosted, hpc]` runner registered on a login node and adjust
`pipeline/run.py` to `sbatch --wait` the heavy stages. The rest of the
workflow (Zenodo fetch, state diffing, result publishing) stays the
same.

A future port to Galaxy or Nextflow is a separate deliverable and
doesn't affect this layout.
