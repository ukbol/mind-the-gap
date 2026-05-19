# Automated pipeline (GitHub Actions + Zenodo)

This document describes how `mind-the-gap` runs end-to-end on GitHub Actions,
how to update inputs, and what to configure on the GitHub side before the
workflow can actually fire.

## Architecture

```
                +----------------+
                |  Zenodo        |   <-- you push new versions of
                |  (inputs)      |       BOLD_Public.tsv, species_list.tsv, ...
                +-------+--------+
                        |
            cron / manual / repository_dispatch
                        v
       +------------------------------------+
       |   .github/workflows/run-pipeline   |
       |                                    |
       |   1. check: resolve latest DOIs    |
       |      vs pipeline/state.json        |
       |   2. run (Larger Runner, container)|
       |      fetch_zenodo.py --data-dir    |
       |      pipeline/run.sh               |
       |        bold_gene_extract           |
       |        otu_clustering              |
       |        gap_analysis                |
       |      publish_zenodo.py             |
       |   3. commit refreshed state.json   |
       +------------------------------------+
                        |
                +-------+---------+
                |  Zenodo         |   <-- new version of "results" concept DOI
                |  (results)      |       + workflow artifact in GitHub
                +-----------------+
```

## Files in this layout

| Path | Purpose |
|---|---|
| `docker/Dockerfile` | Runtime image (Python 3.11 + VSEARCH + EDirect + deps). Built and pushed to `ghcr.io/<owner>/mind-the-gap:latest` by `build-image.yml`. |
| `requirements.txt` | Python deps (BioPython, pandas, numpy, matplotlib, pyyaml, requests). |
| `pipeline/inputs.yaml` | The **only** file users need to edit to swap inputs. Concept DOIs + optional pinned version DOIs + parameters. |
| `pipeline/state.json` | Bot-managed. Records the Zenodo version DOIs consumed by the last successful run. The `check` job diffs against this to decide whether to fire. |
| `pipeline/fetch_zenodo.py` | Resolves DOIs and downloads inputs with MD5 verification. Also implements `--check-only` for the scheduled diff. |
| `pipeline/run.sh` | Orchestrator. Runs `bold_gene_extract` -> `otu_clustering` -> `gap_analysis` against the fetched inputs. |
| `pipeline/publish_zenodo.py` | Uploads `final_result/run/` to a new version of the results concept DOI. |
| `.github/workflows/run-pipeline.yml` | Main pipeline workflow (manual + scheduled + dispatch). |
| `.github/workflows/build-image.yml` | Builds and pushes the runtime container. |

## One-time setup

1. **Mint Zenodo records** for each input you want to track. Publish a v1
   even if empty, so you have a concept DOI to point at.
2. **Mint a Zenodo record for results.** Publish a v1 with a placeholder
   file. Note the **concept DOI** — that is what you put in
   `inputs.yaml -> results.concept_doi`.
3. **Edit `pipeline/inputs.yaml`** and replace each `PLACEHOLDER_*`
   concept DOI with the real ones from step 1. Leave `pinned_version`
   null unless you want to lock the workflow to a specific historical
   version.
4. **Add a Zenodo personal access token** as the GitHub Actions secret
   `ZENODO_TOKEN`. Required scope: `deposit:write` (and `deposit:actions`
   if your token UI exposes it separately).
5. **Enable a Larger Runner.** The workflow targets
   `ubuntu-latest-32-cores`; change this label in
   `.github/workflows/run-pipeline.yml` to whichever Larger Runner SKU
   your organisation has enabled. See
   <https://docs.github.com/en/actions/using-github-hosted-runners/about-larger-runners>.
6. **Allow workflow writes.** Repo Settings -> Actions -> General ->
   "Workflow permissions" -> *Read and write* (needed for the bot
   commit of `state.json`). If you'd rather not, drop the commit step
   and let `state.json` be manually refreshed.
7. **Build the container at least once** by triggering
   `build-image.yml` (it auto-fires on changes to `docker/` and
   `requirements.txt`, or run it manually via *Actions -> Build runtime
   container -> Run workflow*).

## Updating inputs

The user-facing surface for "I want to re-run with new data" is
`pipeline/inputs.yaml`. Two patterns:

- **You control the Zenodo deposit.** Just publish a new version of the
  existing concept DOI on Zenodo. The next scheduled poll (or a manual
  *Run workflow* click) will pick it up automatically. No code change
  needed.
- **You want to pin to a specific historical version.** Set
  `pinned_version` to that version DOI in `inputs.yaml`, commit, push.
  The next run will use the pinned version. Clear `pinned_version` (set
  back to `null`) to resume tracking latest.

## Triggers

- **Scheduled poll** — `cron: '0 6 * * 1'` (weekly on Mondays at 06:00
  UTC). The `check` job is the only thing that runs when there's
  nothing new; it costs almost nothing (`ubuntu-latest`, < 1 min).
- **Manual** — *Actions -> Run gap-analysis pipeline -> Run workflow*.
  Optional inputs let you pin a specific Zenodo version per-run, or
  force a re-run with `force=true` even if `state.json` is current.
- **External webhook (optional)** — anything that can POST to
  `https://api.github.com/repos/<owner>/mind-the-gap/dispatches` with
  body `{"event_type":"zenodo-new-version"}` will fire the workflow
  immediately. Zenodo itself doesn't emit webhooks today, but your
  deposit tooling can.

## Costs

GitHub-hosted Larger Runners are billed per minute. Rough numbers (as
of writing — verify current pricing):

| Runner size | Per-minute |
|---|---|
| `ubuntu-latest-4-cores` | ~$0.008 |
| `ubuntu-latest-16-cores` | ~$0.032 |
| `ubuntu-latest-32-cores` | ~$0.064 |
| `ubuntu-latest-64-cores` | ~$0.128 |

A typical full run on 32 cores (20 GB BOLD TSV download + extract +
cluster + gap analysis) is in the 20-40 min range, so ~$1-3 per run.
Set a spending limit in GitHub billing settings to cap exposure.

The `check` job runs on the free standard runner, so polling weekly
with no new versions is effectively free.

## Troubleshooting

- **`check` exits 2 with "placeholder DOIs"**: You haven't yet replaced
  the placeholder concept DOIs in `pipeline/inputs.yaml`.
- **`fetch_zenodo.py` reports a checksum mismatch**: the file on
  Zenodo changed mid-download. Re-run the workflow.
- **`publish_zenodo.py` says "skipping"**: either `results.concept_doi`
  is `null` in `inputs.yaml`, or `ZENODO_TOKEN` is not set. Decide
  which you want; if you want publishing, fix both.
- **Workflow keeps re-triggering on schedule**: confirm the bot commit
  to `state.json` is actually landing (check repo Settings -> Actions ->
  General -> Workflow permissions).

## Alternative compute backends

The recommendation here is GitHub-hosted Larger Runners because no HPC
node is available. If that changes, you can swap the `runs-on:` label
to a `[self-hosted, hpc]` runner registered on a login node and adjust
`pipeline/run.sh` to `sbatch --wait` the heavy stages. The rest of the
workflow (Zenodo fetch, state diffing, result publishing) stays the
same.

A future port to Galaxy or Nextflow is a separate deliverable and
doesn't affect this layout.
