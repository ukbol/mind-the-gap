#!/usr/bin/env python3
"""Publish a finished pipeline run as a new version of its Zenodo record.

Usage:
    ZENODO_TOKEN=... pipeline/publish_zenodo.py --pipeline <name> --out-dir <dir>

Reads `pipelines.<name>.results.concept_doi` from pipeline/inputs.yaml.
If the concept DOI is a placeholder or unset, prints a skip message and
exits 0 so the workflow stays green during dry-runs.

Uploads every file in <out-dir>/<pipeline>/ to a new version of the
concept DOI and publishes it. Emits the new version DOI to
GITHUB_OUTPUT as `result_doi`.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import requests
import yaml

ZENODO_API = "https://zenodo.org/api"
REPO_ROOT = Path(__file__).resolve().parent.parent
INPUTS_FILE = REPO_ROOT / "pipeline" / "inputs.yaml"
STATE_FILE = REPO_ROOT / "pipeline" / "state.json"
PLACEHOLDER_TOKEN = "PLACEHOLDER"


def _auth(token: str) -> dict:
    return {"Authorization": f"Bearer {token}"}


def _doi_to_recid(doi: str) -> str:
    return doi.rsplit("zenodo.", 1)[1].strip()


def _new_version(token: str, concept_recid: str) -> dict:
    latest = requests.get(f"{ZENODO_API}/records/{concept_recid}", timeout=60)
    latest.raise_for_status()
    latest_recid = latest.json()["id"]
    r = requests.post(
        f"{ZENODO_API}/deposit/depositions/{latest_recid}/actions/newversion",
        headers=_auth(token),
        timeout=60,
    )
    r.raise_for_status()
    draft_url = r.json()["links"]["latest_draft"]
    draft = requests.get(draft_url, headers=_auth(token), timeout=60)
    draft.raise_for_status()
    return draft.json()


def _clear_files(token: str, deposition: dict) -> None:
    for f in deposition.get("files", []):
        requests.delete(
            f"{ZENODO_API}/deposit/depositions/{deposition['id']}/files/{f['id']}",
            headers=_auth(token),
            timeout=60,
        ).raise_for_status()


def _upload(token: str, deposition: dict, path: Path) -> None:
    bucket = deposition["links"]["bucket"]
    with path.open("rb") as fh:
        r = requests.put(f"{bucket}/{path.name}", data=fh, headers=_auth(token), timeout=600)
    r.raise_for_status()


def _update_metadata(token: str, deposition: dict, pipeline_name: str, pipeline_state: dict) -> None:
    metadata_base = deposition.get("metadata", {}) or {}
    title = metadata_base.get("title", f"mind-the-gap results: {pipeline_name}")
    metadata = {
        "title": title,
        "upload_type": metadata_base.get("upload_type", "dataset"),
        "description": (
            f"Automated mind-the-gap run for pipeline <code>{pipeline_name}</code>.<br>"
            f"Input versions:<br><pre>{json.dumps(pipeline_state.get('input_versions', {}), indent=2)}</pre>"
            f"Finished: {pipeline_state.get('last_run', {}).get('finished_utc', 'unknown')}"
        ),
        "creators": metadata_base.get(
            "creators", [{"name": "mind-the-gap automation"}]
        ),
        "version": time.strftime("%Y-%m-%dT%H-%M-%SZ", time.gmtime()),
    }
    r = requests.put(
        f"{ZENODO_API}/deposit/depositions/{deposition['id']}",
        headers={**_auth(token), "Content-Type": "application/json"},
        data=json.dumps({"metadata": metadata}),
        timeout=60,
    )
    r.raise_for_status()


def _publish(token: str, deposition_id: int) -> dict:
    r = requests.post(
        f"{ZENODO_API}/deposit/depositions/{deposition_id}/actions/publish",
        headers=_auth(token),
        timeout=60,
    )
    r.raise_for_status()
    return r.json()


def _emit_gha_output(name: str, value: str) -> None:
    out = os.environ.get("GITHUB_OUTPUT")
    if not out:
        return
    with open(out, "a") as fh:
        fh.write(f"{name}={value}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pipeline", required=True, help="Pipeline name (key under `pipelines:` in inputs.yaml).")
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Root output directory; files in <out-dir>/<pipeline>/ are uploaded.",
    )
    args = parser.parse_args()

    cfg = yaml.safe_load(INPUTS_FILE.read_text())
    pipelines = cfg.get("pipelines") or {}
    if args.pipeline not in pipelines:
        raise SystemExit(f"Pipeline '{args.pipeline}' not in inputs.yaml. Known: {sorted(pipelines)}")
    pipeline_def = pipelines[args.pipeline]
    concept_doi = ((pipeline_def.get("results") or {}).get("concept_doi") or "")

    if not concept_doi or PLACEHOLDER_TOKEN in concept_doi:
        print(f"[publish] '{args.pipeline}': results.concept_doi is unset or placeholder; skipping.")
        return 0

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("[publish] ZENODO_TOKEN not set; skipping.")
        return 0

    pipeline_dir = args.out_dir / args.pipeline
    if not pipeline_dir.exists():
        print(f"[publish] no output directory at {pipeline_dir}; nothing to upload.")
        return 0

    files = sorted(p for p in pipeline_dir.iterdir() if p.is_file())
    if not files:
        print(f"[publish] no files in {pipeline_dir}; nothing to upload.")
        return 0

    state = json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}
    pipeline_state = (state.get("pipelines") or {}).get(args.pipeline, {})

    print(f"[publish] '{args.pipeline}': creating new version of {concept_doi}")
    draft = _new_version(token, _doi_to_recid(concept_doi))
    _clear_files(token, draft)
    for path in files:
        print(f"[publish]   uploading {path.name} ({path.stat().st_size} bytes)")
        _upload(token, draft, path)
    _update_metadata(token, draft, args.pipeline, pipeline_state)
    published = _publish(token, draft["id"])

    new_doi = published.get("doi") or published.get("metadata", {}).get("doi")
    print(f"[publish] '{args.pipeline}' published as {new_doi}")
    _emit_gha_output("result_doi", new_doi or "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
