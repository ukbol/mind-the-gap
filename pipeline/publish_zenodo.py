#!/usr/bin/env python3
"""Publish a completed mind-the-gap run as a new version of a Zenodo record.

Usage:
    ZENODO_TOKEN=... pipeline/publish_zenodo.py --out-dir final_result/run

Reads the results concept DOI from pipeline/inputs.yaml (results.concept_doi).
If unset, prints a skip message and exits 0 so the workflow doesn't fail.

Uploads every file in --out-dir to a new version of the concept DOI and
publishes it. Emits the new version DOI to GITHUB_OUTPUT as `result_doi`.
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


def _auth(token: str) -> dict:
    return {"Authorization": f"Bearer {token}"}


def _doi_to_recid(doi: str) -> str:
    return doi.rsplit("zenodo.", 1)[1].strip()


def _new_version(token: str, concept_recid: str) -> dict:
    """Create a new draft based on the latest published version."""
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


def _update_metadata(token: str, deposition: dict, state: dict) -> None:
    title = deposition.get("metadata", {}).get("title", "mind-the-gap results")
    metadata = {
        "title": title,
        "upload_type": deposition.get("metadata", {}).get("upload_type", "dataset"),
        "description": (
            "Automated mind-the-gap gap-analysis run.<br>"
            f"Input versions:<br><pre>{json.dumps(state.get('input_versions', {}), indent=2)}</pre>"
            f"Finished: {state.get('last_run', {}).get('finished_utc', 'unknown')}"
        ),
        "creators": deposition.get("metadata", {}).get(
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
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory of files to upload.")
    args = parser.parse_args()

    inputs = yaml.safe_load(INPUTS_FILE.read_text())
    concept_doi = (inputs.get("results") or {}).get("concept_doi")
    if not concept_doi:
        print("[publish] results.concept_doi not set in pipeline/inputs.yaml; skipping.")
        return 0

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("[publish] ZENODO_TOKEN not set; skipping.")
        return 0

    files = sorted(p for p in args.out_dir.iterdir() if p.is_file())
    if not files:
        print(f"[publish] no files in {args.out_dir}; nothing to upload.")
        return 0

    state = json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}

    print(f"[publish] creating new version of concept DOI {concept_doi}")
    draft = _new_version(token, _doi_to_recid(concept_doi))
    _clear_files(token, draft)
    for path in files:
        print(f"[publish]   uploading {path.name} ({path.stat().st_size} bytes)")
        _upload(token, draft, path)
    _update_metadata(token, draft, state)
    published = _publish(token, draft["id"])

    new_doi = published.get("doi") or published.get("metadata", {}).get("doi")
    print(f"[publish] published as {new_doi}")
    _emit_gha_output("result_doi", new_doi or "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
