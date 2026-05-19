#!/usr/bin/env python3
"""Resolve and download mind-the-gap inputs from Zenodo.

Reads pipeline/inputs.yaml. For each entry:
  * if `pinned_version` is set, fetches that specific version DOI;
  * otherwise resolves the latest version of the `concept_doi` via the
    Zenodo REST API.

Modes
-----
--check-only
    Print whether any input version has advanced since pipeline/state.json
    was last updated. Sets the GitHub Actions output `should_run=true` (or
    `false`). Exit code 0 always, except 2 when inputs.yaml still contains
    placeholder DOIs (so the scheduled job can no-op cleanly).

--data-dir DIR
    Download all input files into DIR. Skips files already present whose
    MD5 matches the Zenodo checksum. Writes the resolved version DOIs to
    pipeline/state.json on success.
"""

from __future__ import annotations

import argparse
import hashlib
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


def _doi_to_recid(doi: str) -> str:
    """Extract the Zenodo record id from a DOI like '10.5281/zenodo.123456'."""
    if "zenodo." not in doi:
        raise ValueError(f"Not a Zenodo DOI: {doi}")
    return doi.rsplit("zenodo.", 1)[1].strip()


def _get(url: str, **kwargs) -> requests.Response:
    """GET with simple retry on transient errors."""
    last_exc: Exception | None = None
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=60, **kwargs)
            if r.status_code >= 500:
                raise requests.HTTPError(f"{r.status_code} from {url}")
            r.raise_for_status()
            return r
        except (requests.RequestException, requests.HTTPError) as exc:
            last_exc = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Zenodo request failed after retries: {url}") from last_exc


def resolve_latest_version(concept_doi: str) -> dict:
    """Return the Zenodo record JSON for the latest version of a concept DOI."""
    concept_recid = _doi_to_recid(concept_doi)
    # The concept record's `links.latest` redirects to the latest-version record.
    record = _get(f"{ZENODO_API}/records/{concept_recid}").json()
    latest_url = record.get("links", {}).get("latest")
    if not latest_url:
        return record
    return _get(latest_url).json()


def resolve_specific_version(version_doi: str) -> dict:
    return _get(f"{ZENODO_API}/records/{_doi_to_recid(version_doi)}").json()


def _md5(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 - Zenodo metadata uses MD5
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_file(file_entry: dict, dest: Path) -> None:
    url = file_entry["links"]["self"]
    checksum = file_entry.get("checksum", "")  # format: "md5:<hex>"
    expected_md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else None

    if dest.exists() and expected_md5 and _md5(dest) == expected_md5:
        print(f"[fetch] cached: {dest.name}", flush=True)
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[fetch] downloading {dest.name} <- {url}", flush=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
    if expected_md5:
        actual = _md5(tmp)
        if actual != expected_md5:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"Checksum mismatch for {dest.name}: expected {expected_md5}, got {actual}"
            )
    tmp.rename(dest)


def _load_inputs() -> dict:
    with INPUTS_FILE.open() as fh:
        return yaml.safe_load(fh)


def _load_state() -> dict:
    if not STATE_FILE.exists():
        return {"schema_version": 1, "last_run": None, "input_versions": {}}
    with STATE_FILE.open() as fh:
        return json.load(fh)


def _has_placeholder(inputs: dict) -> bool:
    for entry in inputs.get("inputs", {}).values():
        if PLACEHOLDER_TOKEN in (entry.get("concept_doi") or ""):
            return True
    return False


def _resolve_all(inputs: dict) -> dict[str, dict]:
    resolved: dict[str, dict] = {}
    for name, entry in inputs["inputs"].items():
        if entry.get("pinned_version"):
            record = resolve_specific_version(entry["pinned_version"])
        else:
            record = resolve_latest_version(entry["concept_doi"])
        resolved[name] = {
            "version_doi": record.get("doi") or record.get("conceptdoi"),
            "record": record,
            "filename": entry["filename"],
        }
    return resolved


def _emit_gha_output(name: str, value: str) -> None:
    out = os.environ.get("GITHUB_OUTPUT")
    if not out:
        return
    with open(out, "a") as fh:
        fh.write(f"{name}={value}\n")


def cmd_check_only() -> int:
    inputs = _load_inputs()
    if _has_placeholder(inputs):
        print("[check] inputs.yaml still contains placeholder DOIs; refusing to run.")
        _emit_gha_output("should_run", "false")
        return 2

    state = _load_state()
    resolved = _resolve_all(inputs)

    advanced: list[str] = []
    for name, info in resolved.items():
        prev = state["input_versions"].get(name)
        current = info["version_doi"]
        if prev != current:
            advanced.append(f"{name}: {prev} -> {current}")

    if advanced:
        print("[check] new input versions detected:")
        for line in advanced:
            print(f"  - {line}")
        _emit_gha_output("should_run", "true")
    else:
        print("[check] all inputs already at the last-run versions; no run needed.")
        _emit_gha_output("should_run", "false")
    return 0


def cmd_download(data_dir: Path) -> int:
    inputs = _load_inputs()
    if _has_placeholder(inputs):
        print("[fetch] inputs.yaml still contains placeholder DOIs; aborting.")
        return 2

    resolved = _resolve_all(inputs)
    data_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, str] = {}
    for name, info in resolved.items():
        target_name = info["filename"]
        files = info["record"].get("files", [])
        match = next((f for f in files if f.get("key") == target_name), None)
        if match is None:
            raise RuntimeError(
                f"File '{target_name}' not found in Zenodo record for input '{name}'. "
                f"Available: {[f.get('key') for f in files]}"
            )
        dest = data_dir / target_name
        _download_file(match, dest)
        paths[name] = str(dest)

    # Refresh state.json with the versions we just downloaded.
    state = _load_state()
    state["input_versions"] = {name: info["version_doi"] for name, info in resolved.items()}
    state["last_run"] = {
        "completed": False,  # run.sh will flip this to True on success
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with STATE_FILE.open("w") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
        fh.write("\n")

    # Emit resolved paths for downstream steps to consume.
    manifest = data_dir / "manifest.json"
    with manifest.open("w") as fh:
        json.dump({"paths": paths, "versions": state["input_versions"]}, fh, indent=2)
    print(f"[fetch] wrote {manifest}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--check-only", action="store_true", help="Only check for new versions; set GHA output should_run.")
    g.add_argument("--data-dir", type=Path, help="Download all inputs into this directory.")
    args = parser.parse_args()

    if args.check_only:
        return cmd_check_only()
    return cmd_download(args.data_dir)


if __name__ == "__main__":
    sys.exit(main())
