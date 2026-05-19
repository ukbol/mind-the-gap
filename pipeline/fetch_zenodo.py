#!/usr/bin/env python3
"""Resolve and download mind-the-gap inputs from Zenodo (multi-pipeline).

Reads pipeline/inputs.yaml (schema v2: a map of named pipelines).

Modes
-----
--list-pipelines-needing-run [--only A,B,...] [--force]
    Walk every (or the --only-restricted set of) pipeline. For each,
    resolve its inputs' current latest version DOIs and compare with
    pipeline/state.json. Apply dependency closure: if a producer
    pipeline (referenced via depends_on) is in the to-run set, every
    consumer of that producer is added. Emit the resulting list as
    JSON in topological order. With --force, every selected pipeline
    is unconditionally emitted.
    Writes to GITHUB_OUTPUT as `pipelines_to_run=<json>`.

--pipeline <name> --data-dir <dir>
    Download the named pipeline's raw inputs into <dir>/<name>/.
    Decompresses .gz / .tar.gz inputs after download. Resolves
    `depends_on` entries by pulling the upstream producer's latest
    published `results.primary_file` from Zenodo. Writes a per-pipeline
    manifest.json mapping logical input names to on-disk paths and the
    version DOIs consumed.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
from collections import deque
from pathlib import Path

import requests
import yaml

ZENODO_API = "https://zenodo.org/api"
REPO_ROOT = Path(__file__).resolve().parent.parent
INPUTS_FILE = REPO_ROOT / "pipeline" / "inputs.yaml"
STATE_FILE = REPO_ROOT / "pipeline" / "state.json"
PLACEHOLDER_TOKEN = "PLACEHOLDER"


# ----------------------------- Zenodo HTTP -----------------------------------


def _get(url: str, **kwargs) -> requests.Response:
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


def _doi_to_recid(doi: str) -> str:
    if "zenodo." not in doi:
        raise ValueError(f"Not a Zenodo DOI: {doi}")
    return doi.rsplit("zenodo.", 1)[1].strip()


def resolve_latest_version(concept_doi: str) -> dict:
    record = _get(f"{ZENODO_API}/records/{_doi_to_recid(concept_doi)}").json()
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
    checksum = file_entry.get("checksum", "")
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


# ----------------------------- Config loading --------------------------------


def _load_inputs() -> dict:
    return yaml.safe_load(INPUTS_FILE.read_text())


def _load_state() -> dict:
    if not STATE_FILE.exists():
        return {"schema_version": 2, "pipelines": {}}
    state = json.loads(STATE_FILE.read_text())
    state.setdefault("schema_version", 2)
    state.setdefault("pipelines", {})
    return state


def _pipeline_has_placeholder(pipeline_def: dict) -> bool:
    for entry in (pipeline_def.get("inputs") or {}).values():
        if PLACEHOLDER_TOKEN in (entry.get("concept_doi") or ""):
            return True
    return False


def _all_pipeline_names(cfg: dict) -> list[str]:
    return list((cfg.get("pipelines") or {}).keys())


# ----------------------------- Diff / DAG ------------------------------------


def _resolve_input_versions(pipeline_def: dict) -> dict[str, str]:
    """Return {logical_input_name: current_version_doi} for raw inputs."""
    versions: dict[str, str] = {}
    for name, entry in (pipeline_def.get("inputs") or {}).items():
        if entry.get("pinned_version"):
            record = resolve_specific_version(entry["pinned_version"])
        else:
            record = resolve_latest_version(entry["concept_doi"])
        versions[name] = record.get("doi") or record.get("conceptdoi")
    return versions


def _resolve_dependency_versions(
    pipeline_def: dict, cfg: dict
) -> dict[str, str]:
    """Return {dep_logical_name: producer_results_version_doi}.

    Each depends_on entry points at a producer pipeline by name. The
    'version' a consumer cares about is the producer's published
    results record (results.concept_doi resolved to its latest
    version). If the producer hasn't published yet, returns '' which
    means 'unknown' and will trigger a re-run.
    """
    out: dict[str, str] = {}
    for dep_name, producer_name in (pipeline_def.get("depends_on") or {}).items():
        producer = (cfg.get("pipelines") or {}).get(producer_name)
        if producer is None:
            raise KeyError(
                f"depends_on points at unknown pipeline '{producer_name}'."
            )
        results = producer.get("results") or {}
        results_doi = results.get("concept_doi", "")
        if not results_doi or PLACEHOLDER_TOKEN in results_doi:
            out[dep_name] = ""
            continue
        record = resolve_latest_version(results_doi)
        out[dep_name] = record.get("doi") or record.get("conceptdoi") or ""
    return out


def _pipeline_needs_run(
    pipeline_name: str,
    pipeline_def: dict,
    cfg: dict,
    state: dict,
) -> tuple[bool, dict[str, str]]:
    """Return (needs_run, resolved_versions) for one pipeline.

    `resolved_versions` covers both raw inputs and depends_on producers.
    """
    if _pipeline_has_placeholder(pipeline_def):
        # Placeholder concept DOIs mean the pipeline isn't configured yet.
        return False, {}

    resolved = _resolve_input_versions(pipeline_def)
    resolved.update(_resolve_dependency_versions(pipeline_def, cfg))

    prev = state.get("pipelines", {}).get(pipeline_name, {})
    prev_versions = prev.get("input_versions") or {}

    if not prev_versions:
        return True, resolved
    for k, v in resolved.items():
        if prev_versions.get(k) != v:
            return True, resolved
    # Extra safety: if a key disappeared from current resolution, still consider it changed.
    for k in prev_versions:
        if k not in resolved:
            return True, resolved
    return False, resolved


def _topo_sort(cfg: dict, pipeline_names: list[str]) -> list[str]:
    """Topological sort of a subset of pipelines by depends_on edges.

    Producers come before consumers. Pipelines outside `pipeline_names`
    are dropped (they're not being run).
    """
    in_set = set(pipeline_names)
    pipelines = cfg.get("pipelines") or {}
    # Build adjacency limited to in_set.
    indeg = {name: 0 for name in in_set}
    children: dict[str, list[str]] = {name: [] for name in in_set}
    for name in in_set:
        for dep_producer in (pipelines[name].get("depends_on") or {}).values():
            if dep_producer in in_set:
                children[dep_producer].append(name)
                indeg[name] += 1

    queue = deque(sorted(n for n, d in indeg.items() if d == 0))
    out: list[str] = []
    while queue:
        n = queue.popleft()
        out.append(n)
        for c in sorted(children[n]):
            indeg[c] -= 1
            if indeg[c] == 0:
                queue.append(c)

    if len(out) != len(in_set):
        cycle = sorted(in_set - set(out))
        raise RuntimeError(f"Cycle in pipeline depends_on graph involving: {cycle}")
    return out


def _apply_dependency_closure(cfg: dict, to_run: set[str]) -> set[str]:
    """Add every consumer of an in-set producer to the in-set."""
    pipelines = cfg.get("pipelines") or {}
    closed = set(to_run)
    changed = True
    while changed:
        changed = False
        for name, defn in pipelines.items():
            if name in closed:
                continue
            for producer in (defn.get("depends_on") or {}).values():
                if producer in closed:
                    closed.add(name)
                    changed = True
                    break
    return closed


# ----------------------------- Decompression ---------------------------------


def _is_gzipped(entry: dict, filename: str) -> bool:
    explicit = entry.get("gzipped")
    if isinstance(explicit, bool):
        return explicit
    return filename.endswith(".gz") or filename.endswith(".tar.gz") or filename.endswith(".tgz")


def _decompressed_path(compressed: Path) -> Path:
    """For .tsv.gz -> .tsv. For .tar.gz -> sibling directory."""
    name = compressed.name
    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        stem = name[:-7] if name.endswith(".tar.gz") else name[:-4]
        return compressed.parent / stem
    if name.endswith(".gz"):
        return compressed.parent / name[:-3]
    return compressed


def _decompress(compressed: Path) -> Path:
    """Decompress `compressed` in-place; return the path to the decompressed
    artifact (a file for single-file gzip, a directory for tarballs).
    Idempotent: if the decompressed artifact already exists and looks
    fresher than the compressed file, return it unchanged.
    """
    name = compressed.name
    out = _decompressed_path(compressed)

    if out.exists() and out.stat().st_mtime >= compressed.stat().st_mtime:
        print(f"[fetch] decompressed cache hit: {out.name}", flush=True)
        return out

    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        print(f"[fetch] extracting tarball: {name}", flush=True)
        out.mkdir(parents=True, exist_ok=True)
        with tarfile.open(compressed, "r:gz") as tar:
            # `filter="data"` is the safe default in Python 3.12+, omitted on 3.11.
            try:
                tar.extractall(out, filter="data")
            except TypeError:
                tar.extractall(out)
    elif name.endswith(".gz"):
        print(f"[fetch] gunzipping: {name}", flush=True)
        tmp = out.with_suffix(out.suffix + ".part")
        with gzip.open(compressed, "rb") as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)
        tmp.rename(out)
    else:
        # Shouldn't happen; _is_gzipped would have returned False.
        return compressed
    return out


# ----------------------------- Commands --------------------------------------


def _emit_gha_output(name: str, value: str) -> None:
    out = os.environ.get("GITHUB_OUTPUT")
    if not out:
        return
    with open(out, "a") as fh:
        fh.write(f"{name}={value}\n")


def cmd_list(only: list[str] | None, force: bool) -> int:
    cfg = _load_inputs()
    state = _load_state()
    all_names = _all_pipeline_names(cfg)
    if not all_names:
        print("[check] inputs.yaml defines no pipelines.")
        _emit_gha_output("pipelines_to_run", "[]")
        return 0

    if only:
        unknown = [n for n in only if n not in all_names]
        if unknown:
            raise SystemExit(f"Unknown pipeline(s) in --only: {unknown}. Known: {all_names}")
        candidates = only
    else:
        candidates = all_names

    needs_run: set[str] = set()
    placeholder_skipped: list[str] = []
    for name in candidates:
        defn = cfg["pipelines"][name]
        if _pipeline_has_placeholder(defn):
            placeholder_skipped.append(name)
            continue
        if force:
            needs_run.add(name)
            continue
        try:
            run, _ = _pipeline_needs_run(name, defn, cfg, state)
        except Exception as exc:
            print(f"[check] {name}: error resolving versions: {exc}")
            run = False
        if run:
            needs_run.add(name)
            print(f"[check] {name}: NEW versions detected")
        else:
            print(f"[check] {name}: up to date")

    if placeholder_skipped:
        print(
            "[check] skipping pipelines with placeholder DOIs: "
            + ", ".join(placeholder_skipped)
        )

    # Dependency closure: any consumer of an in-set producer must also run.
    needs_run = _apply_dependency_closure(cfg, needs_run)
    # Filter out placeholder pipelines from the closure.
    needs_run -= set(placeholder_skipped)

    # Order topologically (producer-first).
    ordered = _topo_sort(cfg, sorted(needs_run))

    print(f"[check] pipelines_to_run = {ordered}")
    _emit_gha_output("pipelines_to_run", json.dumps(ordered))
    return 0


def cmd_fetch(pipeline_name: str, data_dir: Path) -> int:
    cfg = _load_inputs()
    pipelines = cfg.get("pipelines") or {}
    if pipeline_name not in pipelines:
        raise SystemExit(f"Pipeline '{pipeline_name}' not in inputs.yaml. Known: {sorted(pipelines)}")
    pipeline_def = pipelines[pipeline_name]

    if _pipeline_has_placeholder(pipeline_def):
        print(f"[fetch] pipeline '{pipeline_name}' has placeholder DOIs; aborting.")
        return 2

    target_dir = data_dir / pipeline_name
    target_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, str] = {}
    versions: dict[str, str] = {}

    # Raw inputs.
    for name, entry in (pipeline_def.get("inputs") or {}).items():
        if entry.get("pinned_version"):
            record = resolve_specific_version(entry["pinned_version"])
        else:
            record = resolve_latest_version(entry["concept_doi"])
        versions[name] = record.get("doi") or record.get("conceptdoi") or ""

        target_name = entry["filename"]
        match = next((f for f in record.get("files", []) if f.get("key") == target_name), None)
        if match is None:
            raise RuntimeError(
                f"File '{target_name}' not found in Zenodo record for input '{name}'. "
                f"Available: {[f.get('key') for f in record.get('files', [])]}"
            )
        dest = target_dir / target_name
        _download_file(match, dest)

        if _is_gzipped(entry, target_name):
            decompressed = _decompress(dest)
            paths[name] = str(decompressed)
        else:
            paths[name] = str(dest)

    # Dependencies: fetch the producer's latest published primary_file.
    for dep_name, producer_name in (pipeline_def.get("depends_on") or {}).items():
        producer = pipelines.get(producer_name)
        if producer is None:
            raise KeyError(f"depends_on points at unknown pipeline '{producer_name}'.")
        results = producer.get("results") or {}
        results_doi = results.get("concept_doi", "")
        primary_file = results.get("primary_file")
        if not results_doi or PLACEHOLDER_TOKEN in results_doi or not primary_file:
            raise RuntimeError(
                f"Pipeline '{pipeline_name}' depends on '{producer_name}' but that "
                f"producer has no published results record (concept_doi={results_doi!r}, "
                f"primary_file={primary_file!r}). Run the producer first."
            )
        record = resolve_latest_version(results_doi)
        versions[dep_name] = record.get("doi") or record.get("conceptdoi") or ""
        match = next((f for f in record.get("files", []) if f.get("key") == primary_file), None)
        if match is None:
            raise RuntimeError(
                f"primary_file '{primary_file}' missing from producer '{producer_name}' record. "
                f"Available: {[f.get('key') for f in record.get('files', [])]}"
            )
        dest = target_dir / primary_file
        _download_file(match, dest)
        paths[dep_name] = str(dest)

    # State patch for this pipeline (written by run.py post-run; we just
    # ship a manifest the orchestrator and publisher will consume).
    manifest = {
        "pipeline": pipeline_name,
        "paths": paths,
        "versions": versions,
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (target_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"[fetch] wrote {target_dir / 'manifest.json'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument(
        "--list-pipelines-needing-run",
        action="store_true",
        help="Emit pipelines_to_run JSON (topologically ordered) to GITHUB_OUTPUT.",
    )
    g.add_argument(
        "--pipeline",
        type=str,
        help="Download inputs for the named pipeline.",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated subset of pipelines (used with --list-pipelines-needing-run).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Treat every selected pipeline as needing a run.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Output directory for --pipeline downloads (defaults to ./data).",
    )
    args = parser.parse_args()

    if args.list_pipelines_needing_run:
        only = [s.strip() for s in args.only.split(",")] if args.only else None
        return cmd_list(only, args.force)

    data_dir = args.data_dir or (REPO_ROOT / "data")
    return cmd_fetch(args.pipeline, data_dir)


if __name__ == "__main__":
    sys.exit(main())
