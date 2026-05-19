#!/usr/bin/env python3
"""Per-pipeline orchestrator for the multi-pipeline mind-the-gap suite.

Usage:
    DATA_DIR=... OUT_DIR=... python pipeline/run.py --pipeline <name>

Reads:
    pipeline/inputs.yaml                - pipeline definition (stages, args).
    $DATA_DIR/<pipeline>/manifest.json  - input paths produced by fetch_zenodo.py.

Writes:
    $OUT_DIR/<pipeline>/<NN>_<stage>.tsv  - per-stage outputs.
    pipeline/state.json                   - records the version DOIs and
                                            finish timestamp for this pipeline.

The chosen pipeline's stage list is executed in order. Each stage's
output is wired as the primary input to the next stage. Auxiliary
inputs (the species list, DToL/mito metadata) are resolved from the
manifest.

This script does NOT push to Zenodo; that's publish_zenodo.py.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUTS_FILE = REPO_ROOT / "pipeline" / "inputs.yaml"
STATE_FILE = REPO_ROOT / "pipeline" / "state.json"

sys.path.insert(0, str(REPO_ROOT / "pipeline"))
import stages as stages_module  # noqa: E402


@dataclass
class RunContext:
    pipeline_name: str
    pipeline_def: dict
    data_dir: Path           # data root (per-pipeline subdir lives at data_dir / pipeline_name)
    out_dir: Path            # outputs root (per-pipeline subdir lives at out_dir / pipeline_name)
    manifest: dict           # {"paths": {...}, "versions": {...}}
    env_state: dict[str, str] = field(default_factory=dict)  # env vars produced by script_main stages
    prev_output: Path | None = None

    @property
    def pipeline_data_dir(self) -> Path:
        return self.data_dir / self.pipeline_name

    @property
    def pipeline_out_dir(self) -> Path:
        return self.out_dir / self.pipeline_name


def _load_inputs() -> dict:
    return yaml.safe_load(INPUTS_FILE.read_text())


def _load_manifest(ctx_dir: Path) -> dict:
    path = ctx_dir / "manifest.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Expected manifest.json at {path}. "
            "Run pipeline/fetch_zenodo.py --pipeline <name> --data-dir <dir> first."
        )
    return json.loads(path.read_text())


def _input_path(ctx: RunContext, logical_name: str) -> Path:
    paths = ctx.manifest.get("paths", {})
    if logical_name not in paths:
        raise KeyError(
            f"Pipeline '{ctx.pipeline_name}' has no input or dependency named "
            f"'{logical_name}' in manifest. Known keys: {sorted(paths)}"
        )
    return Path(paths[logical_name])


def _pipeline_primary_input(ctx: RunContext) -> Path:
    """The pipeline-level input that flows into the first 'chain' stage.

    Convention: pipelines have exactly one such input, conventionally
    keyed `raw`, but some terminal pipelines key it differently
    (e.g. `dtol_metadata`, `mito_metadata`). We pick whichever input
    isn't a depends_on dependency. If multiple remain, callers should
    use named lookups instead.
    """
    inputs = ctx.pipeline_def.get("inputs", {}) or {}
    depends_on = ctx.pipeline_def.get("depends_on", {}) or {}
    candidates = [k for k in inputs.keys() if k not in depends_on]
    if not candidates:
        raise RuntimeError(
            f"Pipeline '{ctx.pipeline_name}' has no non-dependency inputs; "
            "cannot determine the primary raw input."
        )
    if len(candidates) > 1:
        raise RuntimeError(
            f"Pipeline '{ctx.pipeline_name}' has multiple raw inputs "
            f"{candidates}; specify which is primary via stage definition "
            "(this pipeline likely needs a script_main stage instead)."
        )
    return _input_path(ctx, candidates[0])


def _render_args(args: dict[str, object], args_map: dict[str, str]) -> list[str]:
    """Translate {yaml_key: value} into ['--flag', 'value', ...] using args_map.

    Booleans render as bare flags (omitted when false). Multi-valued
    args (list/tuple) repeat the flag. Unknown keys raise.
    """
    out: list[str] = []
    for key, value in (args or {}).items():
        flag = args_map.get(key)
        if flag is None:
            raise KeyError(f"Stage arg '{key}' has no mapping in stage's args_map.")
        if isinstance(value, bool):
            if value:
                out.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            for v in value:
                out.extend([flag, str(v)])
            continue
        out.extend([flag, str(value)])
    return out


def _stage_output_path(ctx: RunContext, idx: int, stage_name: str) -> Path:
    # 1-based numbering matches the `primary_file` values in inputs.yaml
    # and reads naturally when listing the per-pipeline output directory.
    return ctx.pipeline_out_dir / f"{idx + 1:02d}_{stage_name}.tsv"


def _run_subprocess(cmd: list[str], env_extra: dict[str, str] | None = None) -> None:
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    print(f"[run] $ {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, cwd=REPO_ROOT)
    if proc.returncode != 0:
        raise SystemExit(f"Stage failed (exit {proc.returncode}): {cmd[0:3]}...")


def _run_positional(spec: dict, stage: stages_module.Stage, ctx: RunContext, idx: int) -> Path:
    if idx == 0:
        input_path = _pipeline_primary_input(ctx)
    else:
        assert ctx.prev_output is not None
        input_path = ctx.prev_output
    output_path = _stage_output_path(ctx, idx, stage.name)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [sys.executable, str(REPO_ROOT / stage.script)]
    if stage.accepts_input_dir and input_path.is_dir():
        cmd.extend(["-i", str(input_path), "-o", str(output_path)])
    else:
        cmd.extend([str(input_path), str(output_path)])
    cmd.extend(_render_args(spec.get("args") or {}, stage.args_map))

    _run_subprocess(cmd)
    return output_path


def _run_flagged(spec: dict, stage: stages_module.Stage, ctx: RunContext, idx: int) -> Path:
    if idx == 0:
        # First stage: the script's primary input flag points at the
        # pipeline-level primary raw input.
        primary_value = _pipeline_primary_input(ctx)
    else:
        assert ctx.prev_output is not None
        primary_value = ctx.prev_output
    output_path = _stage_output_path(ctx, idx, stage.name)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [sys.executable, str(REPO_ROOT / stage.script)]
    if stage.primary_input_flag:
        cmd.extend([stage.primary_input_flag, str(primary_value)])
    if stage.output_flag:
        cmd.extend([stage.output_flag, str(output_path)])
    for logical_name, flag in stage.input_flags.items():
        cmd.extend([flag, str(_input_path(ctx, logical_name))])
    cmd.extend(_render_args(spec.get("args") or {}, stage.args_map))

    _run_subprocess(cmd)
    return output_path


def _run_script_main(spec: dict, stage: stages_module.Stage, ctx: RunContext, idx: int) -> Path:
    env_extra: dict[str, str] = {}

    # Inputs: either pipeline-level (manifest) or an env var produced by a
    # previous stage (e.g. uksi_export reads UKSI_DB set by uksi_import).
    for logical_name, env_var in stage.env_inputs.items():
        if logical_name in ctx.manifest.get("paths", {}):
            env_extra[env_var] = str(_input_path(ctx, logical_name))
        elif env_var in ctx.env_state:
            env_extra[env_var] = ctx.env_state[env_var]
        else:
            raise KeyError(
                f"Stage '{stage.name}' input '{logical_name}' (env {env_var}) "
                "is neither a pipeline input nor produced by a prior stage."
            )

    # Outputs: every entry becomes a path under pipeline_out_dir and is
    # exported as an env var. Track them so later stages can consume.
    ctx.pipeline_out_dir.mkdir(parents=True, exist_ok=True)
    primary_output: Path | None = None
    for key, (env_var, filename) in stage.env_outputs.items():
        path = ctx.pipeline_out_dir / filename
        env_extra[env_var] = str(path)
        if key == stage.primary_output_key:
            primary_output = path

    # Carry forward (and into the process env) for later stages.
    ctx.env_state.update(env_extra)

    # Set env, then run the script as a file in a subprocess so its
    # module-level path globals are computed against the env we just set.
    # `stage.module` here is a dotted path mirroring the file location
    # (e.g. uksi_processing.uksi_db.uksi_import). Translate to a path.
    script_rel = Path(*stage.module.split(".")).with_suffix(".py")
    cmd = [sys.executable, str(REPO_ROOT / script_rel)]
    _run_subprocess(cmd, env_extra=env_extra)

    if primary_output is None:
        raise RuntimeError(
            f"Stage '{stage.name}' did not declare a primary_output_key."
        )
    return primary_output


def run_pipeline(pipeline_name: str, data_dir: Path, out_dir: Path) -> Path:
    cfg = _load_inputs()
    pipelines = cfg.get("pipelines") or {}
    if pipeline_name not in pipelines:
        raise SystemExit(f"Pipeline '{pipeline_name}' not in inputs.yaml. Known: {sorted(pipelines)}")
    pipeline_def = pipelines[pipeline_name]

    ctx = RunContext(
        pipeline_name=pipeline_name,
        pipeline_def=pipeline_def,
        data_dir=data_dir,
        out_dir=out_dir,
        manifest=_load_manifest(data_dir / pipeline_name),
    )
    ctx.pipeline_out_dir.mkdir(parents=True, exist_ok=True)

    stage_specs = pipeline_def.get("stages") or []
    if not stage_specs:
        raise SystemExit(f"Pipeline '{pipeline_name}' has no stages defined.")

    for idx, spec in enumerate(stage_specs):
        stage = stages_module.get(spec["name"])
        print(
            f"[run] === stage {idx + 1}/{len(stage_specs)}: {stage.name} "
            f"({stage.invocation}) ===",
            flush=True,
        )
        if stage.invocation == "positional_inout":
            ctx.prev_output = _run_positional(spec, stage, ctx, idx)
        elif stage.invocation == "flagged":
            ctx.prev_output = _run_flagged(spec, stage, ctx, idx)
        elif stage.invocation == "script_main":
            ctx.prev_output = _run_script_main(spec, stage, ctx, idx)
        else:
            raise SystemExit(f"Unknown invocation '{stage.invocation}' for stage '{stage.name}'.")

    final_output = ctx.prev_output
    assert final_output is not None

    _update_state(pipeline_name, ctx.manifest, final_output)
    print(f"[run] complete. Final output: {final_output}")
    return final_output


def _update_state(pipeline_name: str, manifest: dict, final_output: Path) -> None:
    state = json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}
    state.setdefault("schema_version", 2)
    state.setdefault("pipelines", {})
    state["pipelines"][pipeline_name] = {
        "input_versions": manifest.get("versions", {}),
        "last_run": {
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "final_output": str(final_output.relative_to(REPO_ROOT)) if final_output.is_relative_to(REPO_ROOT) else str(final_output),
        },
    }
    STATE_FILE.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pipeline", required=True, help="Pipeline name (key under `pipelines:` in inputs.yaml).")
    parser.add_argument("--data-dir", type=Path, default=Path(os.environ.get("DATA_DIR", REPO_ROOT / "data")))
    parser.add_argument("--out-dir", type=Path, default=Path(os.environ.get("OUT_DIR", REPO_ROOT / "final_result")))
    args = parser.parse_args()

    run_pipeline(args.pipeline, args.data_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
