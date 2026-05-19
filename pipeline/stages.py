"""Stage dispatch table for the multi-pipeline orchestrator.

Each entry describes how to invoke one of the existing analysis scripts.
Three invocation conventions cover every script we wrap:

- positional_inout:
    python <script> INPUT OUTPUT [flags...]
    The orchestrator passes the previous stage's output (or the
    pipeline's primary raw input, for the first stage) as INPUT,
    and writes OUTPUT to <out_dir>/<NN>_<stage>.tsv.

- flagged:
    python <script> --primary-flag INPUT --output OUTPUT [aux flags] [args]
    Used when scripts expect named flags rather than positional args.
    Auxiliary inputs (species_list, dtol_metadata, mito_metadata) are
    looked up in the per-pipeline manifest produced by fetch_zenodo.py.

- script_main:
    For scripts with no CLI today (uksi_import, uksi_export). The
    orchestrator sets module-level path globals via env vars, then
    imports the module and calls its main(). The matching env-var
    fallbacks live in the scripts themselves (see the uksi_db patches).
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Stage:
    name: str
    invocation: str
    # positional_inout / flagged share this:
    script: str | None = None
    # flagged-specific:
    primary_input_flag: str | None = None
    output_flag: str | None = None
    # Maps a logical input name (key in pipeline.inputs or pipeline.depends_on)
    # to the CLI flag the script expects. The first key in this dict is
    # treated as the "primary" record input for first-stage wiring of
    # flagged invocations that have no `primary_input_flag` set (i.e.
    # when the first input *is* the primary record stream).
    input_flags: dict[str, str] = field(default_factory=dict)
    # YAML arg key -> CLI flag (or "" for a positional). Bool values render
    # as bare `--flag` when truthy and are omitted when falsy.
    args_map: dict[str, str] = field(default_factory=dict)
    # If true, the script accepts a directory via -i/--input-dir when the
    # manifest entry points to a directory (e.g. extracted NCBI tarballs).
    accepts_input_dir: bool = False
    # script_main-specific:
    module: str | None = None
    # Logical pipeline input name -> env var to set before main().
    env_inputs: dict[str, str] = field(default_factory=dict)
    # Output key -> env var. Values are paths set to <out_dir>/<filename>.
    env_outputs: dict[str, tuple[str, str]] = field(default_factory=dict)
    # The output key (within env_outputs) that names the "primary" file
    # this stage produces, used to chain to subsequent stages.
    primary_output_key: str | None = None


# Stage name -> Stage definition. Keys here must match `stages[*].name`
# in pipeline/inputs.yaml.
STAGES: dict[str, Stage] = {
    # -----------------------------------------------------------------
    # UKSI species-list pipeline.
    # -----------------------------------------------------------------
    "uksi_import": Stage(
        name="uksi_import",
        invocation="script_main",
        module="uksi_processing.uksi_db.uksi_import",
        env_inputs={
            "uksi_names": "UKSI_NAMES_TSV",
            "uksi_taxa": "UKSI_TAXA_TSV",
            "pantheon": "UKSI_PANTHEON_TSV",
            "jncc": "UKSI_JNCC_TSV",
            "freshbase": "UKSI_FRESHBASE_TSV",
            "ukceh_freshwater": "UKSI_UKCEH_FRESHWATER_TSV",
        },
        env_outputs={
            "db": ("UKSI_DB", "uksi.db"),
            "log": ("UKSI_IMPORT_LOG", "uksi_import.log"),
        },
        primary_output_key="db",
    ),
    "uksi_export": Stage(
        name="uksi_export",
        invocation="script_main",
        module="uksi_processing.uksi_db.uksi_export",
        env_inputs={
            # `db` here is *not* a pipeline-level input; it's the output of
            # the previous uksi_import stage. The orchestrator stitches it
            # up by carrying over env vars produced by earlier stages.
        },
        env_outputs={
            "valid": ("UKSI_VALID_OUT", "uksi_species_export.tsv"),
            "invalid": ("UKSI_INVALID_OUT", "uksi_invalid_species_export.tsv"),
            "log": ("UKSI_EXPORT_LOG", "uksi_export.log"),
        },
        primary_output_key="valid",
    ),

    # -----------------------------------------------------------------
    # Per-source preprocessors and clustering.
    # All take INPUT OUTPUT positional args.
    # -----------------------------------------------------------------
    "bold_gene_extract": Stage(
        name="bold_gene_extract",
        invocation="positional_inout",
        script="bold_processing/bold_gene_extract/bold_gene_extract.py",
        args_map={"gene": "-g", "substring": "--substring"},
    ),
    "process_midori": Stage(
        name="process_midori",
        invocation="positional_inout",
        script="midori_processing/process_midori.py",
    ),
    "process_unite": Stage(
        name="process_unite",
        invocation="positional_inout",
        script="unite_processing/process_unite.py",
    ),
    "ncbi_gb_extract": Stage(
        name="ncbi_gb_extract",
        invocation="positional_inout",
        script="ncbi_processing/ncbi_gb_extract/ncbi_gb_extract.py",
        args_map={"gene": "-g"},
        accepts_input_dir=True,
    ),
    "otu_clustering": Stage(
        name="otu_clustering",
        invocation="positional_inout",
        script="otu_clustering/otu_clustering.py",
        args_map={
            "threshold": "-t",
            "threads": "--threads",
            "strand": "--strand",
            "min_length": "--min-length",
        },
    ),

    # -----------------------------------------------------------------
    # Flagged terminal stages (consume the species list).
    # -----------------------------------------------------------------
    "gap_analysis": Stage(
        name="gap_analysis",
        invocation="flagged",
        script="gap_analysis/gap_analysis.py",
        primary_input_flag="--records",
        output_flag="--output",
        input_flags={"species_list": "--species-list"},
        args_map={
            "workers": "--workers",
            "batch_size": "--batch-size",
            "cluster_column": "--cluster-column",
            "bold": "--bold",
            "marker": "--marker",
            "filter_kingdom": "--filter-kingdom",
            "no_filter_kingdom": "--no-filter-kingdom",
            "kingdom_list": "--kingdom-list",
        },
    ),
    "dtol_status": Stage(
        name="dtol_status",
        invocation="flagged",
        script="dtol_processing/dtol_status.py",
        primary_input_flag="--dtol-metadata",
        output_flag="--output",
        input_flags={"species_list": "--species-list"},
    ),
    "mito_status": Stage(
        name="mito_status",
        invocation="flagged",
        script="ena_processing/mito_status.py",
        primary_input_flag="--mito-metadata",
        output_flag="--output",
        input_flags={"species_list": "--species-list"},
    ),
}


def get(stage_name: str) -> Stage:
    try:
        return STAGES[stage_name]
    except KeyError as exc:
        raise KeyError(
            f"Unknown stage '{stage_name}'. Known stages: {sorted(STAGES)}"
        ) from exc
