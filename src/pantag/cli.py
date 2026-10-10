"""Configure Pan-TAG and run selected downstream workflows."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

from . import __version__


class ConfigureError(ValueError):
    """Raised when CLI configuration input is invalid."""


def positive_int(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {value!r}")
    return result


def nonnegative_int(value: str) -> int:
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError(f"expected an integer >= 0, got {value!r}")
    return result


def positive_float(value: str) -> float:
    result = float(value)
    if result <= 0:
        raise argparse.ArgumentTypeError(f"expected a number > 0, got {value!r}")
    return result


def fraction(value: str) -> float:
    result = float(value)
    if not 0 <= result <= 1:
        raise argparse.ArgumentTypeError(f"expected a number from 0 to 1, got {value!r}")
    return result


def percentage(value: str) -> float:
    result = float(value)
    if not 0 <= result <= 100:
        raise argparse.ArgumentTypeError(f"expected a number from 0 to 100, got {value!r}")
    return result


def comma_list(value: str) -> list[str]:
    result = [item.strip() for item in value.split(",") if item.strip()]
    if not result:
        raise argparse.ArgumentTypeError("expected a non-empty comma-separated list")
    return result


def _read_samples(path: Path, annotation_mode: str) -> list[dict[str, str]]:
    if not path.is_file():
        raise ConfigureError(f"--samples file does not exist: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["sample_id", "fasta", "gff"]
        if reader.fieldnames != expected:
            raise ConfigureError(
                "--samples header must contain exactly, in order: sample_id, fasta, gff"
            )
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    if not rows:
        raise ConfigureError("--samples must contain a reference sample row")

    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        sample_id = row["sample_id"]
        if not sample_id:
            raise ConfigureError(f"--samples line {line_number}: empty sample_id")
        if sample_id in seen:
            raise ConfigureError(f"--samples contains duplicate sample_id: {sample_id}")
        seen.add(sample_id)
        fasta = row["fasta"]
        if not fasta:
            raise ConfigureError(f"--samples line {line_number}: empty FASTA path")
        if not Path(fasta).is_file():
            raise ConfigureError(f"--samples FASTA does not exist for {sample_id}: {fasta}")
        gff = row["gff"]
        if gff and not Path(gff).is_file():
            raise ConfigureError(f"--samples GFF does not exist for {sample_id}: {gff}")

    if not rows[0]["gff"]:
        raise ConfigureError("the first reference row in --samples must provide a GFF")
    targets = rows[1:]
    if annotation_mode == "provided":
        missing = [row["sample_id"] for row in targets if not row["gff"]]
        if missing:
            raise ConfigureError(
                "--annotation-mode provided requires target GFFs: " + ", ".join(missing)
            )
    else:
        supplied = [row["sample_id"] for row in targets if row["gff"]]
        if supplied:
            raise ConfigureError(
                "--annotation-mode liftoff requires empty target GFF fields: "
                + ", ".join(supplied)
            )
    return rows


def set_nested(config: dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    target = config
    for key in path[:-1]:
        section = target.setdefault(key, {})
        if not isinstance(section, dict):
            raise ConfigureError(f"config template field {'.'.join(path[:-1])} must be a mapping")
        target = section
    target[path[-1]] = value


OVERRIDES: dict[str, tuple[str, ...]] = {
    "results_dir": ("results_dir",), "logs_dir": ("logs_dir",), "work_dir": ("work_dir",),
    "diamond_bin": ("software", "diamond"), "gffread_bin": ("software", "gffread"),
    "gene_feature_types": ("annotation", "gene_feature_types"),
    "pseudogene_feature_types": ("annotation", "pseudogene_feature_types"),
    "transcript_feature_types": ("annotation", "transcript_feature_types"),
    "gene_id_keys": ("annotation", "gene_id_keys"), "gene_name_keys": ("annotation", "gene_name_keys"),
    "transcript_id_keys": ("annotation", "transcript_id_keys"),
    "transcript_parent_gene_keys": ("annotation", "transcript_parent_gene_keys"),
    "biotype_keys": ("annotation", "biotype_keys"), "coding_filter_mode": ("annotation", "coding_filter_mode"),
    "coding_values": ("annotation", "coding_values"),
    "excluded_attribute_terms": ("annotation", "excluded_attribute_terms"),
    "pseudogene_values": ("annotation", "pseudogene_values"),
    "strip_id_prefixes": ("annotation", "strip_id_prefixes"),
    "fallback_gene_name_to_id": ("annotation", "fallback_gene_name_to_id"),
    "liftoff_environment": ("liftoff", "pixi_environment"),
    "liftoff_bin": ("liftoff", "executable"), "liftoff_threads": ("liftoff", "threads"),
    "diamond_threads": ("tag", "diamond_threads"),
    "diamond_max_target_seqs": ("tag", "diamond_max_target_seqs"),
    "blastp_identity": ("tag", "blastp_identity"),
    "blastp_query_coverage": ("tag", "blastp_query_coverage"),
    "blastp_subject_coverage": ("tag", "blastp_subject_coverage"),
    "pseudogene_enabled": ("pseudogene", "enabled"),
    "pseudogene_flank_bp": ("pseudogene", "search_flank_bp"),
    "blastx_identity": ("pseudogene", "blastx_identity"),
    "blastx_query_coverage": ("pseudogene", "blastx_query_coverage"),
    "blastx_subject_coverage": ("pseudogene", "blastx_subject_coverage"),
    "flank_gene_count": ("flanks", "gene_count"),
    "member_overlap_threshold": ("merge", "member_overlap_threshold"),
    "require_same_seqid": ("merge", "require_same_seqid"), "debug_outputs": ("outputs", "debug"),
    "refine_threads": ("refine", "threads"),
    "flank_overlap_threshold": ("refine", "flank_overlap_threshold"),
    "extract_threads": ("downstream", "extract", "threads"),
    "cluster_threads": ("downstream", "cluster", "threads"),
    "cluster_k": ("downstream", "cluster", "k"),
    "ubiquitous_threshold": ("downstream", "cluster", "ubiquitous_threshold"),
    "max_clusters": ("downstream", "cluster", "max_clusters"),
    "silhouette_tolerance": ("downstream", "cluster", "silhouette_tolerance"),
    "min_sequences": ("downstream", "cluster", "min_sequences"),
    "pcoa": ("downstream", "cluster", "pcoa"), "seed": ("downstream", "cluster", "seed"),
    "pcoa_jitter_frac": ("downstream", "cluster", "pcoa_jitter_frac"),
    "pcoa_point_size": ("downstream", "cluster", "pcoa_point_size"),
    "pcoa_alpha": ("downstream", "cluster", "pcoa_alpha"),
}

def configure(args: argparse.Namespace) -> Path:
    template_path = Path(args.template)
    if not template_path.is_file():
        raise ConfigureError(f"--template file does not exist: {template_path}")
    with template_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ConfigureError("--template must contain a YAML mapping")

    annotation_mode = args.annotation_mode
    if annotation_mode not in {"liftoff", "provided"}:
        raise ConfigureError(
            f"--annotation-mode/template annotation.mode has invalid value {annotation_mode!r}"
        )
    samples_path = Path(args.samples)
    rows = _read_samples(samples_path, annotation_mode)
    reference = rows[0]
    requested_reference = getattr(args, "reference_sample", None)
    if requested_reference and requested_reference != reference["sample_id"]:
        raise ConfigureError(
            f"--reference-sample {requested_reference!r} does not match the first --samples row "
            f"{reference['sample_id']!r}"
        )

    config["reference"] = {"sample_id": reference["sample_id"], "fasta": reference["fasta"], "gff": reference["gff"]}
    config["samples"] = str(samples_path)
    set_nested(config, ("annotation", "mode"), annotation_mode)
    liftoff_config = config.get("liftoff")
    if isinstance(liftoff_config, dict):
        liftoff_config.pop("feature_types", None)
    for destination, path in OVERRIDES.items():
        value = getattr(args, destination, None)
        if value is not None:
            set_nested(config, path, value)

    annotation = config.get("annotation")
    if not isinstance(annotation, dict):
        raise ConfigureError("template annotation must be a mapping")
    required_annotation_keys = (
        "gene_feature_types", "pseudogene_feature_types", "transcript_feature_types",
        "gene_id_keys", "gene_name_keys", "transcript_id_keys",
        "transcript_parent_gene_keys", "biotype_keys", "coding_filter_mode",
        "coding_values", "excluded_attribute_terms", "pseudogene_values",
        "strip_id_prefixes", "fallback_gene_name_to_id",
    )
    missing_annotation = [
        key for key in required_annotation_keys if key not in annotation
    ]
    if missing_annotation:
        raise ConfigureError(
            "template annotation is missing required config value(s): "
            + ", ".join(f"annotation.{key}" for key in missing_annotation)
        )
    for key in ("gene_feature_types", "pseudogene_feature_types"):
        configured = annotation[key]
        if isinstance(configured, str):
            configured = [configured]
        if not isinstance(configured, list) or not any(
            str(value).strip() for value in configured
        ):
            raise ConfigureError(
                f"annotation.{key} must contain at least one feature type"
            )
        annotation[key] = configured

    tag = config.get("tag")
    if not isinstance(tag, dict) or "diamond_max_target_seqs" not in tag:
        raise ConfigureError(
            "template is missing required config value: tag.diamond_max_target_seqs"
        )

    try:
        refine_threshold = float(config["refine"]["flank_overlap_threshold"])
    except (TypeError, ValueError) as error:
        raise ConfigureError("refine.flank_overlap_threshold must be a number from 0 to 1") from error
    if not 0 <= refine_threshold <= 1:
        raise ConfigureError("refine.flank_overlap_threshold must be between 0 and 1")
    for path, label in (
        (("tag", "diamond_max_target_seqs"), "tag.diamond_max_target_seqs"),
        (("downstream", "extract", "threads"), "downstream.extract.threads"),
        (("downstream", "cluster", "threads"), "downstream.cluster.threads"),
    ):
        value = nested_value(config, path)
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ConfigureError(f"{label} must be a positive integer")
    tolerance = nested_value(config, ("downstream", "cluster", "silhouette_tolerance"))
    if (
        not isinstance(tolerance, (int, float))
        or isinstance(tolerance, bool)
        or not 0 <= tolerance <= 1
    ):
        raise ConfigureError(
            "downstream.cluster.silhouette_tolerance must be between 0 and 1"
        )
    directories = [
        str(config[key]).strip()
        for key in ("results_dir", "logs_dir", "work_dir")
    ]
    if any(not value for value in directories) or len(set(directories)) != 3:
        raise ConfigureError(
            "--results-dir, --logs-dir and --work-dir must be three non-empty distinct paths"
        )
    executables = {
        "--diamond-bin": str(config["software"]["diamond"]),
        "--gffread-bin": str(config["software"]["gffread"]),
        "--liftoff-bin": str(config["liftoff"]["executable"]),
    }
    for option, executable in executables.items():
        candidate = Path(executable)
        if (candidate.is_absolute() or candidate.parent != Path(".")) and not candidate.is_file():
            raise ConfigureError(f"{option} path does not exist: {executable}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        yaml.safe_dump(config, handle, sort_keys=False, allow_unicode=True)
    return output_path


def add_bool(group: argparse._ArgumentGroup, name: str, dest: str, help_text: str) -> None:
    group.add_argument(name, dest=dest, action=argparse.BooleanOptionalAction, default=None, help=help_text)


def nested_value(config: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = config
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def display_default(value: Any) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, list):
        return "[" + ", ".join(str(item) for item in value) + "]"
    return str(value)


def append_config_defaults(parser: argparse.ArgumentParser, config: dict[str, Any]) -> None:
    excluded = {"samples", "annotation_mode", "pseudogene_enabled", "reference_sample"}
    for action in parser._actions:
        if action.dest in excluded or action.dest not in OVERRIDES:
            continue
        value = nested_value(config, OVERRIDES[action.dest])
        if value is not None:
            action.help = (action.help or "") + f" (config default: {display_default(value)})"


def build_parser() -> argparse.ArgumentParser:
    default_config_path = Path("config/config.yaml")
    default_config: dict[str, Any] = {}
    if default_config_path.is_file():
        with default_config_path.open(encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
            if isinstance(loaded, dict):
                default_config = loaded
    parser = argparse.ArgumentParser(prog="pantag configure", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)

    inputs = parser.add_argument_group("Input and sample options")
    inputs.add_argument("--samples", required=True, help="TSV with sample_id, fasta and gff; first row is reference (required)")
    inputs.add_argument("--reference-sample", help="Validate that this ID is the first samples.tsv row")
    inputs.add_argument("--template", default="config/config.yaml", help="Existing config template whose unspecified values are retained (default: config/config.yaml)")
    inputs.add_argument("--output", "--output-config", dest="output", default="config/config.yaml", help="Generated config.yaml path (default: config/config.yaml)")

    outputs = parser.add_argument_group("Output options")
    outputs.add_argument("--results-dir", help="Results directory")
    outputs.add_argument("--logs-dir", help="Log directory")
    outputs.add_argument("--work-dir", help="Intermediate work directory")
    add_bool(outputs, "--debug-outputs", "debug_outputs", "Publish debug mapping outputs")

    annotation = parser.add_argument_group("Reference and annotation options")
    annotation.add_argument("--annotation-mode", choices=("liftoff", "provided"), required=True, help="Annotate all targets or use all provided GFFs (required)")
    for option, destination, description in (
        ("--gene-feature-types", "gene_feature_types", "Comma-separated gene feature types"),
        ("--pseudogene-feature-types", "pseudogene_feature_types", "Comma-separated gene-level feature types eligible for pseudogene checking"),
        ("--transcript-feature-types", "transcript_feature_types", "Comma-separated transcript feature types"),
        ("--gene-id-keys", "gene_id_keys", "Ordered comma-separated gene ID attribute keys"),
        ("--gene-name-keys", "gene_name_keys", "Ordered comma-separated gene name keys"),
        ("--transcript-id-keys", "transcript_id_keys", "Ordered transcript ID keys"),
        ("--transcript-parent-gene-keys", "transcript_parent_gene_keys", "Ordered transcript parent-gene keys"),
        ("--biotype-keys", "biotype_keys", "Ordered biotype attribute keys"),
        ("--coding-values", "coding_values", "Biotype values treated as coding"),
        ("--excluded-attribute-terms", "excluded_attribute_terms", "Gene attribute terms to exclude"),
        ("--pseudogene-values", "pseudogene_values", "Biotype values treated as pseudogenes"),
        ("--strip-id-prefixes", "strip_id_prefixes", "ID prefixes removed during normalization"),
    ):
        annotation.add_argument(option, dest=destination, type=comma_list, help=description)
    annotation.add_argument("--coding-filter-mode", choices=("auto", "biotype", "translated_protein"), help="Coding-gene filter strategy")
    add_bool(annotation, "--fallback-gene-name-to-id", "fallback_gene_name_to_id", "Use gene ID when gene name is absent")
    annotation.add_argument("--liftoff-environment", help="Pixi environment containing Liftoff")
    annotation.add_argument("--liftoff-bin", help="Liftoff executable name or path")
    annotation.add_argument("--liftoff-threads", type=positive_int, help="Threads per Liftoff job")

    software = parser.add_argument_group("Software and TAG options")
    software.add_argument("--diamond-bin", help="DIAMOND executable name or path")
    software.add_argument("--gffread-bin", help="gffread executable name or path")
    software.add_argument("--diamond-threads", type=positive_int, help="Threads per DIAMOND job")
    software.add_argument("--diamond-max-target-seqs", type=positive_int, help="Maximum non-self DIAMOND targets reported per query")
    software.add_argument("--blastp-identity", type=percentage, help="Minimum blastp identity percentage")
    software.add_argument("--blastp-query-coverage", type=percentage, help="Minimum blastp query coverage percentage")
    software.add_argument("--blastp-subject-coverage", type=percentage, help="Minimum blastp subject coverage percentage")

    pseudo = parser.add_argument_group("Pseudogene and flank options")
    pseudo.add_argument(
        "--pseudogene-enabled",
        dest="pseudogene_enabled",
        action=argparse.BooleanOptionalAction,
        required=True,
        help="Explicitly enable or disable pseudogene recall (required)",
    )
    pseudo.add_argument("--pseudogene-flank-bp", type=nonnegative_int, help="Pseudogene search distance in bp")
    pseudo.add_argument("--blastx-identity", type=percentage, help="Minimum blastx identity percentage")
    pseudo.add_argument("--blastx-query-coverage", type=percentage, help="Minimum blastx query coverage percentage")
    pseudo.add_argument("--blastx-subject-coverage", type=percentage, help="Minimum blastx subject coverage percentage")
    pseudo.add_argument("--flank-gene-count", type=positive_int, help="Required genes on each TAG flank")

    merge_group = parser.add_argument_group("Merge options")
    merge_group.add_argument("--member-overlap-threshold", type=fraction, help="Minimum local/global member overlap")
    add_bool(merge_group, "--require-same-seqid", "require_same_seqid", "Require identical local/global seq_id before merging")

    refinement = parser.add_argument_group("Core refinement options")
    refinement.add_argument("--refine-threads", type=positive_int, help="Workers for core refine")
    refinement.add_argument("--flank-overlap-threshold", type=fraction, help="Minimum fraction of a candidate block's sample flanks matching global flanks")

    downstream = parser.add_argument_group("Optional downstream options")
    downstream.add_argument("--extract-threads", type=positive_int, help="Concurrent sample FASTA workers for sequence extraction")
    downstream.add_argument("--cluster-threads", type=positive_int, help="Workers for k-mer generation within each Pan-cluster")
    downstream.add_argument("--cluster-k", type=positive_int, help="k-mer length")
    downstream.add_argument("--ubiquitous-threshold", type=fraction, help="Ubiquitous k-mer frequency threshold")
    downstream.add_argument("--max-clusters", type=nonnegative_int, help="Maximum clusters; 0 selects automatically")
    downstream.add_argument("--silhouette-tolerance", type=fraction, help="Allowed loss from the best silhouette when preferring fewer groups")
    downstream.add_argument("--min-sequences", type=positive_int, help="Minimum sequences needed for clustering")
    add_bool(downstream, "--pcoa", "pcoa", "Generate PCoA output")
    downstream.add_argument("--seed", type=int, help="Random seed")
    downstream.add_argument("--pcoa-jitter-frac", type=nonnegative_int_or_float, help="PCoA jitter fraction")
    downstream.add_argument("--pcoa-point-size", type=positive_float, help="PCoA point size")
    downstream.add_argument("--pcoa-alpha", type=fraction, help="PCoA point alpha")

    append_config_defaults(parser, default_config)
    return parser


def build_downstream_parser(command: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=f"pantag {command}",
        description=(
            "Extract sequences for selected or all Pan-clusters."
            if command == "extract"
            else "Run sequence clustering for selected or all Pan-clusters."
        ),
    )
    parser.add_argument("--config", required=True, help="Run-specific config.yaml")
    targets = parser.add_mutually_exclusive_group(required=True)
    targets.add_argument(
        "--tag", action="append", dest="tags", metavar="GLOBAL_TAG_ID",
        help="Pan-cluster to process; repeat this option for multiple TAGs",
    )
    targets.add_argument("--all", action="store_true", dest="all_tags", help="Process all Pan-clusters")
    parser.add_argument("--cores", type=positive_int, required=True, help="Total Snakemake cores")
    parser.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Show the jobs Snakemake would run without executing them",
    )
    parser.add_argument(
        "--printshellcmds", "-p", action="store_true",
        help="Print shell commands for scheduled jobs",
    )
    return parser


def build_extract_parser() -> argparse.ArgumentParser:
    return build_downstream_parser("extract")


def build_cluster_parser() -> argparse.ArgumentParser:
    return build_downstream_parser("cluster")


def run_downstream(args: argparse.Namespace, target: str) -> int:
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        raise ConfigureError(f"--config file does not exist: {config_path}")
    project_root = Path(__file__).resolve().parents[2]
    snakefile = project_root / "workflow" / "Snakefile"
    if not snakefile.is_file():
        raise ConfigureError(f"Pan-TAG Snakefile does not exist: {snakefile}")
    snakemake = shutil.which("snakemake")
    if not snakemake:
        raise ConfigureError("snakemake is not available in the current Pixi environment")
    tags = list(dict.fromkeys(args.tags or []))
    invalid = [tag for tag in tags if not re.fullmatch(r"[A-Za-z0-9_.-]+", tag)]
    if invalid:
        raise ConfigureError(
            "Invalid Pan-cluster ID(s): " + ", ".join(invalid)
        )
    runtime_targets = "__ALL__" if args.all_tags else ",".join(tags)
    command = [
        snakemake, target,
        "--directory", str(Path.cwd().resolve()),
        "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(args.cores),
        "--rerun-incomplete",
        "--config",
        f"_pantag_downstream_targets={runtime_targets}",
    ]
    if target == "downstream_cluster":
        command[command.index("--config"):command.index("--config")] = [
            "--allowed-rules",
            "downstream_extract",
            "downstream_publish_sequence",
            "selected_clusters_ready",
            "downstream_cluster_tag",
            "build_structure_gene_index",
            "downstream_structure_tag",
            "downstream_structure_collapsed_tag",
            "downstream_cluster_summary",
            "downstream_cluster",
        ]
    if getattr(args, "dry_run", False):
        command.append("--dry-run")
    if getattr(args, "printshellcmds", False):
        command.append("--printshellcmds")
    return subprocess.run(command, check=False).returncode


def run_extract(args: argparse.Namespace) -> int:
    return run_downstream(args, "downstream_sequences")


def run_cluster(args: argparse.Namespace) -> int:
    return run_downstream(args, "downstream_cluster")


def build_root_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pantag", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("command", nargs="?", choices=("configure", "extract", "cluster"))
    return parser


def nonnegative_int_or_float(value: str) -> float:
    result = float(value)
    if result < 0:
        raise argparse.ArgumentTypeError(f"expected a number >= 0, got {value!r}")
    return result


def main(argv: list[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    if not raw or raw[0] in {"-h", "--help"}:
        build_root_parser().print_help()
        return 0
    if raw[0] in {"extract", "cluster"}:
        command = raw[0]
        parser = build_extract_parser() if command == "extract" else build_cluster_parser()
        args = parser.parse_args(raw[1:])
        try:
            return run_extract(args) if command == "extract" else run_cluster(args)
        except (ConfigureError, OSError) as error:
            parser.error(str(error))
    if raw[0] == "configure":
        raw = raw[1:]
    elif raw[0] == "--version":
        build_root_parser().parse_args(raw)
        return 0
    else:
        build_root_parser().error(f"unknown command: {raw[0]}")
    parser = build_parser()
    args = parser.parse_args(raw)
    try:
        output = configure(args)
    except (ConfigureError, OSError, yaml.YAMLError, ValueError) as error:
        parser.error(str(error))
    print(f"[Pan-TAG] Configuration written: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
