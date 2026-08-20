#!/usr/bin/env python3
"""Prepare Pan-TAG inputs.

This module currently provides only two commands:

1. validate-inputs
   Check config.yaml, samples.tsv and referenced FASTA/GFF files.

2. tag-inputs
   Generate the MCScanX four-column GFF, gene-name map and longest-protein
   FASTA. This command preserves the already validated pipeline1/pipeline2
   behavior.

Liftoff, gffread, DIAMOND and MCScanX are called separately by Snakemake.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import yaml
from Bio import SeqIO

from .gff import (
    gene_id,
    gene_name,
    load_schema,
    normalize_id,
    parse_attributes,
    transcript_ids,
    transcript_parent,
)


VERSION = "0.1.0"


class PanTagPrepareError(RuntimeError):
    """A user-facing Pan-TAG input or processing error."""


def ensure_parent(path: str | Path) -> Path:
    """Create an output file's parent directory and return the Path."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


# ---------------------------------------------------------------------------
# validate-inputs
# ---------------------------------------------------------------------------

def command_validate_inputs(args: argparse.Namespace) -> int:
    """Check the basic Pan-TAG configuration and input files."""
    config_path = Path(args.config)
    samples_path = Path(args.samples)

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise PanTagPrepareError("config.yaml must contain a YAML mapping.")

    reference = config.get("reference")
    if not isinstance(reference, dict):
        raise PanTagPrepareError("config.yaml is missing the 'reference' section.")

    problems: list[str] = []
    annotation = config.get("annotation", {})
    annotation_mode = str(annotation.get("mode", "")).strip().lower()
    if annotation_mode not in {"liftoff", "provided"}:
        problems.append("annotation.mode must be 'liftoff' or 'provided'")

    refine = config.get("refine", {})
    try:
        flank_overlap_threshold = float(refine["flank_overlap_threshold"])
    except KeyError:
        problems.append("Missing config value: refine.flank_overlap_threshold")
    except (TypeError, ValueError):
        problems.append("refine.flank_overlap_threshold must be a number from 0 to 1")
    else:
        if not 0 <= flank_overlap_threshold <= 1:
            problems.append("refine.flank_overlap_threshold must be between 0 and 1")

    reference_id = str(reference.get("sample_id", "")).strip()
    reference_fasta = str(reference.get("fasta", "")).strip()
    reference_gff = str(reference.get("gff", "")).strip()

    for key, value in (
        ("reference.sample_id", reference_id),
        ("reference.fasta", reference_fasta),
        ("reference.gff", reference_gff),
    ):
        if not value:
            problems.append(f"Missing config value: {key}")

    with samples_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = ["sample_id", "fasta", "gff"]
        if reader.fieldnames != required_columns:
            raise PanTagPrepareError(
                "samples.tsv must contain exactly three columns in this order: "
                "sample_id, fasta, gff"
            )

        samples: list[dict[str, str]] = []
        seen_sample_ids: set[str] = set()

        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            fasta = (row.get("fasta") or "").strip()
            gff = (row.get("gff") or "").strip()

            if not sample_id:
                problems.append(f"samples.tsv line {line_number}: empty sample_id")
                continue
            if line_number == 2 and sample_id != reference_id:
                problems.append(
                    "The first samples.tsv row must match reference.sample_id"
                )
            if sample_id in seen_sample_ids:
                problems.append(f"Duplicate sample_id: {sample_id}")
            seen_sample_ids.add(sample_id)
            if not fasta:
                problems.append(
                    f"samples.tsv line {line_number}: empty FASTA path"
                )
            if sample_id == reference_id:
                if fasta and fasta != reference_fasta:
                    problems.append(
                        "Reference FASTA differs between config.yaml and samples.tsv"
                    )
                if gff and gff != reference_gff:
                    problems.append(
                        "Reference GFF differs between config.yaml and samples.tsv"
                    )
                continue
            if annotation_mode == "liftoff" and gff:
                problems.append(
                    f"samples.tsv line {line_number}: gff must be empty in liftoff mode"
                )
            if annotation_mode == "provided" and not gff:
                problems.append(
                    f"samples.tsv line {line_number}: gff is required in provided mode"
                )

            samples.append({"sample_id": sample_id, "fasta": fasta, "gff": gff})

    files_to_check = [
        ("reference FASTA", reference_fasta),
        ("reference GFF3", reference_gff),
    ]
    files_to_check.extend(
        (f"sample FASTA [{row['sample_id']}]", row["fasta"])
        for row in samples
        if row["fasta"]
    )
    files_to_check.extend(
        (f"sample GFF3 [{row['sample_id']}]", row["gff"])
        for row in samples
        if row["gff"]
    )

    liftoff_samples = samples if annotation_mode == "liftoff" else []
    if annotation_mode == "liftoff":
        liftoff = config.get("liftoff")
        if not isinstance(liftoff, dict):
            problems.append(
                "config.yaml requires a 'liftoff' section because some samples have no GFF"
            )
        else:
            for key in ("feature_types",):
                value = str(liftoff.get(key, "")).strip()
                if not value:
                    problems.append(f"Missing config value: liftoff.{key}")
                elif not Path(value).is_file():
                    problems.append(f"Liftoff input does not exist: {value}")

    for label, filename in files_to_check:
        if filename and not Path(filename).is_file():
            problems.append(f"{label} does not exist: {filename}")

    report = {
        "status": "PASS" if not problems else "FAIL",
        "reference_sample": reference_id,
        "annotation_mode": annotation_mode,
        "target_sample_count": len(samples),
        "provided_gff_count": sum(bool(row["gff"]) for row in samples),
        "liftoff_sample_count": len(liftoff_samples),
        "problems": problems,
    }

    output = ensure_parent(args.output)
    output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    if problems:
        raise PanTagPrepareError(
            "Input validation failed:\n- " + "\n- ".join(problems)
        )

    print(
        f"[Pan-TAG] Input validation passed: {len(samples)} target sample(s)."
    )
    return 0


# ---------------------------------------------------------------------------
# tag-inputs
# ---------------------------------------------------------------------------
def command_tag_inputs(args: argparse.Namespace) -> int:
    schema = load_schema(args.config)
    gene_features = {value.lower() for value in schema.gene_feature_types}
    transcript_features = {value.lower() for value in schema.transcript_feature_types}
    genes: dict[str, dict[str, Any]] = {}
    transcript_to_gene: dict[str, str] = {}
    fallback_name_count = 0
    coding_values = {value.lower() for value in schema.coding_values}
    excluded_terms = {
        value.lower() for value in schema.excluded_attribute_terms if value
    }

    with Path(args.annotation_gff).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                raise PanTagPrepareError(
                    f"Malformed GFF row {line_number}: {line.rstrip()}"
                )
            feature = parts[2].lower()
            attributes = parse_attributes(parts[8])
            if feature in gene_features:
                identifier = gene_id(attributes, schema)
                if not identifier:
                    continue
                name, used_fallback = gene_name(attributes, schema, identifier)
                if not name:
                    continue
                fallback_name_count += int(used_fallback)
                genes[identifier] = {
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "gid": identifier,
                    "gname": name,
                    "biotypes": {
                        value.strip().lower()
                        for key in schema.biotype_keys
                        for value in attributes.get(key, "").split(",")
                        if value.strip()
                    },
                    "excluded_by_attribute": any(
                        term in parts[8].lower() for term in excluded_terms
                    ),
                }
            elif feature in transcript_features:
                parent = transcript_parent(attributes, schema)
                if parent:
                    for identifier in transcript_ids(attributes, schema):
                        transcript_to_gene[identifier] = parent

    coding_annotation_available = any(
        bool(gene["biotypes"] & coding_values) for gene in genes.values()
    )
    apply_biotype_filter = schema.coding_filter_mode == "biotype" or (
        schema.coding_filter_mode == "auto" and coding_annotation_available
    )
    excluded_attribute_genes = sum(
        bool(gene["excluded_by_attribute"]) for gene in genes.values()
    )
    noncoding_biotype_genes = sum(
        apply_biotype_filter and not bool(gene["biotypes"] & coding_values)
        for gene in genes.values()
    )

    best_proteins: dict[str, Any] = {}
    unmatched_proteins = 0
    with Path(args.raw_proteins).open(encoding="utf-8") as protein_handle:
        for record in SeqIO.parse(protein_handle, "fasta"):
            identifier = normalize_id(record.id, schema.strip_id_prefixes)
            parent = transcript_to_gene.get(identifier)
            if parent is None and identifier in genes:
                parent = identifier
            if parent not in genes:
                unmatched_proteins += 1
                continue
            clean_sequence = str(record.seq).replace("*", "").replace(".", "")
            if not clean_sequence:
                continue
            record.seq = record.seq.__class__(clean_sequence)
            record.id = parent
            record.description = ""
            if parent not in best_proteins or len(record.seq) > len(best_proteins[parent].seq):
                best_proteins[parent] = record

    translated_genes = [
        genes[identifier]
        for identifier in best_proteins
        if not genes[identifier]["excluded_by_attribute"]
        and (
            not apply_biotype_filter
            or bool(genes[identifier]["biotypes"] & coding_values)
        )
    ]
    translated_genes.sort(key=lambda item: (item["chrom"], item["start"], -item["end"]))
    retained_genes: list[dict[str, Any]] = []
    nested_skipped = 0
    last_by_chrom: dict[str, dict[str, Any]] = {}
    for gene in translated_genes:
        previous = last_by_chrom.get(gene["chrom"])
        if previous is not None and gene["end"] <= previous["end"]:
            nested_skipped += 1
            continue
        retained_genes.append(gene)
        last_by_chrom[gene["chrom"]] = gene

    retained_ids = {str(gene["gid"]) for gene in retained_genes}
    with (
        ensure_parent(args.output_mcscanx_gff).open("w", encoding="utf-8") as gff_handle,
        ensure_parent(args.output_gene_map).open("w", encoding="utf-8") as map_handle,
        ensure_parent(args.output_proteins).open("w", encoding="utf-8") as protein_handle,
    ):
        map_handle.write("Gene_ID\tGene_Name\n")
        for gene in retained_genes:
            gff_handle.write(
                f"{gene['chrom']}\t{gene['gid']}\t{gene['start']}\t{gene['end']}\n"
            )
            map_handle.write(f"{gene['gid']}\t{gene['gname']}\n")
        SeqIO.write(
            (best_proteins[str(gene["gid"])] for gene in retained_genes),
            protein_handle,
            "fasta",
        )

    retained = len(retained_genes)
    protein_count = len(retained_ids)

    if args.output_qc:
        qc_path = ensure_parent(args.output_qc)
        with qc_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["metric", "value"])
            writer.writerow(["sample_id", args.sample_id])
            writer.writerow(["retained_mcscanx_genes", retained])
            writer.writerow(["nested_genes_skipped", nested_skipped])
            writer.writerow(["translated_protein_genes", len(best_proteins)])
            writer.writerow(["coding_filter_mode", schema.coding_filter_mode])
            writer.writerow(["biotype_filter_applied", str(apply_biotype_filter).lower()])
            writer.writerow(["noncoding_biotype_genes_filtered", noncoding_biotype_genes])
            writer.writerow(["excluded_attribute_genes", excluded_attribute_genes])
            writer.writerow(["gene_name_fallback_to_id", fallback_name_count])
            writer.writerow(["unmatched_protein_records", unmatched_proteins])
            writer.writerow(["longest_proteins_written", protein_count])

    print(
        f"[Pan-TAG] Prepared TAG inputs for {args.sample_id}: "
        f"{retained} MCScanX gene(s), {protein_count} protein(s)."
    )
    return 0


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pantag.prepare",
        description="Validate Pan-TAG inputs and prepare MCScanX input files.",
    )
    parser.add_argument("--version", action="version", version=VERSION)

    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser(
        "validate-inputs",
        help="Check config.yaml, samples.tsv and referenced input files.",
    )
    validate.add_argument("--config", required=True)
    validate.add_argument("--samples", required=True)
    validate.add_argument("--output", required=True)
    validate.set_defaults(func=command_validate_inputs)

    tag_inputs = subparsers.add_parser(
        "tag-inputs",
        help="Generate MCScanX GFF/map and longest-protein FASTA.",
    )
    tag_inputs.add_argument("--sample-id", required=True)
    tag_inputs.add_argument("--config", required=True)
    tag_inputs.add_argument("--annotation-gff", required=True)
    tag_inputs.add_argument("--raw-proteins", required=True)
    tag_inputs.add_argument("--output-mcscanx-gff", required=True)
    tag_inputs.add_argument("--output-gene-map", required=True)
    tag_inputs.add_argument("--output-proteins", required=True)
    tag_inputs.add_argument("--output-qc")
    tag_inputs.set_defaults(func=command_tag_inputs)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        return int(args.func(args))
    except PanTagPrepareError as exc:
        print(f"[Pan-TAG ERROR] {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"[Pan-TAG ERROR] File not found: {exc.filename}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
