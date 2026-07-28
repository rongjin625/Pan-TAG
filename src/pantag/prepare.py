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
        required_columns = {"sample_id", "fasta"}
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            raise PanTagPrepareError(
                "samples.tsv is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        samples: list[dict[str, str]] = []
        seen_sample_ids: set[str] = set()

        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            fasta = (row.get("fasta") or "").strip()

            if not sample_id:
                problems.append(f"samples.tsv line {line_number}: empty sample_id")
                continue
            if sample_id in seen_sample_ids:
                problems.append(f"Duplicate sample_id: {sample_id}")
            if not fasta:
                problems.append(
                    f"samples.tsv line {line_number}: empty FASTA path"
                )

            seen_sample_ids.add(sample_id)
            samples.append({"sample_id": sample_id, "fasta": fasta})

    if reference_id and reference_id in seen_sample_ids:
        problems.append(
            "The reference sample must not be repeated in samples.tsv."
        )

    files_to_check = [
        ("reference FASTA", reference_fasta),
        ("reference GFF3", reference_gff),
    ]
    files_to_check.extend(
        (f"sample FASTA [{row['sample_id']}]", row["fasta"])
        for row in samples
        if row["fasta"]
    )

    for label, filename in files_to_check:
        if filename and not Path(filename).is_file():
            problems.append(f"{label} does not exist: {filename}")

    report = {
        "status": "PASS" if not problems else "FAIL",
        "reference_sample": reference_id,
        "target_sample_count": len(samples),
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
# tag-inputs: faithful packaging of pipeline1 + pipeline2
# ---------------------------------------------------------------------------

def legacy_prepare_gff_filtered(
    input_gff: str | Path,
    output_gff: str | Path,
    map_out: str | Path,
) -> tuple[int, int]:
    """Faithful implementation of pipeline1_get_clean_gff.py."""
    genes: list[dict[str, Any]] = []

    with Path(input_gff).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                raise PanTagPrepareError(
                    f"Malformed GFF line in {input_gff}: {line.rstrip()}"
                )
            if parts[2] != "gene":
                continue

            attr = parts[8]
            if "gene_type=protein_coding" in attr and "readthrough_gene" not in attr:
                info: dict[str, str] = {}
                for item in attr.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info[key] = value

                gid = info.get("ID") or info.get("gene_id")
                gname = info.get("gene_name") or gid
                if gid:
                    genes.append(
                        {
                            "chrom": parts[0],
                            "start": int(parts[3]),
                            "end": int(parts[4]),
                            "gid": gid,
                            "gname": gname,
                        }
                    )

    genes.sort(key=lambda item: (item["chrom"], item["start"], -item["end"]))

    filtered_genes: list[dict[str, Any]] = []
    nested_skipped = 0
    if genes:
        last_gene = genes[0]
        filtered_genes.append(last_gene)
        for curr in genes[1:]:
            if curr["chrom"] != last_gene["chrom"]:
                filtered_genes.append(curr)
                last_gene = curr
                continue
            if curr["end"] <= last_gene["end"]:
                sys.stderr.write(
                    f"[Nested] 跳过被包含基因: {curr['gid']} "
                    f"(在 {last_gene['gid']} 内部)\n"
                )
                nested_skipped += 1
                continue
            filtered_genes.append(curr)
            last_gene = curr

    output_gff_path = ensure_parent(output_gff)
    map_out_path = ensure_parent(map_out)
    with (
        output_gff_path.open("w", encoding="utf-8") as gff_handle,
        map_out_path.open("w", encoding="utf-8") as map_handle,
    ):
        map_handle.write("Gene_ID\tGene_Name\n")
        for gene in filtered_genes:
            gff_handle.write(
                f"{gene['chrom']}\t{gene['gid']}\t{gene['start']}\t{gene['end']}\n"
            )
            map_handle.write(f"{gene['gid']}\t{gene['gname']}\n")

    return len(filtered_genes), nested_skipped


def legacy_filter_longest_and_rename(
    raw_faa: str | Path,
    filtered_gff: str | Path,
    original_gff: str | Path,
    output_faa: str | Path,
) -> tuple[int, int]:
    """Faithful implementation of pipeline2_filter_longest_and_rename.py."""
    whitelist_gids: set[str] = set()
    with Path(filtered_gff).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                raise PanTagPrepareError(
                    f"Malformed four-column GFF line: {line.rstrip()}"
                )
            whitelist_gids.add(parts[1])

    print(f"[*] 已加载白名单基因数量: {len(whitelist_gids)}")

    tid_to_gid: dict[str, str] = {}
    with Path(original_gff).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                raise PanTagPrepareError(
                    f"Malformed GFF line in {original_gff}: {line.rstrip()}"
                )
            if parts[2] in ["mRNA", "transcript"]:
                attr = parts[8]
                tid = ""
                gid = ""
                for item in attr.split(";"):
                    if item.startswith("ID="):
                        tid = item.split("=", 1)[1]
                    if item.startswith("gene_id="):
                        gid = item.split("=", 1)[1]
                if tid and gid in whitelist_gids:
                    tid_to_gid[tid] = gid

    gene_best_records: dict[str, Any] = {}
    for record in SeqIO.parse(str(raw_faa), "fasta"):
        tid = record.id
        gid = tid_to_gid.get(tid)
        if not gid:
            continue

        clean_seq = str(record.seq).replace("*", "").replace(".", "")
        record.seq = record.seq.__class__(clean_seq)
        seq_len = len(clean_seq)

        if gid not in gene_best_records or seq_len > len(gene_best_records[gid].seq):
            record.id = gid
            record.description = ""
            gene_best_records[gid] = record

    output_path = ensure_parent(output_faa)
    with output_path.open("w", encoding="utf-8") as handle:
        SeqIO.write(gene_best_records.values(), handle, "fasta")

    print(f"[*] 最终提取蛋白数量: {len(gene_best_records)}")
    return len(whitelist_gids), len(gene_best_records)


def command_tag_inputs(args: argparse.Namespace) -> int:
    retained, nested_skipped = legacy_prepare_gff_filtered(
        input_gff=args.annotation_gff,
        output_gff=args.output_mcscanx_gff,
        map_out=args.output_gene_map,
    )

    whitelist_count, protein_count = legacy_filter_longest_and_rename(
        raw_faa=args.raw_proteins,
        filtered_gff=args.output_mcscanx_gff,
        original_gff=args.annotation_gff,
        output_faa=args.output_proteins,
    )

    if args.output_qc:
        qc_path = ensure_parent(args.output_qc)
        with qc_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["metric", "value"])
            writer.writerow(["sample_id", args.sample_id])
            writer.writerow(["retained_mcscanx_genes", retained])
            writer.writerow(["nested_genes_skipped", nested_skipped])
            writer.writerow(["protein_whitelist_genes", whitelist_count])
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