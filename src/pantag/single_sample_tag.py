#!/usr/bin/env python3
"""Identify tandemly arrayed gene clusters in one sample.

This module packages the DIAMOND/MCScanX section of pipeline.sh together with
the validated behavior of pipeline3_tandem_processor.py.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


class SingleSampleTagError(RuntimeError):
    """A user-facing single-sample TAG error."""


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def run_command(command: list[str], *, stdout_path: Path | None = None) -> None:
    print("[Pan-TAG]", " ".join(command), file=sys.stderr)
    try:
        if stdout_path is None:
            subprocess.run(command, check=True)
        else:
            ensure_parent(stdout_path)
            with stdout_path.open("w", encoding="utf-8") as handle:
                subprocess.run(command, check=True, stdout=handle)
    except FileNotFoundError as exc:
        raise SingleSampleTagError(f"Executable not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise SingleSampleTagError(
            f"Command failed with exit code {exc.returncode}: {' '.join(command)}"
        ) from exc


def filter_blastp(
    raw_path: str | Path,
    output_path: str | Path,
    min_identity: float,
    min_query_coverage: float,
    min_subject_coverage: float,
) -> int:
    """Filter DIAMOND hits by identity and both sequence coverages."""
    retained = 0
    output = ensure_parent(output_path)
    with Path(raw_path).open("r", encoding="utf-8") as source, output.open(
        "w", encoding="utf-8"
    ) as target:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                raise SingleSampleTagError(
                    f"DIAMOND row {line_number} has {len(fields)} columns; expected 14"
                )
            if (
                float(fields[2]) >= min_identity
                and float(fields[12]) >= min_query_coverage
                and float(fields[13]) >= min_subject_coverage
            ):
                target.write(line if line.endswith("\n") else line + "\n")
                retained += 1
    return retained


def family_prefix(name: str, strategy: str, pattern: str = "") -> str:
    if strategy == "none":
        return str(name)
    if strategy == "first_digit_prefix":
        return re.split(r"\d", str(name))[0].rstrip("-")
    if strategy == "regex":
        match = re.search(pattern, str(name))
        if match is None:
            return str(name)
        return match.group(1) if match.groups() else match.group(0)
    raise SingleSampleTagError(f"Unknown family merge strategy: {strategy}")


def format_member_names(ids: list[str], id_to_name: dict[str, str]) -> str:
    names = [id_to_name.get(gene_id, gene_id) for gene_id in ids]
    counts = Counter(names)
    seen: set[str] = set()
    result: list[str] = []
    for name in names:
        if name in seen:
            continue
        result.append(f"{name}({counts[name]})" if counts[name] > 1 else name)
        seen.add(name)
    return ",".join(result)


def load_gene_map(path: str | Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"Gene_ID", "Gene_Name"}.issubset(
            reader.fieldnames
        ):
            raise SingleSampleTagError(
                f"{path} must contain Gene_ID and Gene_Name columns"
            )
        for row in reader:
            mapping[row["Gene_ID"]] = row["Gene_Name"]
    return mapping


def build_clusters(
    mcscanx_gff: str | Path,
    tandem_file: str | Path,
    gene_map: str | Path,
    output_file: str | Path,
    family_merge_strategy: str = "none",
    family_merge_regex: str = "",
) -> int:
    """Build local TAGs using pipeline3's connected-component behavior."""
    id_to_name = load_gene_map(gene_map)
    gene_coords: dict[str, tuple[str, int, int]] = {}
    chrom_genes: dict[str, list[str]] = defaultdict(list)

    with Path(mcscanx_gff).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise SingleSampleTagError(
                    f"Malformed MCScanX GFF row {line_number}: {line.rstrip()}"
                )
            chrom, gene_id = fields[0], fields[1]
            start, end = int(fields[2]), int(fields[3])
            gene_coords[gene_id] = (chrom, start, end)
            chrom_genes[chrom].append(gene_id)

    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda gene_id: gene_coords[gene_id][1])

    adjacency: dict[str, set[str]] = defaultdict(set)
    nodes: set[str] = set()
    tandem_path = Path(tandem_file)
    if tandem_path.exists():
        with tandem_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.strip().replace(" ", "").split(",")
                if len(fields) < 2:
                    continue
                left, right = fields[0], fields[1]
                if left in gene_coords and right in gene_coords:
                    nodes.update((left, right))
                    adjacency[left].add(right)
                    adjacency[right].add(left)

    initial: list[dict[str, object]] = []
    visited: set[str] = set()
    for node in sorted(
        nodes, key=lambda gene_id: (gene_coords[gene_id][0], gene_coords[gene_id][1])
    ):
        if node in visited:
            continue
        stack = [node]
        visited.add(node)
        members: list[str] = []
        while stack:
            current = stack.pop()
            members.append(current)
            for neighbor in sorted(adjacency[current], reverse=True):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        chrom = gene_coords[members[0]][0]
        ordered = sorted(members, key=lambda gene_id: gene_coords[gene_id][1])
        initial.append(
            {
                "ids": ordered,
                "chrom": chrom,
                "start": gene_coords[ordered[0]][1],
                "end": gene_coords[ordered[-1]][2],
            }
        )

    initial.sort(key=lambda cluster: (str(cluster["chrom"]), int(cluster["start"])))
    final: list[dict[str, object]] = []
    if initial:
        current = initial[0]
        for following in initial[1:]:
            current_ids = list(current["ids"])
            following_ids = list(following["ids"])
            current_prefixes = {
                family_prefix(
                    id_to_name.get(gene_id, gene_id),
                    family_merge_strategy,
                    family_merge_regex,
                )
                for gene_id in current_ids
            }
            following_prefixes = {
                family_prefix(
                    id_to_name.get(gene_id, gene_id),
                    family_merge_strategy,
                    family_merge_regex,
                )
                for gene_id in following_ids
            }
            if (
                current["chrom"] == following["chrom"]
                and not current_prefixes.isdisjoint(following_prefixes)
            ):
                merged_ids = sorted(
                    set(current_ids + following_ids),
                    key=lambda gene_id: gene_coords[gene_id][1],
                )
                current["ids"] = merged_ids
                current["start"] = min(
                    int(current["start"]), int(following["start"])
                )
                current["end"] = max(int(current["end"]), int(following["end"]))
            else:
                final.append(current)
                current = following
        final.append(current)

    output = ensure_parent(output_file)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "Cluster_ID",
                "Size",
                "Chrom",
                "Start",
                "End",
                "Length",
                "Upstream_5th_Gene",
                "Downstream_5th_Gene",
                "Members_Annotated",
                "Members_IDs",
            ]
        )
        for index, cluster in enumerate(final, start=1):
            ids = list(cluster["ids"])
            chrom = str(cluster["chrom"])
            ordered_chrom = chrom_genes[chrom]
            first_index = ordered_chrom.index(ids[0])
            last_index = ordered_chrom.index(ids[-1])
            upstream = (
                id_to_name.get(
                    ordered_chrom[first_index - 5],
                    ordered_chrom[first_index - 5],
                )
                if first_index >= 5
                else "START"
            )
            downstream = (
                id_to_name.get(
                    ordered_chrom[last_index + 5],
                    ordered_chrom[last_index + 5],
                )
                if last_index + 5 < len(ordered_chrom)
                else "END"
            )
            start, end = int(cluster["start"]), int(cluster["end"])
            writer.writerow(
                [
                    f"Cluster_{index}",
                    len(ids),
                    chrom,
                    start,
                    end,
                    end - start + 1,
                    upstream,
                    downstream,
                    format_member_names(ids, id_to_name),
                    ",".join(ids),
                ]
            )
    return len(final)


def command_run(args: argparse.Namespace) -> int:
    if args.output_final and not all((args.config, args.gff, args.fasta)):
        raise SingleSampleTagError(
            "--output-final requires --config, --gff and --fasta"
        )
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    database_prefix = workdir / f"{args.sample_id}_db"
    mcscanx_prefix = workdir / args.sample_id

    expected_gff = Path(f"{mcscanx_prefix}.gff")
    expected_blast = Path(f"{mcscanx_prefix}.blast")
    if Path(args.mcscanx_gff).resolve() != expected_gff.resolve():
        shutil.copyfile(args.mcscanx_gff, expected_gff)

    run_command(
        [
            args.diamond_bin,
            "makedb",
            "--in",
            args.proteins,
            "--db",
            str(database_prefix),
            "--quiet",
        ]
    )
    run_command(
        [
            args.diamond_bin,
            "blastp",
            "--query",
            args.proteins,
            "--db",
            str(database_prefix),
            "--out",
            args.output_blast_raw,
            "--evalue",
            "1e-5",
            "--max-hsps",
            "1",
            "--max-target-seqs",
            "5",
            "--threads",
            str(args.threads),
            "--quiet",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qcovhsp",
            "scovhsp",
        ]
    )
    retained = filter_blastp(
        args.output_blast_raw,
        args.output_blast,
        args.min_identity,
        args.min_query_coverage,
        args.min_subject_coverage,
    )
    if Path(args.output_blast).resolve() != expected_blast.resolve():
        shutil.copyfile(args.output_blast, expected_blast)

    mcscanx_stdout = workdir / f"{args.sample_id}.mcscanx.stdout.log"
    run_command([args.mcscanx_bin, str(mcscanx_prefix)], stdout_path=mcscanx_stdout)
    generated_tandem = Path(f"{mcscanx_prefix}.tandem")
    if not generated_tandem.exists():
        generated_tandem.touch()
    if generated_tandem.resolve() != Path(args.output_tandem).resolve():
        ensure_parent(args.output_tandem)
        shutil.copyfile(generated_tandem, args.output_tandem)

    clusters = build_clusters(
        args.mcscanx_gff,
        args.output_tandem,
        args.gene_map,
        args.output_cluster,
        args.family_merge_strategy,
        args.family_merge_regex,
    )
    if args.output_final:
        from .gff import load_schema
        from .retrieve_pseudogenes import (
            add_flanks,
            finalize_tag_metrics,
            load_gene_records,
            read_table,
            strict_qc,
            write_table,
        )

        schema = load_schema(args.config)
        coordinates, chrom_genes, positions, names, pseudogenes = load_gene_records(
            args.gff, schema
        )
        _, rows = read_table(args.output_cluster)
        for row in rows:
            row["Pseudos_Annotated"] = "None"
            row["Pseudos_IDs"] = "None"
        add_flanks(rows, chrom_genes, positions, names, args.flank_gene_count)
        finalize_tag_metrics(
            rows,
            coordinates,
            chrom_genes,
            {str(item["id"]) for item in pseudogenes},
        )
        retained_rows, qc_statuses = strict_qc(
            rows, coordinates, args.fasta, args.flank_gene_count
        )
        write_table(
            args.output_final,
            [
                "Cluster_ID", "Size", "Chromosome", "Start", "End", "Length",
                "Gene_Density_per_kb", "Members_Annotated", "Members_IDs",
                "Upstream_Genes", "Downstream_Genes",
            ],
            retained_rows,
        )
        print(
            "[Pan-TAG QC] "
            + " ".join(f"{key}={value}" for key, value in sorted(qc_statuses.items()))
        )
    print(
        f"[Pan-TAG] {args.sample_id}: {retained} blastp hit(s), "
        f"{clusters} local TAG(s)."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--mcscanx-bin", required=True)
    parser.add_argument("--mcscanx-gff", required=True)
    parser.add_argument("--gene-map", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-identity", type=float, default=20.0)
    parser.add_argument("--min-query-coverage", type=float, default=20.0)
    parser.add_argument("--min-subject-coverage", type=float, default=20.0)
    parser.add_argument("--diamond-bin", default="diamond")
    parser.add_argument(
        "--family-merge-strategy",
        choices=("none", "first_digit_prefix", "regex"),
        default="none",
    )
    parser.add_argument("--family-merge-regex", default="")
    parser.add_argument("--output-blast-raw", required=True)
    parser.add_argument("--output-blast", required=True)
    parser.add_argument("--output-tandem", required=True)
    parser.add_argument("--output-cluster", required=True)
    parser.add_argument("--output-final")
    parser.add_argument("--config")
    parser.add_argument("--gff")
    parser.add_argument("--fasta")
    parser.add_argument("--flank-gene-count", type=int, default=10)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return command_run(args)
    except (SingleSampleTagError, OSError, ValueError) as exc:
        print(f"[Pan-TAG ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
