#!/usr/bin/env python3
"""Index GFF genes once and draw per-sample TAG arrangements."""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow

from pantag.gff import gene_id, gene_name, load_schema, normalize_id, parse_attributes


class StructureError(RuntimeError):
    """User-facing structure plotting failure."""


@dataclass(frozen=True)
class Gene:
    seqid: str
    gene_id: str
    gene_name: str
    start: int
    end: int
    strand: str


def resolve_path(value: str, run_dir: Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (run_dir / path).resolve()


def normalize_gene_identifier(value: str, sample: str, prefixes: tuple[str, ...]) -> str:
    identifier = normalize_id(value, prefixes)
    if identifier.startswith(f"{sample}_"):
        identifier = identifier[len(sample) + 1:]
    return normalize_id(identifier, prefixes)


def read_samples(path: str, run_dir: Path, annotation_mode: str, liftoff_dir: Path):
    rows = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "fasta", "gff"]:
            raise StructureError("samples.tsv must contain exactly sample_id, fasta and gff")
        for index, row in enumerate(reader):
            sample = (row.get("sample_id") or "").strip()
            supplied = (row.get("gff") or "").strip()
            gff = (
                resolve_path(supplied, run_dir)
                if index == 0 or annotation_mode == "provided"
                else liftoff_dir / f"{sample}.liftoff.gff3"
            )
            rows.append((sample, gff))
    return rows


def read_matrix(path: str, tag: str, sample_ids: list[str]):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        row = next(
            (item for item in csv.DictReader(handle, delimiter="\t")
             if item.get("global_tag_id") == tag),
            None,
        )
    if row is None:
        raise StructureError(f"Global TAG absent from matrix: {tag}")
    selected = []
    for sample in sample_ids:
        value = (row.get(sample) or "").strip()
        if re.fullmatch(r"[1-9]\d*", value):
            selected.append((int(value), sample))
    return sorted(selected, key=lambda item: (item[0], item[1]))


def read_groups(path: str, tag: str, prefixes: tuple[str, ...]):
    required = {
        "global_tag_id", "sample_id", "gene_id", "gene_name",
        "sequence_group_id", "sequence_group_color",
    }
    mapping = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise StructureError("PCoA coordinates missing: " + ", ".join(sorted(missing)))
        for row in reader:
            if row["global_tag_id"] == tag:
                sample = row["sample_id"]
                identifier = normalize_gene_identifier(row["gene_id"], sample, prefixes)
                if identifier:
                    mapping[(sample, identifier)] = (
                        row["sequence_group_id"], row["sequence_group_color"]
                    )
    if not mapping:
        raise StructureError(f"No usable sequence-group records for {tag}")
    return mapping


def read_gff(path: Path, sample: str, schema):
    genes = {}
    allowed = set(schema.gene_feature_types) | set(schema.pseudogene_feature_types)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in allowed:
                continue
            attributes = parse_attributes(fields[8])
            identifier = normalize_gene_identifier(
                gene_id(attributes, schema), sample, schema.strip_id_prefixes
            )
            name, _ = gene_name(attributes, schema, identifier)
            if identifier and name:
                genes.setdefault(fields[0], []).append(Gene(
                    fields[0], identifier, name,
                    int(fields[3]), int(fields[4]), fields[6],
                ))
    for records in genes.values():
        records.sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
    return genes


def index_sample(sample: str, gff: Path, schema):
    rows = [
        (sample, seqid, order, gene.gene_id, gene.gene_name,
         gene.start, gene.end, gene.strand)
        for seqid, genes in read_gff(gff, sample, schema).items()
        for order, gene in enumerate(genes)
    ]
    return sample, rows


def build_index(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise StructureError("--threads must be >= 1")
    run_dir = Path(args.run_dir).resolve()
    schema = load_schema(args.config)
    assets = read_samples(
        args.samples, run_dir, schema.mode, Path(args.liftoff_dir).resolve()
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    connection = sqlite3.connect(temporary)
    try:
        connection.execute(
            "CREATE TABLE genes (sample_id TEXT, seqid TEXT, gene_order INTEGER, "
            "gene_id TEXT, gene_name TEXT, start INTEGER, end INTEGER, strand TEXT)"
        )
        workers = min(args.threads, len(assets))
        completed = 0
        with ThreadPoolExecutor(max_workers=workers) as executor:
            pending_assets = iter(assets)
            futures = set()
            for _ in range(workers):
                sample, gff = next(pending_assets)
                futures.add(executor.submit(index_sample, sample, gff, schema))
            while futures:
                finished, futures = wait(futures, return_when=FIRST_COMPLETED)
                for future in finished:
                    sample, rows = future.result()
                    connection.executemany(
                        "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows
                    )
                    connection.commit()
                    completed += 1
                    print(
                        f"[structure-index] completed={completed}/{len(assets)} "
                        f"sample={sample} genes={len(rows)}",
                        flush=True,
                    )
                    asset = next(pending_assets, None)
                    if asset is not None:
                        futures.add(executor.submit(index_sample, *asset, schema))
        connection.execute("CREATE INDEX genes_by_id ON genes(sample_id, gene_id)")
        connection.execute(
            "CREATE INDEX genes_by_order ON genes(sample_id, seqid, gene_order)"
        )
        connection.commit()
    finally:
        connection.close()
    temporary.replace(output)
    print(f"[structure-index] DONE samples={len(assets)} output={output}", flush=True)


def indexed_region(connection, sample: str, target_ids: set[str], flank: int):
    placeholders = ",".join("?" for _ in target_ids)
    rows = connection.execute(
        "SELECT seqid, gene_order, gene_id, gene_name, start, end, strand "
        f"FROM genes WHERE sample_id=? AND gene_id IN ({placeholders}) "
        "ORDER BY seqid, gene_order",
        (sample, *sorted(target_ids)),
    ).fetchall()
    if not rows:
        raise StructureError(f"No indexed target genes for {sample}")
    seqids = {row[0] for row in rows}
    if len(seqids) != 1:
        raise StructureError(f"Target genes occur on multiple seqids for {sample}")
    seqid = rows[0][0]
    orders = [row[1] for row in rows]
    first, last = min(orders), max(orders)

    def genes_between(start: int, end: int):
        if end < start:
            return []
        return [
            Gene(seqid, gene_id, name, gene_start, gene_end, strand)
            for gene_id, name, gene_start, gene_end, strand in connection.execute(
                "SELECT gene_id, gene_name, start, end, strand FROM genes "
                "WHERE sample_id=? AND seqid=? AND gene_order BETWEEN ? AND ? "
                "ORDER BY gene_order",
                (sample, seqid, start, end),
            )
        ]

    targets = [Gene(row[0], row[2], row[3], row[4], row[5], row[6]) for row in rows]
    return (
        genes_between(max(0, first - flank), first - 1),
        targets,
        genes_between(last + 1, last + flank),
    )


def draw_gene(axis, x: float, y: float, gene: Gene, color: str, label: str | None = None) -> None:
    dx = 0.8 if gene.strand == "+" else -0.8
    origin = x if gene.strand == "+" else x + 0.8
    axis.add_patch(FancyArrow(
        origin, y, dx, 0, width=0.18, head_width=0.35, head_length=0.18,
        length_includes_head=True, color=color,
    ))
    axis.text(x + 0.4, y + 0.28, label if label is not None else gene.gene_name,
              ha="center", va="bottom", fontsize=6, rotation=45)


def short_group_label(group_id: str) -> str:
    """Display only the SG suffix while retaining full IDs in data files."""
    match = re.search(r"(?:^|-)SG(\d+)$", str(group_id))
    return f"SG{match.group(1)}" if match else str(group_id)


def draw_collapsed(
    output: Path,
    compositions_output: Path,
    tag: str,
    records: list[tuple[int, str, list[Gene], list[Gene], dict[str, str]]],
) -> None:
    """Draw one row per identical (copy number, ordered group) composition."""
    grouped: dict[tuple[int, tuple[str, ...]], list[str]] = {}
    colors: dict[str, str] = {}
    for copy_number, sample, _left, targets, group_map in records:
        composition = tuple(group_map[gene.gene_id][0] for gene in targets)
        grouped.setdefault((copy_number, composition), []).append(sample)
        colors.update({group: color for group, color in group_map.values()})

    compositions_output.parent.mkdir(parents=True, exist_ok=True)
    with compositions_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["copy_number", "composition", "count", "frequency", "samples"])
        total = len(records)
        for (copy_number, composition), samples in sorted(
            grouped.items(), key=lambda item: (-item[0][0], sorted(item[1])[0])
        ):
            names = sorted(samples)
            writer.writerow([
                copy_number, ",".join(composition), len(names),
                f"{len(names) / total:.6f}", ",".join(names),
            ])

    figure, axis = plt.subplots(figsize=(14, max(4, len(grouped) * 0.7)))
    y_step = 1.5
    y = max(len(grouped) - 1, 0) * y_step
    max_x = 0.0
    used: dict[str, str] = {}
    total = len(records)
    for (copy_number, composition), samples in sorted(
        grouped.items(), key=lambda item: (-item[0][0], sorted(item[1])[0])
    ):
        representative = next(record for record in records if record[1] == sorted(samples)[0])
        _cn, _sample, _left, targets, group_map = representative
        count = len(samples)
        axis.text(-1.5, y, f"CN={copy_number}  n={count}  ({count / total:.1%})",
                  ha="right", va="center", fontsize=9)
        x = 0.0
        for gene, group in zip(targets, composition):
            color = next(
                (color for gene_id, (group_id, color) in group_map.items() if group_id == group),
                colors.get(group, "0.65"),
            )
            used[group] = color
            draw_gene(axis, x, y, gene, color, label=short_group_label(group))
            x += 1
        axis.plot([0, x], [y, y], color="black", linewidth=0.5, zorder=0)
        max_x = max(max_x, x)
        y -= y_step

    handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=8,
                      markerfacecolor=color, markeredgecolor="none", label=group)
               for group, color in sorted(used.items())]
    axis.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                frameon=False, fontsize=8, title="Sequence group")
    axis.set_xlim(-4, max_x + 1)
    axis.set_ylim(y - y_step, max(len(grouped), 1) * y_step)
    axis.set_title(f"{tag} (collapsed compositions)", fontsize=13, pad=20)
    axis.axis("off")
    figure.subplots_adjust(left=0.10, right=0.84, top=0.95, bottom=0.05)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)


def run(args: argparse.Namespace) -> None:
    if args.flank < 0:
        raise StructureError("--flank must be >= 0")
    if bool(args.collapsed_output) != bool(args.compositions_output):
        raise StructureError("--collapsed-output and --compositions-output must be provided together")
    run_dir = Path(args.run_dir).resolve()
    schema = load_schema(args.config)
    assets = read_samples(
        args.samples, run_dir, schema.mode, Path(args.liftoff_dir).resolve()
    )
    selected = read_matrix(args.matrix, args.global_tag_id, [row[0] for row in assets])
    groups = read_groups(args.coordinates, args.global_tag_id, schema.strip_id_prefixes)
    targets_by_sample = {}
    for sample, identifier in groups:
        targets_by_sample.setdefault(sample, set()).add(identifier)

    records = []
    connection = sqlite3.connect(
        f"file:{Path(args.gene_index).resolve()}?mode=ro", uri=True
    )
    try:
        for copy_number, sample in selected:
            target_ids = targets_by_sample.get(sample)
            if not target_ids:
                raise StructureError(f"No sequence-group records for {sample}")
            left, targets, right = indexed_region(connection, sample, target_ids, args.flank)
            group_map = {
                gene.gene_id: groups[(sample, gene.gene_id)][0]
                for gene in targets
            }
            color_map = {
                groups[(sample, gene.gene_id)][0]: groups[(sample, gene.gene_id)][1]
                for gene in targets
            }
            records.append((copy_number, sample, left, targets, right, group_map, color_map))
    finally:
        connection.close()

    if not args.output:
        if not args.collapsed_output:
            raise StructureError("at least one structure output is required")
        draw_collapsed(
            Path(args.collapsed_output), Path(args.compositions_output),
            args.global_tag_id,
            [
                (copy_number, sample, left, targets,
                 {gene.gene_id: (group_map[gene.gene_id], color_map[group_map[gene.gene_id]])
                  for gene in targets})
                for copy_number, sample, left, targets, _right, group_map, color_map in records
            ],
        )
        print(f"[structure] tag={args.global_tag_id} collapsed_only=true", flush=True)
        return

    figure, axis = plt.subplots(figsize=(16, max(4, len(records) * 0.7)))
    used_groups = {}
    flanks_plotted = False
    plotted = 0
    max_x = 0.0
    y_step = 1.5
    y = max(len(selected) - 1, 0) * y_step
    for copy_number, sample, left, targets, right, _group_map, _color_map in records:
            flanks_plotted = flanks_plotted or bool(left or right)
            axis.text(-1.5, y, f"CN={copy_number}  {sample}",
                      ha="right", va="center", fontsize=9)
            x = 0.0
            for gene in left:
                draw_gene(axis, x, y, gene, "0.65")
                x += 1
            x += 0.7
            for gene in targets:
                group_id, color = groups[(sample, gene.gene_id)]
                used_groups[group_id] = color
                draw_gene(axis, x, y, gene, color)
                x += 1
            x += 0.7
            for gene in right:
                draw_gene(axis, x, y, gene, "0.65")
                x += 1
            axis.plot([0, x], [y, y], color="black", linewidth=0.5, zorder=0)
            max_x = max(max_x, x)
            y -= y_step
            plotted += 1
    if not plotted:
        plt.close(figure)
        raise StructureError(f"No eligible sample structures could be drawn for {args.global_tag_id}")
    handles = [
        Line2D([0], [0], marker="s", linestyle="none", markersize=8,
               markerfacecolor=color, markeredgecolor="none", label=group)
        for group, color in sorted(used_groups.items())
    ]
    if flanks_plotted:
        handles.append(Line2D([0], [0], marker="s", linestyle="none", markersize=8,
                              markerfacecolor="0.65", markeredgecolor="none",
                              label="Flanking gene"))
    axis.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                frameon=False, fontsize=8, title="Sequence group")
    axis.set_xlim(-4, max_x + 1)
    axis.set_ylim(y - y_step, max(plotted, 1) * y_step)
    axis.set_title(args.global_tag_id, fontsize=13, pad=20)
    axis.axis("off")
    figure.subplots_adjust(left=0.10, right=0.84, top=0.95, bottom=0.05)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
    if args.collapsed_output:
        collapsed_records = [
            (copy_number, sample, left, targets,
             {gene.gene_id: (group_map[gene.gene_id], color_map[group_map[gene.gene_id]])
              for gene in targets})
            for copy_number, sample, left, targets, _right, group_map, color_map in records
        ]
        draw_collapsed(
            Path(args.collapsed_output), Path(args.compositions_output),
            args.global_tag_id, collapsed_records,
        )
    print(f"[structure] tag={args.global_tag_id} samples={plotted}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--global-tag-id", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--gene-index", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--liftoff-dir", required=True)
    parser.add_argument("--output")
    parser.add_argument("--collapsed-output")
    parser.add_argument("--compositions-output")
    parser.add_argument("--flank", type=int, default=3, help="Flanking genes shown per side")
    return parser


def build_index_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build the structure-plot GFF index")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--liftoff-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, required=True)
    return parser


def main() -> int:
    try:
        if len(sys.argv) > 1 and sys.argv[1] == "index":
            build_index(build_index_parser().parse_args(sys.argv[2:]))
        else:
            run(build_parser().parse_args())
    except (OSError, sqlite3.Error, ValueError, StructureError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
