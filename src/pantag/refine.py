#!/usr/bin/env python3
"""Core single-pass refinement of candidate Pan-TAGs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import pysam

from pantag.gff import AnnotationSchema, gene_id, gene_name, load_schema, parse_attributes


class RefineError(RuntimeError):
    """User-facing refine failure."""


@dataclass(frozen=True)
class SampleAsset:
    sample_id: str
    fasta: str
    gff: str
    fasta_index: str = ""


@dataclass(frozen=True)
class GeneLite:
    gene_id: str
    gene_name: str
    cleaned_name: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class TagHit:
    gene_id: str
    gene_name: str
    cleaned_name: str
    seqid: str
    start: int
    end: int
    strand: str
    gene_order_index: int


COUNT_FIELDS = ["candidate_global_tag_id", "sample_id", "matrix_cell", "status"]
SPAN_FIELDS = [
    "candidate_global_tag_id", "sample_id", "valid_block_count", "seqids",
    "block_coordinates", "refined_count", "reason",
]
MANIFEST_FIELDS = [
    "record_id", "global_tag_id", "sample_id", "gene_id", "gene_name",
    "seqid", "start", "end", "strand", "copy_index", "block_id",
    "block_start", "block_end", "refined_count", "qc_status",
    "extraction_eligible",
]

_GENE_TO_CLUSTERS: dict[str, tuple[str, ...]] | None = None
_CLUSTER_TO_FLANKS: dict[str, frozenset[str]] | None = None
_FLANK_COUNT = 0
_FLANK_THRESHOLD = 0.0
_SCHEMA: AnnotationSchema | None = None


def clean_name(value: str) -> str:
    return re.sub(r"(_\d+|\.\d+|\(\d+\))+$", "", (value or "").strip())


def split_values(value: str) -> list[str]:
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:+-]+", "_", str(value).strip()) or "NA"


def atomic_tsv(path: str | Path, fields: list[str], rows) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_candidates(path: str):
    required = {
        "candidate_global_tag_id", "gene_names", "reference_seqid",
        "left_flank_genes", "right_flank_genes",
    }
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = required.difference(fields)
        if missing:
            raise RefineError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        order: list[str] = []
        rows: dict[str, dict[str, str]] = {}
        gene_map: dict[str, set[str]] = defaultdict(set)
        flank_map: dict[str, frozenset[str]] = {}
        for row in reader:
            candidate_id = row["candidate_global_tag_id"].strip()
            if not candidate_id or candidate_id in rows:
                raise RefineError(f"Invalid or duplicate candidate ID: {candidate_id!r}")
            seeds = {clean_name(value) for value in split_values(row["gene_names"])}
            seeds.discard("")
            if not seeds:
                raise RefineError(f"{candidate_id} has no usable gene_names")
            flanks = {
                clean_name(value)
                for key in ("left_flank_genes", "right_flank_genes")
                for value in split_values(row[key])
            }
            flanks.discard("")
            order.append(candidate_id)
            rows[candidate_id] = row
            flank_map[candidate_id] = frozenset(flanks)
            for seed in seeds:
                gene_map[seed].add(candidate_id)
    return fields, order, rows, {
        name: tuple(sorted(ids)) for name, ids in gene_map.items()
    }, flank_map


def read_samples(
    path: str, reference_id: str, reference_fasta: str, reference_gff: str,
    target_gff_dir: str, annotation_mode: str,
) -> list[SampleAsset]:
    assets: list[SampleAsset] = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "fasta", "gff"]:
            raise RefineError(f"{path} must contain exactly sample_id, fasta and gff")
        for row_number, row in enumerate(reader, 1):
            sample = (row.get("sample_id") or "").strip()
            if row_number == 1 and sample != reference_id:
                raise RefineError("The first samples.tsv row must be the reference")
            supplied_gff = (row.get("gff") or "").strip()
            gff = (
                reference_gff if sample == reference_id else
                supplied_gff if annotation_mode == "provided" else
                str(Path(target_gff_dir) / sample / f"{sample}.liftoff.gff3")
            )
            fasta = reference_fasta if sample == reference_id else (row.get("fasta") or "").strip()
            assets.append(SampleAsset(sample, fasta, gff))
    if not assets:
        raise RefineError("samples.tsv is empty")
    ids = [asset.sample_id for asset in assets]
    if len(ids) != len(set(ids)):
        raise RefineError("sample_id values must be unique")
    return assets


def initialize_worker(gene_map, flank_map, flank_count, threshold, schema) -> None:
    global _GENE_TO_CLUSTERS, _CLUSTER_TO_FLANKS, _FLANK_COUNT, _FLANK_THRESHOLD, _SCHEMA
    _GENE_TO_CLUSTERS = gene_map
    _CLUSTER_TO_FLANKS = flank_map
    _FLANK_COUNT = flank_count
    _FLANK_THRESHOLD = threshold
    _SCHEMA = schema


def read_gff_once(path: str, schema: AnnotationSchema, gene_map):
    chrom_genes: dict[str, list[GeneLite]] = defaultdict(list)
    allowed = {value.lower() for value in schema.gene_feature_types}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2].lower() not in allowed:
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError as exc:
                raise RefineError(f"{path}:{line_number}: invalid gene coordinates") from exc
            if start < 1 or end < start:
                raise RefineError(f"{path}:{line_number}: invalid gene interval {start}-{end}")
            attrs = parse_attributes(parts[8])
            identifier = gene_id(attrs, schema)
            name, _ = gene_name(attrs, schema, identifier)
            cleaned = clean_name(name)
            if identifier and name and cleaned:
                chrom_genes[parts[0]].append(
                    GeneLite(identifier, name, cleaned, start, end, parts[6])
                )

    sample_tag: dict[str, dict[str, list[TagHit]]] = defaultdict(lambda: defaultdict(list))
    for seqid, genes in chrom_genes.items():
        genes.sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
        for index, gene in enumerate(genes):
            for candidate_id in gene_map.get(gene.cleaned_name, ()):
                sample_tag[candidate_id][seqid].append(
                    TagHit(
                        gene.gene_id, gene.gene_name, gene.cleaned_name, seqid,
                        gene.start, gene.end, gene.strand, index,
                    )
                )
    return dict(chrom_genes), {
        candidate_id: dict(groups) for candidate_id, groups in sample_tag.items()
    }


def flank_score(genes, hits, cluster_flank, flank_count: int) -> float:
    first, last = hits[0].gene_order_index, hits[-1].gene_order_index
    flank_genes = (
        genes[max(0, first - flank_count):first]
        + genes[last + 1:last + 1 + flank_count]
    )
    sample_flank = {gene.cleaned_name for gene in flank_genes}
    return len(sample_flank & cluster_flank) / len(cluster_flank) if cluster_flank else -1.0


def scan_sample(asset: SampleAsset, workdir: str):
    if _GENE_TO_CLUSTERS is None or _CLUSTER_TO_FLANKS is None or _SCHEMA is None:
        raise RefineError("Refine worker was not initialized")
    started = time.monotonic()
    for path in (asset.fasta, asset.fasta_index, asset.gff):
        if not Path(path).is_file():
            raise RefineError(f"Missing input for {asset.sample_id}: {path}")
    chrom_genes, sample_tag = read_gff_once(asset.gff, _SCHEMA, _GENE_TO_CLUSTERS)
    count_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    span_rows: list[dict[str, object]] = []
    fasta = pysam.FastaFile(asset.fasta, filepath_index=asset.fasta_index)
    try:
        missing = sorted(set(chrom_genes).difference(fasta.references))
        if missing:
            raise RefineError(
                f"{asset.sample_id}: GFF seqid absent from FASTA: " + ", ".join(missing[:10])
            )
        for candidate_id in sorted(sample_tag):
            groups = sample_tag[candidate_id]
            valid = []
            for seqid in sorted(groups):
                hits = groups[seqid]
                hits.sort(key=lambda hit: hit.gene_order_index)
                score = flank_score(
                    chrom_genes[seqid], hits, _CLUSTER_TO_FLANKS[candidate_id], _FLANK_COUNT
                )
                if score >= _FLANK_THRESHOLD:
                    valid.append((
                        seqid, hits, min(hit.start for hit in hits),
                        max(hit.end for hit in hits), score,
                    ))
            if not valid:
                count_rows.append({
                    "candidate_global_tag_id": candidate_id, "sample_id": asset.sample_id,
                    "matrix_cell": "NA", "status": "no_valid_block",
                })
                continue
            total = sum(len(block[1]) for block in valid)
            if len(valid) > 1:
                count_rows.append({
                    "candidate_global_tag_id": candidate_id, "sample_id": asset.sample_id,
                    "matrix_cell": f"{total}(span)", "status": "span",
                })
                span_rows.append({
                    "candidate_global_tag_id": candidate_id, "sample_id": asset.sample_id,
                    "valid_block_count": len(valid),
                    "seqids": ",".join(block[0] for block in valid),
                    "block_coordinates": ";".join(
                        f"{block[0]}:{block[2]}-{block[3]}" for block in valid
                    ),
                    "refined_count": total,
                    "reason": "MULTIPLE_FLANK_SUPPORTED_BLOCKS",
                })
                continue
            seqid, hits, block_start, block_end, _ = valid[0]
            if "N" in fasta.fetch(seqid, block_start - 1, block_end).upper():
                count_rows.append({
                    "candidate_global_tag_id": candidate_id, "sample_id": asset.sample_id,
                    "matrix_cell": f"{total}(gap)", "status": "gap",
                })
                continue
            count_rows.append({
                "candidate_global_tag_id": candidate_id, "sample_id": asset.sample_id,
                "matrix_cell": str(total), "status": "pass",
            })
            for copy_index, hit in enumerate(hits, 1):
                manifest_rows.append({
                    "record_id": "", "global_tag_id": candidate_id,
                    "sample_id": asset.sample_id, "gene_id": hit.gene_id,
                    "gene_name": hit.gene_name, "seqid": seqid,
                    "start": hit.start, "end": hit.end,
                    "strand": hit.strand if hit.strand in {"+", "-"} else "+",
                    "copy_index": copy_index,
                    "block_id": f"{asset.sample_id}:{candidate_id}:block1",
                    "block_start": block_start, "block_end": block_end,
                    "refined_count": total, "qc_status": "pass",
                    "extraction_eligible": "true",
                })
    finally:
        fasta.close()

    root = Path(workdir)
    count_path = root / f"{asset.sample_id}.counts.tsv"
    manifest_path = root / f"{asset.sample_id}.manifest.tsv"
    span_path = root / f"{asset.sample_id}.span.tsv"
    atomic_tsv(count_path, COUNT_FIELDS, count_rows)
    atomic_tsv(manifest_path, MANIFEST_FIELDS, manifest_rows)
    atomic_tsv(span_path, SPAN_FIELDS, span_rows)
    elapsed = time.monotonic() - started
    print(
        f"[refine] DONE sample={asset.sample_id} hit_clusters={len(sample_tag)} "
        f"pass_records={len(manifest_rows)} spans={len(span_rows)} elapsed_sec={elapsed:.2f}",
        flush=True,
    )
    return asset.sample_id, str(count_path), str(manifest_path), str(span_path), {
        "hit_clusters": len(sample_tag), "manifest_records": len(manifest_rows),
        "span_clusters": len(span_rows), "elapsed_sec": elapsed,
    }


def read_counts(path: str) -> dict[str, str]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return {
            row["candidate_global_tag_id"]: row["matrix_cell"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def summarize(cells: list[str]) -> dict[str, object]:
    detected = sum(cell != "NA" for cell in cells)
    values = [int(cell) for cell in cells if cell.isdigit() and int(cell) > 0]
    mean = statistics.mean(values) if values else 0.0
    std = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "detected_sample_count": detected,
        "gap_free_sample_count": len(values),
        "gap_free_sample_fraction": f"{len(values) / detected if detected else 0.0:.2f}",
        "gene_count_mean": f"{mean:.2f}", "gene_count_std": f"{std:.2f}",
        "gene_count_cv": f"{std / mean if mean else 0.0:.2f}",
    }


def combine_manifests(path, assets, manifest_paths, id_map) -> None:
    def rows():
        for asset in assets:
            with Path(manifest_paths[asset.sample_id]).open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    candidate_id = row["global_tag_id"]
                    final_id = id_map.get(candidate_id)
                    if not final_id:
                        continue
                    row["global_tag_id"] = final_id
                    row["block_id"] = row["block_id"].replace(candidate_id, final_id, 1)
                    row["record_id"] = (
                        f"{safe_token(final_id)}|sample={safe_token(row['sample_id'])}"
                        f"|gene_name={safe_token(row['gene_name'])}"
                        f"|gene_id={safe_token(row['gene_id'])}|copy={row['copy_index']}"
                    )
                    yield row
    atomic_tsv(path, MANIFEST_FIELDS, rows())


def run(args: argparse.Namespace) -> None:
    schema = load_schema(args.config)
    assets = read_samples(
        args.samples, args.reference_sample, args.reference_fasta,
        args.reference_gff, args.target_gff_dir, schema.mode,
    )
    fasta_indexes = getattr(args, "fasta_indexes", None)
    if fasta_indexes is None:
        fasta_indexes = [f"{asset.fasta}.fai" for asset in assets]
    if len(fasta_indexes) != len(assets):
        raise RefineError("--fasta-indexes must provide one index per samples.tsv row")
    assets = [
        SampleAsset(asset.sample_id, asset.fasta, asset.gff, index)
        for asset, index in zip(assets, fasta_indexes)
    ]
    fields, order, rows, gene_map, flank_map = read_candidates(args.tags)
    if not order:
        raise RefineError("No candidate Pan-TAG records found")
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    executor = ProcessPoolExecutor if os.name == "posix" else ThreadPoolExecutor
    count_paths, manifest_paths, span_paths = {}, {}, {}
    with executor(
        max_workers=min(max(1, args.threads), len(assets)),
        initializer=initialize_worker,
        initargs=(gene_map, flank_map, args.flank_gene_count, args.flank_overlap_threshold, schema),
    ) as pool:
        futures = {pool.submit(scan_sample, asset, args.workdir): asset for asset in assets}
        for future in as_completed(futures):
            sample, counts, manifest, spans, _ = future.result()
            count_paths[sample], manifest_paths[sample], span_paths[sample] = counts, manifest, spans

    span_rows = []
    for asset in assets:
        with Path(span_paths[asset.sample_id]).open(encoding="utf-8", newline="") as handle:
            span_rows.extend(csv.DictReader(handle, delimiter="\t"))
    rank = {candidate_id: index for index, candidate_id in enumerate(order)}
    span_rows.sort(key=lambda row: (rank[row["candidate_global_tag_id"]], row["sample_id"]))
    removed = {row["candidate_global_tag_id"] for row in span_rows}
    id_map, map_rows = {}, []
    for candidate_id in order:
        if candidate_id in removed:
            map_rows.append({
                "candidate_global_tag_id": candidate_id, "global_tag_id": "",
                "status": "removed_span",
            })
        else:
            final_id = f"Pan_cluster_{len(id_map) + 1:04d}"
            id_map[candidate_id] = final_id
            map_rows.append({
                "candidate_global_tag_id": candidate_id,
                "global_tag_id": final_id, "status": "retained",
            })

    sample_ids = [asset.sample_id for asset in assets]
    counts_by_sample = {sample: read_counts(count_path) for sample, count_path in count_paths.items()}
    final_rows, matrix_rows = [], []
    for candidate_id in order:
        final_id = id_map.get(candidate_id)
        if not final_id:
            continue
        raw = dict(rows[candidate_id])
        raw.pop("candidate_global_tag_id", None)
        final_rows.append({"global_tag_id": final_id, **raw})
        cells = [counts_by_sample[sample].get(candidate_id, "NA") for sample in sample_ids]
        if any("(span)" in cell or cell == "0" for cell in cells):
            raise RefineError(f"Invalid final matrix state for {candidate_id}: {cells}")
        matrix_rows.append({
            "global_tag_id": final_id, "gene_names": raw["gene_names"],
            "reference_seqid": raw["reference_seqid"],
            **summarize(cells), **dict(zip(sample_ids, cells)),
        })

    atomic_tsv(
        args.output_tags,
        ["global_tag_id"] + [field for field in fields if field != "candidate_global_tag_id"],
        final_rows,
    )
    atomic_tsv(args.output_matrix, [
        "global_tag_id", "gene_names", "reference_seqid", "detected_sample_count",
        "gap_free_sample_count", "gap_free_sample_fraction", "gene_count_mean",
        "gene_count_std", "gene_count_cv", *sample_ids,
    ], matrix_rows)
    atomic_tsv(args.output_id_map, ["candidate_global_tag_id", "global_tag_id", "status"], map_rows)
    atomic_tsv(args.output_span_audit, SPAN_FIELDS, span_rows)
    combine_manifests(args.output_manifest, assets, manifest_paths, id_map)
    print(
        f"[refine] COMPLETE candidates={len(order)} retained={len(final_rows)} "
        f"removed_span={len(removed)}", flush=True,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in (
        "tags", "config", "samples", "reference-sample", "reference-fasta",
        "reference-gff", "target-gff-dir", "workdir", "output-tags",
        "output-matrix", "output-manifest", "output-id-map", "output-span-audit",
    ):
        parser.add_argument(f"--{option}", required=True)
    parser.add_argument("--flank-gene-count", type=int, default=10)
    parser.add_argument("--flank-overlap-threshold", type=float, default=0.20)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--fasta-indexes", nargs="+", required=True)
    return parser


def main() -> int:
    try:
        args = build_parser().parse_args()
        if args.flank_gene_count < 1 or args.threads < 1:
            raise RefineError("flank count and threads must be positive")
        if not 0 <= args.flank_overlap_threshold <= 1:
            raise RefineError("flank overlap threshold must be between 0 and 1")
        run(args)
    except (OSError, ValueError, RefineError, pysam.utils.SamtoolsError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
