#!/usr/bin/env python3
"""Cluster sequences within one Global TAG using canonical k-mers and Mash distance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.spatial.distance import squareform


_CODE = np.full(256, -1, dtype=np.int8)
for base, value in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("a", 0), ("c", 1), ("g", 2), ("t", 3)):
    _CODE[ord(base)] = value
_COMP = (3, 2, 1, 0)


class ClusterError(RuntimeError):
    """User-facing clustering failure."""


def read_fasta(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def read_metadata(path: str, global_tag_id: str) -> dict[str, dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"record_id", "global_tag_id", "sample_id", "gene_id", "gene_name", "copy_index"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ClusterError(f"Sequence records missing: {', '.join(sorted(missing))}")
        return {
            row["record_id"]: row
            for row in reader
            if row["global_tag_id"] == global_tag_id
        }


def seq_kmer_codes(sequence: str, k: int) -> set[int]:
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    code = reverse = length = 0
    kmers: set[int] = set()
    values = _CODE[np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)]
    for raw in values:
        value = int(raw)
        if value < 0:
            code = reverse = length = 0
            continue
        code = ((code << 2) | value) & mask
        reverse = (reverse >> 2) | (_COMP[value] << shift)
        length += 1
        if length >= k:
            kmers.add(min(code, reverse))
    return kmers


def build_kmer_matrix(
    records: list[tuple[str, str]], k: int, ubiquitous_threshold: float, threads: int,
) -> tuple[csr_matrix, int, int]:
    codes: list[np.ndarray] = []
    rows: list[np.ndarray] = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        kmer_sets = executor.map(
            lambda record: seq_kmer_codes(record[1], k), records
        )
        indexed_kmers = enumerate(kmer_sets)
        for index, values in indexed_kmers:
            if values:
                if k <= 31:
                    encoded = np.fromiter(values, dtype=np.uint64, count=len(values))
                else:
                    # Python integers preserve exact 2-bit encodings beyond 64 bits.
                    encoded = np.asarray(list(values), dtype=object)
                codes.append(encoded)
                rows.append(np.full(len(values), index, dtype=np.int32))
    if not codes:
        return csr_matrix((len(records), 0), dtype=np.float32), 0, 0
    all_codes = np.concatenate(codes)
    all_rows = np.concatenate(rows)
    unique, columns = np.unique(all_codes, return_inverse=True)
    matrix = csr_matrix(
        (np.ones(all_rows.size, dtype=np.float32), (all_rows, columns.astype(np.int64))),
        shape=(len(records), unique.size),
    )
    counts = np.asarray((matrix > 0).sum(axis=0)).ravel()
    keep = counts <= ubiquitous_threshold * len(records)
    return matrix[:, keep], int(unique.size), int(keep.sum())


def mash_distance(matrix: csr_matrix, k: int) -> np.ndarray:
    binary = (matrix > 0).astype(np.float32)
    intersection = (binary @ binary.T).toarray().astype(np.float32, copy=False)
    sizes = np.asarray(binary.sum(axis=1)).ravel().astype(np.float32)
    union = sizes[:, None] + sizes[None, :] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0)
    np.fill_diagonal(jaccard, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = -(1.0 / float(k)) * np.log((2.0 * jaccard) / (1.0 + jaccard))
    finite_positive = distance[np.isfinite(distance) & (distance > 0)]
    cap = float(finite_positive.max()) if finite_positive.size else 1.0
    distance[~np.isfinite(distance)] = cap
    distance[distance < 0] = 0.0
    np.fill_diagonal(distance, 0.0)
    return ((distance + distance.T) / 2.0).astype(np.float32, copy=False)


def silhouette_score(distance: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels)
    unique, inverse = np.unique(labels, return_inverse=True)
    if unique.size < 2 or unique.size >= len(labels):
        return float("nan")
    sizes = np.bincount(inverse)
    sums = np.zeros((len(labels), unique.size), dtype=np.float64)
    for cluster in range(unique.size):
        sums[:, cluster] = distance[:, inverse == cluster].sum(axis=1)
    means = sums / sizes[np.newaxis, :]
    within = np.zeros(len(labels), dtype=np.float64)
    for cluster in range(unique.size):
        members = inverse == cluster
        if sizes[cluster] > 1:
            within[members] = sums[members, cluster] / (sizes[cluster] - 1)
    means[np.arange(len(labels)), inverse] = np.inf
    nearest = means.min(axis=1)
    denominator = np.maximum(within, nearest)
    scores = np.zeros(len(labels), dtype=np.float64)
    valid = np.isfinite(nearest) & (denominator > 0)
    scores[valid] = (nearest[valid] - within[valid]) / denominator[valid]
    scores[sizes[inverse] <= 1] = 0.0
    return float(scores.mean())


def choose_cluster_count(scores: dict[int, float], tolerance: float) -> int:
    """Choose the simplest solution within tolerance of the best silhouette."""
    best_score = max(scores.values())
    return min(
        count for count, score in scores.items()
        if score >= best_score - tolerance
    )


def adjusted_rand_score(true_labels: list[str], labels: np.ndarray) -> float:
    if len(true_labels) < 2 or len(set(true_labels)) < 2:
        return float("nan")
    true_values = {value: i for i, value in enumerate(sorted(set(true_labels)))}
    pred_values = {int(value): i for i, value in enumerate(sorted(set(labels)))}
    table = np.zeros((len(true_values), len(pred_values)), dtype=np.float64)
    for true, pred in zip(true_labels, labels):
        table[true_values[true], pred_values[int(pred)]] += 1
    comb2 = lambda x: x * (x - 1.0) / 2.0
    sum_cells = comb2(table).sum()
    sum_rows = comb2(table.sum(axis=1)).sum()
    sum_cols = comb2(table.sum(axis=0)).sum()
    total = comb2(float(len(true_labels)))
    expected = sum_rows * sum_cols / total
    maximum = 0.5 * (sum_rows + sum_cols)
    if maximum == expected:
        return 1.0 if sum_cells == maximum else 0.0
    return float((sum_cells - expected) / (maximum - expected))


def pcoa_top2(distance: np.ndarray) -> tuple[np.ndarray | None, np.ndarray | None]:
    n = distance.shape[0]
    if n < 3:
        return None, None
    squared = distance.astype(np.float64) ** 2
    def matvec(vector: np.ndarray) -> np.ndarray:
        centered = vector - vector.mean()
        product = squared @ centered
        return -0.5 * (product - product.mean())
    try:
        values, vectors = eigsh(LinearOperator((n, n), matvec=matvec, dtype=np.float64), k=min(4, n - 1), which="LA")
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
    except Exception:
        centering = np.eye(n) - np.ones((n, n)) / n
        gram = -0.5 * centering.dot(squared).dot(centering)
        values, vectors = np.linalg.eigh((gram + gram.T) / 2.0)
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
    positive = values > 1e-9
    values, vectors = values[positive], vectors[:, positive]
    if values.size == 0:
        return None, None
    if values.size == 1:
        coordinates = np.column_stack(
            (vectors[:, 0] * np.sqrt(values[0]), np.zeros(n, dtype=float))
        )
        return coordinates, np.asarray([1.0, 0.0])
    return vectors[:, :2] * np.sqrt(values[:2]), values[:2] / values.sum()


PCOA_FIELDS = [
    "global_tag_id", "record_id", "sample_id", "gene_id", "gene_name",
    "copy_index", "sequence_group_id", "sequence_group_color", "PCo1", "PCo2",
    "PCo1_explained_fraction", "PCo2_explained_fraction",
]


def sequence_group_metadata(
    tag: str,
    records: list[tuple[str, str]],
    labels: np.ndarray,
) -> tuple[list[str], dict[str, str]]:
    """Convert scipy labels to reproducible public IDs and colors."""
    members: dict[int, list[str]] = {}
    for (record_id, _), label in zip(records, labels):
        members.setdefault(int(label), []).append(record_id)
    ordered = sorted(members, key=lambda label: min(members[label]))
    width = max(2, len(str(len(ordered))))
    label_to_group = {
        label: f"{tag}-SG{index:0{width}d}"
        for index, label in enumerate(ordered, start=1)
    }
    group_ids = [label_to_group[int(label)] for label in labels]
    tab20 = plt.get_cmap("tab20")
    extended = plt.get_cmap("hsv")
    colors = {}
    for index, label in enumerate(ordered):
        color = (
            tab20(index)
            if index < tab20.N
            else extended((index - tab20.N) / max(len(ordered) - tab20.N, 1))
        )
        colors[label_to_group[label]] = to_hex(color, keep_alpha=False)
    return group_ids, colors


def write_pcoa(
    path: str,
    coordinates_path: str,
    distance: np.ndarray,
    records: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
    labels: np.ndarray,
    genes: list[str],
    tag: str,
    args: argparse.Namespace,
) -> None:
    coordinates, explained = pcoa_top2(distance)
    if coordinates is None:
        return
    group_ids, group_colors = sequence_group_metadata(tag, records, labels)
    coordinate_rows = []
    for (record_id, _), group_id, point in zip(records, group_ids, coordinates):
        meta = metadata[record_id]
        coordinate_rows.append({
            "global_tag_id": tag,
            "record_id": record_id,
            "sample_id": meta["sample_id"],
            "gene_id": meta["gene_id"],
            "gene_name": meta["gene_name"],
            "copy_index": meta["copy_index"],
            "sequence_group_id": group_id,
            "sequence_group_color": group_colors[group_id],
            "PCo1": f"{point[0]:.8g}",
            "PCo2": f"{point[1]:.8g}",
            "PCo1_explained_fraction": f"{explained[0]:.8g}",
            "PCo2_explained_fraction": f"{explained[1]:.8g}",
        })
    write_tsv(coordinates_path, PCOA_FIELDS, coordinate_rows)
    if not path:
        return
    rng = np.random.default_rng(args.seed)
    plotted = coordinates.copy()
    if args.pcoa_jitter_frac > 0:
        scale = max(float(np.ptp(coordinates[:, 0])), float(np.ptp(coordinates[:, 1])), 1e-9)
        plotted += rng.normal(0.0, scale * args.pcoa_jitter_frac, plotted.shape)
    figure, axes = plt.subplots(1, 2, figsize=(16, 7), sharex=True, sharey=True)
    gene_names = sorted(set(genes))
    cmap = plt.get_cmap("tab20")

    for index, gene in enumerate(gene_names):
        selected = np.asarray(genes) == gene
        axes[0].scatter(
            plotted[selected, 0],
            plotted[selected, 1],
            s=args.pcoa_point_size,
            alpha=args.pcoa_alpha,
            color=cmap(index % cmap.N),
            label=gene,
            linewidths=0,
        )
    for group_id in sorted(group_colors):
        selected = np.asarray(group_ids) == group_id
        axes[1].scatter(
            plotted[selected, 0],
            plotted[selected, 1],
            s=args.pcoa_point_size,
            alpha=args.pcoa_alpha,
            color=group_colors[group_id],
            label=group_id,
            linewidths=0,
        )
    for axis in axes:
        axis.set_xlabel(f"PCo1 ({explained[0] * 100:.2f}%)")
        axis.legend(fontsize=7, markerscale=3, loc="best")
    axes[0].set_ylabel(f"PCo2 ({explained[1] * 100:.2f}%)")
    axes[0].set_title("Annotated gene names")
    axes[1].set_title("Sequence-defined groups")
    figure.suptitle(f"{tag}: Sequence Clustering Based on {args.k}-mers")
    figure.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=300)
    plt.close(figure)


def write_unavailable_pcoa(path: str, tag: str, status: str) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(16, 7), sharex=True, sharey=True)
    for axis, title in zip(axes, ("Annotated gene names", "Sequence-defined groups")):
        axis.axis("off")
        axis.set_title(title)
        axis.text(
            0.5, 0.5, f"PCoA unavailable: {status}",
            ha="center", va="center", transform=axis.transAxes,
        )
    figure.suptitle(tag)
    figure.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=300)
    plt.close(figure)


SUMMARY_FIELDS = [
    "global_tag_id", "status", "k", "ubiquitous_threshold", "n_sequences",
    "n_unique_sequences", "n_genes", "n_samples", "n_kmers_total",
    "n_kmers_kept", "best_n_clusters", "best_silhouette", "ARI_vs_gene",
]


def write_tsv(path: str, fields: list[str], rows: list[dict[str, object]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_one(args: argparse.Namespace) -> None:
    if args.k < 1:
        raise ClusterError("k must be a positive integer")
    if not 0 < args.ubiquitous_threshold <= 1:
        raise ClusterError("ubiquitous-threshold must be in (0, 1]")
    if args.threads < 1:
        raise ClusterError("threads must be a positive integer")
    if not 0 <= args.silhouette_tolerance <= 1:
        raise ClusterError("silhouette-tolerance must be between 0 and 1")
    records = read_fasta(args.fasta)
    metadata = read_metadata(args.records, args.global_tag_id)
    if any(header not in metadata for header, _ in records):
        raise ClusterError("FASTA contains record IDs absent from sequence_records.tsv")
    genes = [metadata[header]["gene_name"] for header, _ in records]
    samples = [metadata[header]["sample_id"] for header, _ in records]
    summary: dict[str, object] = {
        "global_tag_id": args.global_tag_id, "status": "ok", "k": args.k,
        "ubiquitous_threshold": args.ubiquitous_threshold,
        "n_sequences": len(records),
        "n_unique_sequences": len({hashlib.md5(seq.encode()).hexdigest() for _, seq in records}),
        "n_genes": len(set(genes)), "n_samples": len(set(samples)),
        "n_kmers_total": 0, "n_kmers_kept": 0, "best_n_clusters": "",
        "best_silhouette": "", "ARI_vs_gene": "",
    }
    write_tsv(args.pcoa_coordinates, PCOA_FIELDS, [])
    if len(records) < args.min_sequences:
        summary["status"] = "too_few_sequences"
    else:
        matrix, total_kmers, kept_kmers = build_kmer_matrix(
            records, args.k, args.ubiquitous_threshold, args.threads
        )
        summary["n_kmers_total"], summary["n_kmers_kept"] = total_kmers, kept_kmers
        if kept_kmers == 0 or np.any(np.asarray(matrix.sum(axis=1)).ravel() == 0):
            summary["status"] = "no_usable_kmers"
        else:
            distance = mash_distance(matrix, args.k)
            tree = linkage(squareform(distance, checks=False), method="average")
            automatic_max = len(set(genes)) + 3
            maximum = args.max_clusters if args.max_clusters > 0 else automatic_max
            maximum = min(maximum, len(records) - 1)
            scores: dict[int, float] = {1: 0.0}
            for count in range(2, maximum + 1):
                candidate = fcluster(tree, t=count, criterion="maxclust")
                score = silhouette_score(distance, candidate)
                if math.isfinite(score):
                    scores[count] = score
            best = choose_cluster_count(scores, args.silhouette_tolerance)
            labels = np.ones(len(records), dtype=int) if best == 1 else fcluster(tree, t=best, criterion="maxclust")
            summary["best_n_clusters"] = len(set(int(label) for label in labels))
            summary["best_silhouette"] = scores[best]
            summary["ARI_vs_gene"] = adjusted_rand_score(genes, labels)
            write_pcoa(
                args.pcoa, args.pcoa_coordinates, distance, records, metadata,
                labels, genes, args.global_tag_id, args,
            )
    for field in (
        "ubiquitous_threshold", "best_silhouette", "ARI_vs_gene",
    ):
        value = summary.get(field, "")
        if value != "" and isinstance(value, (int, float, np.floating)):
            numeric = float(value)
            summary[field] = f"{numeric:.2f}" if math.isfinite(numeric) else "NA"
    if args.pcoa and not Path(args.pcoa).is_file():
        reason = (
            str(summary["status"])
            if summary["status"] != "ok"
            else "insufficient PCoA dimensions"
        )
        write_unavailable_pcoa(args.pcoa, args.global_tag_id, reason)
    write_tsv(args.summary, SUMMARY_FIELDS, [summary])


def merge_summaries(args: argparse.Namespace) -> None:
    rows_by_tag: dict[str, dict[str, str]] = {}
    for path in args.inputs:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(SUMMARY_FIELDS) - set(reader.fieldnames or [])
            if missing:
                raise ClusterError(
                    f"{path} is missing summary column(s): {', '.join(sorted(missing))}"
                )
            for row in reader:
                projected = {field: row.get(field, "") for field in SUMMARY_FIELDS}
                if projected["global_tag_id"]:
                    rows_by_tag[projected["global_tag_id"]] = projected
    rows = sorted(rows_by_tag.values(), key=lambda row: row["global_tag_id"])
    write_tsv(args.output, SUMMARY_FIELDS, rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--global-tag-id", required=True)
    run.add_argument("--fasta", required=True)
    run.add_argument("--records", required=True)
    run.add_argument("--summary", required=True)
    run.add_argument("--pcoa", default="")
    run.add_argument("--pcoa-coordinates", required=True)
    run.add_argument("--k", type=int, required=True)
    run.add_argument("--ubiquitous-threshold", type=float, required=True)
    run.add_argument("--max-clusters", type=int, required=True, help="0 uses gene-name count + 3")
    run.add_argument("--silhouette-tolerance", type=float, required=True)
    run.add_argument("--min-sequences", type=int, required=True)
    run.add_argument("--seed", type=int, required=True)
    run.add_argument("--pcoa-jitter-frac", type=float, required=True)
    run.add_argument("--pcoa-point-size", type=float, required=True)
    run.add_argument("--pcoa-alpha", type=float, required=True)
    run.add_argument("--threads", type=int, required=True)
    merge = subparsers.add_parser("merge")
    merge.add_argument("--inputs", nargs="*", default=[])
    merge.add_argument("--output", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        run_one(args) if args.command == "run" else merge_summaries(args)
    except (OSError, ValueError, ClusterError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
