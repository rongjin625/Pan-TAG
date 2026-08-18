#!/usr/bin/env python3
"""Merge per-sample local TAGs into pangenome-wide TAGs.

The matching behavior is based on pipeline8_flank_based_merge_union.py, with
explicit inputs, dynamic storage, deterministic ordering and evidence tables.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


CLEAN_SUFFIX = re.compile(r"(_\d+|\(\d+\))+$")


class MergeError(RuntimeError):
    """A user-facing global TAG merge error."""


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def clean_values(value: str) -> set[str]:
    result: set[str] = set()
    for item in (value or "").split(","):
        item = item.strip()
        if not item or item.lower() in {"none", "boundary", "start", "end"}:
            continue
        result.add(CLEAN_SUFFIX.sub("", item))
    return result


def raw_values(value: str) -> set[str]:
    result: set[str] = set()
    for item in (value or "").split(","):
        cleaned = item.strip()
        if cleaned and cleaned.lower() not in {"none", "boundary", "start", "end"}:
            result.add(cleaned)
    return result


def sample_from_path(path: str | Path) -> str:
    file_path = Path(path)
    for suffix in ("_TAGs.tsv",):
        if file_path.name.endswith(suffix):
            return file_path.name.removesuffix(suffix)
    if file_path.parent.name:
        return file_path.parent.name
    return file_path.stem


@dataclass
class LocalTag:
    sample: str
    cluster_id: str
    chrom: str
    member_names: set[str]
    member_ids: set[str]
    pseudo_names: set[str]
    pseudo_ids: set[str]
    left_flanks: set[str]
    right_flanks: set[str]
    names: set[str] = field(init=False)
    members: set[str] = field(init=False)
    flanks: set[str] = field(init=False)

    def __post_init__(self) -> None:
        self.names = self.member_names | self.pseudo_names
        self.members = self.member_ids | self.pseudo_ids
        self.flanks = self.left_flanks | self.right_flanks


@dataclass
class GlobalTag:
    chrom: str
    source: str
    names: set[str]
    members: set[str]
    left_flanks: set[str]
    right_flanks: set[str]
    local_tags: list[tuple[str, str]] = field(default_factory=list)
    score_sum: float = 1.0
    match_count: int = 1
    min_score: float = 1.0
    active: bool = True
    flanks: set[str] = field(init=False)

    def __post_init__(self) -> None:
        self.flanks = self.left_flanks | self.right_flanks


def load_local_tags(path: str | Path, sample: str) -> list[LocalTag]:
    required = {
        "Cluster_ID",
        "Chromosome",
        "Upstream_Genes",
        "Downstream_Genes",
        "Members_Annotated",
        "Members_IDs",
    }
    result: list[LocalTag] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise MergeError(
                f"{path} is missing columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            result.append(
                LocalTag(
                    sample=sample,
                    cluster_id=row["Cluster_ID"],
                    chrom=row["Chromosome"],
                    member_names=clean_values(row["Members_Annotated"]),
                    member_ids=raw_values(row["Members_IDs"]),
                    pseudo_names=clean_values(row.get("Pseudos_Annotated", "")),
                    pseudo_ids=raw_values(row.get("Pseudos_IDs", "")),
                    left_flanks=clean_values(row["Upstream_Genes"]),
                    right_flanks=clean_values(row["Downstream_Genes"]),
                )
            )
    return result


def overlap(left: set[str], right: set[str]) -> tuple[float, float, float]:
    intersection = len(left & right)
    forward = intersection / len(left) if left else 0.0
    reverse = intersection / len(right) if right else 0.0
    return forward, reverse, max(forward, reverse)


def flank_jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def index_values(
    index: dict[str, set[int]], values: set[str], global_index: int
) -> None:
    for value in values:
        index[value].add(global_index)


def merge_global_nodes(
    target: GlobalTag,
    source: GlobalTag,
    target_index: int | None = None,
    name_index: dict[str, set[int]] | None = None,
    flank_index: dict[str, set[int]] | None = None,
) -> None:
    new_names = source.names - target.names
    new_flanks = source.flanks - target.flanks
    target.names |= source.names
    target.members |= source.members
    target.left_flanks |= source.left_flanks
    target.right_flanks |= source.right_flanks
    target.flanks |= source.flanks
    target.local_tags.extend(source.local_tags)
    target.score_sum += source.score_sum
    target.match_count += source.match_count
    target.min_score = min(target.min_score, source.min_score)
    source.active = False
    if target_index is not None and name_index is not None and flank_index is not None:
        index_values(name_index, new_names, target_index)
        index_values(flank_index, new_flanks, target_index)


def add_local_tag(
    global_tag: GlobalTag,
    local_tag: LocalTag,
    score: float,
    global_index: int | None = None,
    name_index: dict[str, set[int]] | None = None,
    flank_index: dict[str, set[int]] | None = None,
) -> None:
    new_names = local_tag.names - global_tag.names
    new_flanks = local_tag.flanks - global_tag.flanks
    global_tag.names |= local_tag.names
    global_tag.members |= local_tag.members
    global_tag.left_flanks |= local_tag.left_flanks
    global_tag.right_flanks |= local_tag.right_flanks
    global_tag.flanks |= local_tag.flanks
    global_tag.local_tags.append((local_tag.sample, local_tag.cluster_id))
    global_tag.score_sum += score
    global_tag.match_count += 1
    global_tag.min_score = min(global_tag.min_score, score)
    if global_index is not None and name_index is not None and flank_index is not None:
        index_values(name_index, new_names, global_index)
        index_values(flank_index, new_flanks, global_index)


def indexed_candidates(
    local_tag: LocalTag,
    name_index: dict[str, set[int]],
    flank_index: dict[str, set[int]],
) -> list[int]:
    name_candidates: set[int] = set()
    for name in local_tag.names:
        name_candidates.update(name_index.get(name, ()))
    flank_candidates: set[int] = set()
    for flank in local_tag.flanks:
        flank_candidates.update(flank_index.get(flank, ()))
    return sorted(name_candidates & flank_candidates)


def candidate_indices(
    local_tag: LocalTag,
    global_tags: list[GlobalTag],
    threshold: float,
    require_same_seqid: bool = False,
    indices: list[int] | None = None,
) -> list[dict[str, float | int]]:
    candidates: list[dict[str, float | int]] = []
    search_indices = indices if indices is not None else range(len(global_tags))
    for index in search_indices:
        global_tag = global_tags[index]
        if not global_tag.active:
            continue
        if require_same_seqid and local_tag.chrom != global_tag.chrom:
            continue
        if indices is None and local_tag.flanks.isdisjoint(global_tag.flanks):
            continue
        intersection = len(local_tag.names & global_tag.names)
        if not intersection:
            continue
        forward = intersection / len(local_tag.names) if local_tag.names else 0.0
        reverse = intersection / len(global_tag.names) if global_tag.names else 0.0
        score = max(forward, reverse)
        if score < threshold:
            continue
        candidates.append(
            {
                "index": index,
                "forward": forward,
                "reverse": reverse,
                "score": score,
                "flank_jaccard": flank_jaccard(local_tag.flanks, global_tag.flanks),
            }
        )
    return candidates


def write_tsv(
    path: str | Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> None:
    output = ensure_parent(path)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)


def command_run(args: argparse.Namespace) -> int:
    require_same_seqid = bool(getattr(args, "require_same_seqid", False))
    print(f"[Pan-TAG] require_same_seqid={str(require_same_seqid).lower()}")
    paths_by_sample: dict[str, str] = {}
    for path in args.local_tags:
        sample = sample_from_path(path)
        if sample in paths_by_sample:
            raise MergeError(f"Duplicate local TAG input for sample: {sample}")
        paths_by_sample[sample] = path
    if args.reference_sample not in paths_by_sample:
        raise MergeError(
            f"Reference sample {args.reference_sample!r} is missing from --local-tags"
        )

    samples = [args.reference_sample] + sorted(
        sample for sample in paths_by_sample if sample != args.reference_sample
    )
    globals_: list[GlobalTag] = []
    name_index: dict[str, set[int]] = defaultdict(set)
    flank_index: dict[str, set[int]] = defaultdict(set)
    total_local_tags = 0

    reference_tags = load_local_tags(
        paths_by_sample[args.reference_sample], args.reference_sample
    )
    total_local_tags += len(reference_tags)
    for local_tag in reference_tags:
        global_index = len(globals_)
        globals_.append(
            GlobalTag(
                chrom=local_tag.chrom,
                source=args.reference_sample,
                names=set(local_tag.names),
                members=set(local_tag.members),
                left_flanks=set(local_tag.left_flanks),
                right_flanks=set(local_tag.right_flanks),
                local_tags=[(local_tag.sample, local_tag.cluster_id)],
            )
        )
        index_values(name_index, local_tag.names, global_index)
        index_values(flank_index, local_tag.flanks, global_index)
    del reference_tags

    for sample in samples[1:]:
        sample_tags = load_local_tags(paths_by_sample[sample], sample)
        total_local_tags += len(sample_tags)
        for local_tag in sample_tags:
            search_indices = indexed_candidates(local_tag, name_index, flank_index)
            candidates = candidate_indices(
                local_tag,
                globals_,
                args.member_overlap_threshold,
                require_same_seqid,
                search_indices,
            )
            target_index: int | None = None
            score = 1.0

            if len(candidates) == 1:
                target_index = int(candidates[0]["index"])
                score = float(candidates[0]["score"])
            elif len(candidates) > 1:
                union_names = set().union(
                    *(globals_[int(candidate["index"])].names for candidate in candidates)
                )
                union_forward, union_reverse, union_score = overlap(
                    local_tag.names, union_names
                )
                best_forward = max(float(item["forward"]) for item in candidates)
                best_reverse = max(float(item["reverse"]) for item in candidates)
                if union_forward > best_forward or union_reverse > best_reverse:
                    ordered_indices = sorted(int(item["index"]) for item in candidates)
                    target_index = ordered_indices[0]
                    for source_index in ordered_indices[1:]:
                        if not globals_[source_index].active:
                            continue
                        merge_global_nodes(
                            globals_[target_index],
                            globals_[source_index],
                            target_index,
                            name_index,
                            flank_index,
                        )
                    score = union_score
                else:
                    best = sorted(
                        candidates,
                        key=lambda item: (
                            -float(item["score"]),
                            -float(item["flank_jaccard"]),
                            int(item["index"]),
                        ),
                    )[0]
                    target_index = int(best["index"])
                    score = float(best["score"])

            if target_index is None:
                target_index = len(globals_)
                globals_.append(
                    GlobalTag(
                        chrom=local_tag.chrom,
                        source=sample,
                        names=set(local_tag.names),
                        members=set(local_tag.members),
                        left_flanks=set(local_tag.left_flanks),
                        right_flanks=set(local_tag.right_flanks),
                        local_tags=[(local_tag.sample, local_tag.cluster_id)],
                    )
                )
                index_values(name_index, local_tag.names, target_index)
                index_values(flank_index, local_tag.flanks, target_index)
            else:
                add_local_tag(
                    globals_[target_index],
                    local_tag,
                    score,
                    target_index,
                    name_index,
                    flank_index,
                )
        del sample_tags

    active = [
        (index, global_tag)
        for index, global_tag in enumerate(globals_)
        if global_tag.active
    ]
    active.sort(
        key=lambda item: (
            0 if item[1].source == args.reference_sample else 1,
            item[1].chrom,
            min(
                (
                    local
                    for local in item[1].local_tags
                ),
                default=("", ""),
            ),
            item[0],
        )
    )
    internal_to_global = {
        internal: f"Candidate_cluster_{number:04d}"
        for number, (internal, _) in enumerate(active, start=1)
    }

    tag_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for internal_index, global_tag in active:
        global_id = internal_to_global[internal_index]
        sample_map: dict[str, list[str]] = defaultdict(list)
        for sample, cluster_id in global_tag.local_tags:
            sample_map[sample].append(cluster_id)
            mapping_rows.append(
                {
                    "candidate_global_tag_id": global_id,
                    "sample_id": sample,
                    "local_tag_id": cluster_id,
                }
            )
        tag_rows.append(
            {
                "candidate_global_tag_id": global_id,
                "pan_cluster_size": len(global_tag.names),
                "reference_seqid": global_tag.chrom,
                "source_sample_id": global_tag.source,
                "occupancy": len(sample_map),
                "gene_names": ",".join(sorted(global_tag.names)),
                "gene_ids": ",".join(sorted(global_tag.members)),
                "left_flank_genes": ",".join(sorted(global_tag.left_flanks)),
                "right_flank_genes": ",".join(sorted(global_tag.right_flanks)),
            }
        )

    write_tsv(
        args.output_tags,
        [
            "candidate_global_tag_id",
            "pan_cluster_size",
            "reference_seqid",
            "source_sample_id",
            "occupancy",
            "gene_names",
            "gene_ids",
            "left_flank_genes",
            "right_flank_genes",
        ],
        tag_rows,
    )
    write_tsv(
        args.output_mapping,
        ["candidate_global_tag_id", "sample_id", "local_tag_id"],
        sorted(
            mapping_rows,
            key=lambda row: (
                row["candidate_global_tag_id"],
                row["sample_id"],
                row["local_tag_id"],
            ),
        ),
    )
    print(
        f"[Pan-TAG] Merged {total_local_tags} "
        f"local TAG(s) into {len(active)} global TAG(s)."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--local-tags", nargs="+", required=True)
    parser.add_argument("--member-overlap-threshold", type=float, default=0.1)
    parser.add_argument(
        "--require-same-seqid",
        action="store_true",
        help="Only merge local/global TAGs with identical seq_id values.",
    )
    parser.add_argument("--output-tags", required=True)
    parser.add_argument("--output-mapping", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return command_run(args)
    except (MergeError, OSError, ValueError, KeyError) as exc:
        print(f"[Pan-TAG ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
