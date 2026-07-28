#!/usr/bin/env python3
"""Merge per-sample local TAGs into pangenome-wide TAGs.

The matching behavior is based on pipeline8_flank_based_merge_union.py, with
explicit inputs, dynamic storage, deterministic ordering and evidence tables.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


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
        result.add(re.sub(r"(_\d+|\(\d+\))+$", "", item))
    return result


def raw_values(value: str) -> set[str]:
    return {
        item.strip()
        for item in (value or "").split(",")
        if item.strip()
        and item.strip().lower() not in {"none", "boundary", "start", "end"}
    }


def sample_from_path(path: str | Path) -> str:
    file_path = Path(path)
    if file_path.parent.name:
        return file_path.parent.name
    suffix = "_with_10_flanking_cleaned.txt"
    return file_path.name.removesuffix(suffix)


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

    @property
    def names(self) -> set[str]:
        return self.member_names | self.pseudo_names

    @property
    def members(self) -> set[str]:
        return self.member_ids | self.pseudo_ids

    @property
    def flanks(self) -> set[str]:
        return self.left_flanks | self.right_flanks


@dataclass
class GlobalTag:
    chrom: str
    source: str
    names: set[str]
    members: set[str]
    left_flanks: set[str]
    right_flanks: set[str]
    local_tags: list[LocalTag] = field(default_factory=list)
    score_sum: float = 1.0
    match_count: int = 1
    min_score: float = 1.0
    active: bool = True

    @property
    def flanks(self) -> set[str]:
        return self.left_flanks | self.right_flanks


def load_local_tags(path: str | Path, sample: str) -> list[LocalTag]:
    required = {
        "Cluster_ID",
        "Chrom",
        "Up_5th_ID",
        "Down_5th_ID",
        "Members_Annotated",
        "Members_IDs",
        "Pseudos_Annotated",
        "Pseudos_IDs",
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
                    chrom=row["Chrom"],
                    member_names=clean_values(row["Members_Annotated"]),
                    member_ids=raw_values(row["Members_IDs"]),
                    pseudo_names=clean_values(row["Pseudos_Annotated"]),
                    pseudo_ids=raw_values(row["Pseudos_IDs"]),
                    left_flanks=clean_values(row["Up_5th_ID"]),
                    right_flanks=clean_values(row["Down_5th_ID"]),
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


def merge_global_nodes(target: GlobalTag, source: GlobalTag) -> None:
    target.names |= source.names
    target.members |= source.members
    target.left_flanks |= source.left_flanks
    target.right_flanks |= source.right_flanks
    target.local_tags.extend(source.local_tags)
    target.score_sum += source.score_sum
    target.match_count += source.match_count
    target.min_score = min(target.min_score, source.min_score)
    source.active = False


def add_local_tag(global_tag: GlobalTag, local_tag: LocalTag, score: float) -> None:
    global_tag.names |= local_tag.names
    global_tag.members |= local_tag.members
    global_tag.left_flanks |= local_tag.left_flanks
    global_tag.right_flanks |= local_tag.right_flanks
    global_tag.local_tags.append(local_tag)
    global_tag.score_sum += score
    global_tag.match_count += 1
    global_tag.min_score = min(global_tag.min_score, score)


def candidate_indices(
    local_tag: LocalTag, global_tags: list[GlobalTag], threshold: float
) -> list[dict[str, float | int]]:
    candidates: list[dict[str, float | int]] = []
    for index, global_tag in enumerate(global_tags):
        if not global_tag.active:
            continue
        if local_tag.flanks.isdisjoint(global_tag.flanks):
            continue
        if local_tag.names.isdisjoint(global_tag.names):
            continue
        forward, reverse, score = overlap(local_tag.names, global_tag.names)
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
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def command_run(args: argparse.Namespace) -> int:
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
    local_by_sample = {
        sample: load_local_tags(paths_by_sample[sample], sample) for sample in samples
    }

    globals_: list[GlobalTag] = []
    evidence_rows: list[dict[str, object]] = []
    fusion_rows: list[dict[str, object]] = []
    redirects: dict[int, int] = {}

    for local_tag in local_by_sample[args.reference_sample]:
        globals_.append(
            GlobalTag(
                chrom=local_tag.chrom,
                source=args.reference_sample,
                names=set(local_tag.names),
                members=set(local_tag.members),
                left_flanks=set(local_tag.left_flanks),
                right_flanks=set(local_tag.right_flanks),
                local_tags=[local_tag],
            )
        )

    for sample in samples[1:]:
        for local_tag in local_by_sample[sample]:
            candidates = candidate_indices(
                local_tag, globals_, args.member_overlap_threshold
            )
            target_index: int | None = None
            score = 1.0
            decision = "NEW"

            if len(candidates) == 1:
                target_index = int(candidates[0]["index"])
                score = float(candidates[0]["score"])
                decision = "MATCH"
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
                        fusion_rows.append(
                            {
                                "sample_id": sample,
                                "local_tag_id": local_tag.cluster_id,
                                "target_internal_id": target_index,
                                "source_internal_id": source_index,
                                "union_score": round(union_score, 6),
                            }
                        )
                        merge_global_nodes(globals_[target_index], globals_[source_index])
                        redirects[source_index] = target_index
                    score = union_score
                    decision = "FUSION"
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
                    decision = "BEST_MATCH"

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
                        local_tags=[local_tag],
                    )
                )
            else:
                add_local_tag(globals_[target_index], local_tag, score)

            evidence_rows.append(
                {
                    "sample_id": sample,
                    "local_tag_id": local_tag.cluster_id,
                    "target_internal_id": target_index,
                    "candidate_count": len(candidates),
                    "match_score": round(score, 6),
                    "decision": decision,
                }
            )

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
                    (local.sample, local.cluster_id)
                    for local in item[1].local_tags
                ),
                default=("", ""),
            ),
            item[0],
        )
    )
    internal_to_global = {
        internal: f"MC_{number:04d}"
        for number, (internal, _) in enumerate(active, start=1)
    }

    def resolve_internal(index: int) -> int:
        seen: set[int] = set()
        while index in redirects:
            if index in seen:
                raise MergeError("Internal fusion redirect cycle detected")
            seen.add(index)
            index = redirects[index]
        return index

    tag_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for internal_index, global_tag in active:
        global_id = internal_to_global[internal_index]
        sample_map: dict[str, list[str]] = defaultdict(list)
        for local in global_tag.local_tags:
            sample_map[local.sample].append(local.cluster_id)
            mapping_rows.append(
                {
                    "global_tag_id": global_id,
                    "sample_id": local.sample,
                    "local_tag_id": local.cluster_id,
                }
            )
            for member_id in sorted(local.member_ids):
                member_rows.append(
                    {
                        "global_tag_id": global_id,
                        "sample_id": local.sample,
                        "local_tag_id": local.cluster_id,
                        "member_type": "member",
                        "member_id": member_id,
                    }
                )
            for member_id in sorted(local.pseudo_ids):
                member_rows.append(
                    {
                        "global_tag_id": global_id,
                        "sample_id": local.sample,
                        "local_tag_id": local.cluster_id,
                        "member_type": "pseudogene",
                        "member_id": member_id,
                    }
                )
        tag_rows.append(
            {
                "global_tag_id": global_id,
                "reference_seqid": global_tag.chrom,
                "source_sample_id": global_tag.source,
                "occupancy": len(sample_map),
                "local_tag_count": len(global_tag.local_tags),
                "sample_ids": ",".join(sorted(sample_map)),
                "local_tag_ids": ";".join(
                    f"{sample}:{','.join(sorted(ids))}"
                    for sample, ids in sorted(sample_map.items())
                ),
                "gene_names": ",".join(sorted(global_tag.names)),
                "member_ids": ",".join(sorted(global_tag.members)),
                "member_instance_count": sum(
                    len(local.member_ids) for local in global_tag.local_tags
                ),
                "pseudogene_instance_count": sum(
                    len(local.pseudo_ids) for local in global_tag.local_tags
                ),
                "left_flank_genes": ",".join(sorted(global_tag.left_flanks)),
                "right_flank_genes": ",".join(sorted(global_tag.right_flanks)),
                "mean_match_score": round(
                    global_tag.score_sum / global_tag.match_count, 6
                ),
                "min_match_score": round(global_tag.min_score, 6),
            }
        )

    for row in evidence_rows:
        internal = resolve_internal(int(row.pop("target_internal_id")))
        row["global_tag_id"] = internal_to_global.get(internal, "")
    for row in fusion_rows:
        target = resolve_internal(int(row.pop("target_internal_id")))
        source = resolve_internal(int(row.pop("source_internal_id")))
        row["target_global_tag_id"] = internal_to_global.get(target, "")
        row["source_global_tag_id"] = internal_to_global.get(source, "")

    write_tsv(
        args.output_tags,
        [
            "global_tag_id",
            "reference_seqid",
            "source_sample_id",
            "occupancy",
            "local_tag_count",
            "sample_ids",
            "local_tag_ids",
            "gene_names",
            "member_ids",
            "member_instance_count",
            "pseudogene_instance_count",
            "left_flank_genes",
            "right_flank_genes",
            "mean_match_score",
            "min_match_score",
        ],
        tag_rows,
    )
    write_tsv(
        args.output_members,
        [
            "global_tag_id",
            "sample_id",
            "local_tag_id",
            "member_type",
            "member_id",
        ],
        sorted(
            member_rows,
            key=lambda row: (
                row["global_tag_id"],
                row["sample_id"],
                row["local_tag_id"],
                row["member_type"],
                row["member_id"],
            ),
        ),
    )
    write_tsv(
        args.output_mapping,
        ["global_tag_id", "sample_id", "local_tag_id"],
        sorted(
            mapping_rows,
            key=lambda row: (
                row["global_tag_id"],
                row["sample_id"],
                row["local_tag_id"],
            ),
        ),
    )
    write_tsv(
        args.output_evidence,
        [
            "sample_id",
            "local_tag_id",
            "candidate_count",
            "match_score",
            "decision",
            "global_tag_id",
        ],
        evidence_rows,
    )
    write_tsv(
        args.output_fusions,
        [
            "sample_id",
            "local_tag_id",
            "union_score",
            "target_global_tag_id",
            "source_global_tag_id",
        ],
        fusion_rows,
    )
    print(
        f"[Pan-TAG] Merged {sum(len(tags) for tags in local_by_sample.values())} "
        f"local TAG(s) into {len(active)} global TAG(s)."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--local-tags", nargs="+", required=True)
    parser.add_argument("--member-overlap-threshold", type=float, default=0.1)
    parser.add_argument("--output-tags", required=True)
    parser.add_argument("--output-members", required=True)
    parser.add_argument("--output-mapping", required=True)
    parser.add_argument("--output-evidence", required=True)
    parser.add_argument("--output-fusions", required=True)
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
