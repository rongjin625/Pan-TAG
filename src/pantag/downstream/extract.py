#!/usr/bin/env python3
"""Extract full gene intervals directly from a TAG sequence manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import time
from collections import Counter, OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pysam
from Bio.Seq import Seq


class ExtractError(RuntimeError):
    """User-facing extraction failure."""


def allow_large_csv_fields() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def sample_assets(samples, reference_sample, reference_fasta, fasta_indexes):
    rows = []
    with Path(samples).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "fasta"}.issubset(reader.fieldnames or []):
            raise ExtractError(f"{samples} must contain sample_id and fasta")
        for row in reader:
            sample = (row.get("sample_id") or "").strip()
            fasta = (row.get("fasta") or "").strip()
            if not sample or not fasta:
                raise ExtractError(f"{samples} contains an empty sample_id or fasta")
            rows.append((sample, reference_fasta if sample == reference_sample else fasta))
    if len(rows) != len(fasta_indexes):
        raise ExtractError(
            "--fasta-indexes must provide one index per samples.tsv row "
            f"({len(fasta_indexes)} indexes for {len(rows)} samples)"
        )
    if len({sample for sample, _ in rows}) != len(rows):
        raise ExtractError(f"{samples} contains duplicate sample_id values")
    return {
        sample: (fasta, index)
        for (sample, fasta), index in zip(rows, fasta_indexes)
    }


class WriterCache:
    def __init__(self, directory: Path, max_open: int = 128):
        self.directory = directory
        self.max_open = max_open
        self.handles: OrderedDict[str, object] = OrderedDict()

    def write(self, tag: str, header: str, sequence: str) -> None:
        if tag in self.handles:
            handle = self.handles.pop(tag)
        else:
            if len(self.handles) >= self.max_open:
                _, old = self.handles.popitem(last=False)
                old.close()
            handle = (self.directory / f"{tag}.fa").open("a", encoding="utf-8")
        self.handles[tag] = handle
        handle.write(f">{header}\n")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


SPOOL_FIELDS = [
    "record_id", "global_tag_id", "sample_id", "gene_id", "gene_name",
    "copy_index", "sequence_length", "sequence",
]
RECORD_FIELDS = SPOOL_FIELDS[:-1]
STATE_VERSION = 1


def extract_sample(sample, rows, fasta_path, index_path, spool_path):
    """Random-access one sample FASTA and atomically publish its private spool."""
    started = time.monotonic()
    temporary = spool_path.with_suffix(spool_path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        fasta = pysam.FastaFile(fasta_path, filepath_index=index_path)
    except (OSError, ValueError) as exc:
        raise ExtractError(
            f"{sample}: cannot open FASTA {fasta_path} with index {index_path}: {exc}"
        ) from exc
    try:
        references = set(fasta.references)
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=SPOOL_FIELDS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                seqid = row["seqid"]
                if seqid not in references:
                    raise ExtractError(
                        f"{sample}: seqid {seqid} is absent from FASTA {fasta_path}"
                    )
                start, end = int(row["start"]), int(row["end"])
                if start < 1 or end < start or end > fasta.get_reference_length(seqid):
                    raise ExtractError(
                        f"{sample}: invalid interval for {row['record_id']}: "
                        f"{seqid}:{start}-{end} in {fasta_path}"
                    )
                sequence = fasta.fetch(seqid, start - 1, end).upper()
                if row["strand"] == "-":
                    sequence = str(Seq(sequence).reverse_complement())
                writer.writerow({
                    "record_id": row["record_id"], "global_tag_id": row["global_tag_id"],
                    "sample_id": sample, "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"], "copy_index": row["copy_index"],
                    "sequence_length": len(sequence), "sequence": sequence,
                })
        os.replace(temporary, spool_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        fasta.close()
    return spool_path, len(rows), time.monotonic() - started


def content_fingerprint(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path):
    candidate = Path(path).resolve()
    stat = candidate.stat()
    return {"path": str(candidate), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def extraction_state(args, selected_tags, assets, rows_by_sample):
    return {
        "version": STATE_VERSION,
        "tags": selected_tags,
        "manifest_sha256": content_fingerprint(args.manifest),
        "matrix_sha256": content_fingerprint(args.matrix),
        "samples": [
            {
                "sample_id": sample,
                "expected_records": len(rows_by_sample[sample]),
                "fasta": file_identity(assets[sample][0]),
                "index": file_identity(assets[sample][1]),
            }
            for sample in rows_by_sample
        ],
    }


def read_state(path):
    try:
        with path.open(encoding="utf-8") as handle:
            state = json.load(handle)
    except (OSError, ValueError, TypeError):
        return None
    return state if isinstance(state, dict) else None


def write_state(path, state):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(state, handle, sort_keys=True, indent=2)
        handle.write("\n")
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> None:
    started = time.monotonic()
    if args.threads < 1:
        raise ExtractError("threads must be a positive integer")
    allow_large_csv_fields()
    assets = sample_assets(
        args.samples, args.reference_sample, args.reference_fasta, args.fasta_indexes
    )
    with Path(args.matrix).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "global_tag_id" not in (reader.fieldnames or []):
            raise ExtractError(f"{args.matrix} must contain global_tag_id")
        all_tags = [row["global_tag_id"].strip() for row in reader if row["global_tag_id"].strip()]
    tag_set = set(all_tags)
    if len(all_tags) != len(tag_set):
        raise ExtractError("global_tag_id values in the refined matrix must be unique")
    requested = list(dict.fromkeys(getattr(args, "tags", None) or []))
    missing_tags = sorted(set(requested) - tag_set)
    if missing_tags:
        raise ExtractError(
            "Requested Pan-cluster ID(s) are absent from the refined matrix: "
            + ", ".join(missing_tags)
        )
    selected_tags = all_tags if not requested else [tag for tag in all_tags if tag in requested]
    selected_set = set(selected_tags)

    required = {
        "record_id", "global_tag_id", "sample_id", "gene_id", "gene_name",
        "seqid", "start", "end", "strand", "copy_index", "qc_status",
        "extraction_eligible",
    }
    with Path(args.manifest).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ExtractError(f"Manifest missing columns: {', '.join(sorted(missing))}")
        manifest_rows = list(reader)
    unknown = sorted({row["global_tag_id"] for row in manifest_rows} - tag_set)
    if unknown:
        raise ExtractError(
            "Manifest contains TAGs absent from the final matrix: " + ", ".join(unknown[:10])
        )
    rows = [
        row for row in manifest_rows
        if row["global_tag_id"] in selected_set
        and row["extraction_eligible"].strip().lower() == "true"
    ]
    rows_by_sample = {}
    for row in rows:
        sample = row["sample_id"]
        if sample not in assets:
            raise ExtractError(f"No FASTA and index configured for sample {sample}")
        rows_by_sample.setdefault(sample, []).append(row)
    for sample in rows_by_sample:
        fasta_path, index_path = assets[sample]
        for path, label in ((fasta_path, "FASTA"), (index_path, "FASTA index")):
            if not Path(path).is_file():
                raise ExtractError(f"Missing {label} for {sample}: {path}")

    outdir = Path(args.output_dir)
    temporary = outdir.with_name(outdir.name + ".tmp")
    state = extraction_state(args, selected_tags, assets, rows_by_sample)
    state_path = temporary / ".extract_state.json"
    if temporary.exists() and read_state(state_path) != state:
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True, exist_ok=True)
    spool_directory = temporary / ".sample_spools"
    spool_directory.mkdir(exist_ok=True)
    if not state_path.is_file():
        write_state(state_path, state)

    spool_paths = {}
    work_items = []
    reused = 0
    for ordinal, sample in enumerate(rows_by_sample):
        fasta_path, index_path = assets[sample]
        spool_path = spool_directory / f"{ordinal:06d}.tsv"
        spool_path.with_suffix(spool_path.suffix + ".tmp").unlink(missing_ok=True)
        if spool_path.is_file():
            spool_paths[sample] = spool_path
            reused += 1
        else:
            work_items.append((sample, rows_by_sample[sample], fasta_path, index_path, spool_path))

    total_samples = len(rows_by_sample)
    worker_count = min(args.threads, max(1, len(work_items)))
    print(
        f"[extract] tags={len(selected_tags)} eligible_records={len(rows)} "
        f"samples={total_samples} reused={reused} workers={worker_count}", flush=True,
    )
    completed = reused
    if work_items:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {executor.submit(extract_sample, *item): item[0] for item in work_items}
            for future in as_completed(futures):
                sample = futures[future]
                spool_path, count, elapsed = future.result()
                spool_paths[sample] = spool_path
                completed += 1
                print(
                    f"[extract] completed={completed}/{total_samples} sample={sample} "
                    f"records={count} elapsed_sec={elapsed:.2f}", flush=True,
                )

    print("[extract] consolidating sample spools into per-TAG FASTA files", flush=True)
    for path in temporary.glob("*.fa"):
        path.unlink()
    for tag in selected_tags:
        (temporary / f"{tag}.fa").touch()
    writers = WriterCache(temporary)
    expected = Counter((row["global_tag_id"], row["sample_id"]) for row in rows)
    observed = Counter()
    records_path = Path(args.output_records)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    temp_records = records_path.with_suffix(records_path.suffix + ".tmp")
    try:
        with temp_records.open("w", encoding="utf-8", newline="") as records_handle:
            records_writer = csv.DictWriter(
                records_handle, fieldnames=RECORD_FIELDS, delimiter="\t", lineterminator="\n"
            )
            records_writer.writeheader()
            for sample in rows_by_sample:
                with spool_paths[sample].open(encoding="utf-8", newline="") as spool_handle:
                    for row in csv.DictReader(spool_handle, delimiter="\t"):
                        writers.write(row["global_tag_id"], row["record_id"], row["sequence"])
                        observed[(row["global_tag_id"], row["sample_id"])] += 1
                        records_writer.writerow({field: row[field] for field in RECORD_FIELDS})
    finally:
        writers.close()
    if expected != observed:
        raise ExtractError("Eligible manifest counts do not equal extracted FASTA counts")

    shutil.rmtree(spool_directory)
    state_path.unlink(missing_ok=True)
    if outdir.exists():
        shutil.rmtree(outdir)
    os.replace(temporary, outdir)
    os.replace(temp_records, records_path)
    print(
        f"[extract] DONE tags={len(selected_tags)} sequences={len(rows)} "
        f"elapsed_sec={time.monotonic() - started:.2f}", flush=True,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--fasta-indexes", nargs="+", required=True)
    parser.add_argument("--tag", action="append", dest="tags", metavar="GLOBAL_TAG_ID")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-records", required=True)
    return parser


def main() -> int:
    try:
        run(build_parser().parse_args())
    except (OSError, ValueError, csv.Error, ExtractError, pysam.utils.SamtoolsError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
