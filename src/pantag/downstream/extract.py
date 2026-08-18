#!/usr/bin/env python3
"""Extract full gene intervals directly from a TAG sequence manifest."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from collections import Counter, OrderedDict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


class ExtractError(RuntimeError):
    """User-facing extraction failure."""


def sample_fastas(samples: str, reference_sample: str, reference_fasta: str) -> dict[str, str]:
    result = {reference_sample: reference_fasta}
    with Path(samples).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "fasta"}.issubset(reader.fieldnames or []):
            raise ExtractError(f"{samples} must contain sample_id and fasta")
        for row in reader:
            result[row["sample_id"].strip()] = row["fasta"].strip()
    return result


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
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i : i + 80] + "\n")

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def run(args: argparse.Namespace) -> None:
    fasta_map = sample_fastas(args.samples, args.reference_sample, args.reference_fasta)
    with Path(args.matrix).open(encoding="utf-8", newline="") as handle:
        matrix_reader = csv.DictReader(handle, delimiter="\t")
        if "global_tag_id" not in (matrix_reader.fieldnames or []):
            raise ExtractError(f"{args.matrix} must contain global_tag_id")
        all_tags = [
            row["global_tag_id"].strip()
            for row in matrix_reader
            if row["global_tag_id"].strip()
        ]
    tag_set = set(all_tags)
    if len(all_tags) != len(tag_set):
        raise ExtractError("global_tag_id values in the refined matrix must be unique")
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
    unknown_tags = sorted(
        {row["global_tag_id"] for row in manifest_rows}.difference(tag_set)
    )
    if unknown_tags:
        raise ExtractError(
            "Manifest contains TAGs absent from the final matrix: "
            + ", ".join(unknown_tags[:10])
        )
    rows = [
        row for row in manifest_rows
        if row["extraction_eligible"].strip().lower() == "true"
    ]

    outdir = Path(args.output_dir)
    temporary = outdir.with_name(outdir.name + ".tmp")
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True)
    for tag in all_tags:
        (temporary / f"{tag}.fa").touch()
    writers = WriterCache(temporary)
    expected = Counter((row["global_tag_id"], row["sample_id"]) for row in rows)
    observed: Counter[tuple[str, str]] = Counter()
    records_out: list[dict[str, object]] = []
    try:
        rows_by_sample: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            if row["global_tag_id"] not in tag_set:
                raise ExtractError(f"Manifest TAG absent from matrix: {row['global_tag_id']}")
            rows_by_sample.setdefault(row["sample_id"], []).append(row)
        for sample, sample_rows in rows_by_sample.items():
            if sample not in fasta_map:
                raise ExtractError(f"No FASTA configured for sample {sample}")
            if not Path(fasta_map[sample]).is_file():
                raise ExtractError(f"Missing FASTA for {sample}: {fasta_map[sample]}")
            index = SeqIO.index(fasta_map[sample], "fasta")
            try:
                for row in sample_rows:
                    seqid = row["seqid"]
                    if seqid not in index:
                        raise ExtractError(f"{sample}: seqid {seqid} is absent from FASTA")
                    start, end = int(row["start"]), int(row["end"])
                    if start < 1 or end < start or end > len(index[seqid].seq):
                        raise ExtractError(f"Invalid interval for {row['record_id']}: {seqid}:{start}-{end}")
                    sequence = str(index[seqid].seq[start - 1 : end]).upper()
                    if row["strand"] == "-":
                        sequence = str(Seq(sequence).reverse_complement())
                    writers.write(row["global_tag_id"], row["record_id"], sequence)
                    observed[(row["global_tag_id"], sample)] += 1
                    records_out.append(
                        {
                            "record_id": row["record_id"],
                            "global_tag_id": row["global_tag_id"],
                            "sample_id": sample,
                            "gene_id": row["gene_id"],
                            "gene_name": row["gene_name"],
                            "copy_index": row["copy_index"],
                            "sequence_length": len(sequence),
                        }
                    )
            finally:
                index.close()
    finally:
        writers.close()
    if expected != observed:
        raise ExtractError("Eligible manifest counts do not equal extracted FASTA counts")

    records_path = Path(args.output_records)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    temp_records = records_path.with_suffix(records_path.suffix + ".tmp")
    fields = ["record_id", "global_tag_id", "sample_id", "gene_id", "gene_name", "copy_index", "sequence_length"]
    with temp_records.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records_out)
    if outdir.exists():
        shutil.rmtree(outdir)
    os.replace(temporary, outdir)
    os.replace(temp_records, records_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-records", required=True)
    try:
        run(parser.parse_args())
    except (OSError, ValueError, ExtractError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
