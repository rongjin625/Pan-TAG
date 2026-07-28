#!/usr/bin/env python3
"""Retrieve pseudogenes and finalize one sample's local TAG table.

This module packages the validated behavior of pipeline4--7.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class PseudogeneError(RuntimeError):
    """A user-facing pseudogene retrieval error."""


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def run_command(command: list[str]) -> None:
    print("[Pan-TAG]", " ".join(command), file=sys.stderr)
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise PseudogeneError(f"Executable not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise PseudogeneError(
            f"Command failed with exit code {exc.returncode}: {' '.join(command)}"
        ) from exc


def read_table(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise PseudogeneError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def tag_member_proteins(
    cluster_file: str | Path, proteins: str | Path, output: str | Path
) -> int:
    _, clusters = read_table(cluster_file)
    gene_to_cluster: dict[str, str] = {}
    for row in clusters:
        for gene_id in row["Members_IDs"].split(","):
            if gene_id:
                gene_to_cluster[gene_id] = row["Cluster_ID"]

    count = 0
    output_path = ensure_parent(output)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in SeqIO.parse(str(proteins), "fasta"):
            cluster_id = gene_to_cluster.get(record.id)
            if cluster_id:
                handle.write(f">{cluster_id}|{record.id}\n{record.seq}\n")
                count += 1
    return count


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def load_gene_records(
    gff: str | Path,
) -> tuple[
    dict[str, tuple[str, int, int]],
    dict[str, list[str]],
    dict[str, int],
    dict[str, str],
    list[dict[str, object]],
]:
    coordinates: dict[str, tuple[str, int, int]] = {}
    ordered_data: list[tuple[str, int, str]] = []
    names: dict[str, str] = {}
    pseudogenes: list[dict[str, object]] = []

    with Path(gff).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise PseudogeneError(f"Malformed GFF row {line_number}: {line.rstrip()}")
            if fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("ID") or attributes.get("gene_id")
            if not gene_id:
                continue
            gene_name = attributes.get("gene_name") or gene_id
            chrom, start, end, strand = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )
            coordinates[gene_id] = (chrom, start - 1, end)
            ordered_data.append((chrom, start, gene_id))
            if gene_name:
                names[gene_id] = re.sub(r"\.\d+$", "", gene_name.strip())
            if "pseudogene" in fields[8].lower():
                pseudogenes.append(
                    {
                        "id": gene_id,
                        "name": gene_name,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }
                )

    ordered_data.sort(key=lambda item: (item[0], item[1]))
    chrom_genes: dict[str, list[str]] = defaultdict(list)
    positions: dict[str, int] = {}
    for chrom, _, gene_id in ordered_data:
        positions[gene_id] = len(chrom_genes[chrom])
        chrom_genes[chrom].append(gene_id)
    return coordinates, dict(chrom_genes), positions, names, pseudogenes


def extract_nearby_pseudogenes(
    cluster_file: str | Path,
    pseudogenes: list[dict[str, object]],
    fasta: str | Path,
    output: str | Path,
    flank_bp: int,
) -> int:
    _, cluster_rows = read_table(cluster_file)
    search_regions = [
        {
            "id": row["Cluster_ID"],
            "chrom": row["Chrom"],
            "start": int(row["Start"]) - flank_bp,
            "end": int(row["End"]) + flank_bp,
        }
        for row in cluster_rows
    ]
    selected: list[tuple[str, dict[str, object]]] = []
    for pseudogene in pseudogenes:
        for cluster in search_regions:
            if (
                pseudogene["chrom"] == cluster["chrom"]
                and int(pseudogene["start"]) >= int(cluster["start"])
                and int(pseudogene["end"]) <= int(cluster["end"])
            ):
                selected.append((str(cluster["id"]), pseudogene))
                break

    fasta_index = SeqIO.index(str(fasta), "fasta")
    records: list[SeqRecord] = []
    try:
        for cluster_id, pseudogene in selected:
            chrom = str(pseudogene["chrom"])
            if chrom not in fasta_index:
                raise PseudogeneError(f"GFF seqid not found in FASTA: {chrom}")
            sequence = fasta_index[chrom].seq[
                int(pseudogene["start"]) - 1 : int(pseudogene["end"])
            ]
            if pseudogene["strand"] == "-":
                sequence = sequence.reverse_complement()
            record_id = (
                f"{cluster_id}|{pseudogene['id']}|"
                f"{pseudogene['name']}|{pseudogene['strand']}"
            )
            records.append(SeqRecord(sequence, id=record_id, description=""))
    finally:
        fasta_index.close()

    output_path = ensure_parent(output)
    SeqIO.write(records, str(output_path), "fasta")
    return len(records)


def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    return total + current_end - current_start + 1


def merge_hsp_coverage(raw_hits: str | Path, output: str | Path) -> int:
    pairs: dict[tuple[str, str], dict[str, object]] = {}
    path = Path(raw_hits)
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 14:
                    raise PseudogeneError(
                        f"blastx row {line_number} has {len(fields)} columns"
                    )
                key = (fields[0], fields[1])
                bitscore = float(fields[11])
                start, end = sorted((int(fields[6]), int(fields[7])))
                if key not in pairs:
                    pairs[key] = {
                        "bitscore": bitscore,
                        "row": fields,
                        "qlen": int(fields[12]),
                        "intervals": [(start, end)],
                    }
                else:
                    pairs[key]["intervals"].append((start, end))
                    if bitscore > float(pairs[key]["bitscore"]):
                        pairs[key]["bitscore"] = bitscore
                        pairs[key]["row"] = fields

    output_path = ensure_parent(output)
    with output_path.open("w", encoding="utf-8") as handle:
        for key in sorted(pairs):
            item = pairs[key]
            row = list(item["row"])
            coverage = 100 * merge_intervals(item["intervals"]) / int(item["qlen"])
            row[12] = str(round(coverage, 2))
            handle.write("\t".join(row) + "\n")
    return len(pairs)


def filter_blastx(
    merged_hits: str | Path,
    output: str | Path,
    min_identity: float,
    min_coverage: float,
) -> int:
    count = 0
    output_path = ensure_parent(output)
    with Path(merged_hits).open("r", encoding="utf-8") as source, output_path.open(
        "w", encoding="utf-8"
    ) as target:
        for line in source:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if float(fields[2]) >= min_identity and float(fields[12]) >= min_coverage:
                target.write(line if line.endswith("\n") else line + "\n")
                count += 1
    return count


def attach_pseudogenes(
    cluster_file: str | Path, blast_hits: str | Path
) -> tuple[list[str], list[dict[str, str]]]:
    cluster_to_pseudos: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with Path(blast_hits).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            query = fields[0].split("|")
            subject = fields[1].split("|")
            if len(query) >= 3 and len(subject) >= 2 and query[0] == subject[0]:
                cluster_to_pseudos[query[0]].add((query[1], query[2]))

    header, rows = read_table(cluster_file)
    output_header = header + ["Pseudo_Count", "Pseudos_Annotated", "Pseudos_IDs"]
    output_rows: list[dict[str, str]] = []
    for row in rows:
        pseudos = sorted(cluster_to_pseudos.get(row["Cluster_ID"], set()))
        result = dict(row)
        if pseudos:
            ids = [item[0] for item in pseudos]
            names = [item[1] for item in pseudos]
            counts = Counter(names)
            result["Pseudo_Count"] = str(len(pseudos))
            result["Pseudos_Annotated"] = ",".join(
                f"{name}({count})" if count > 1 else name
                for name, count in counts.items()
            )
            result["Pseudos_IDs"] = ",".join(ids)
        else:
            result.update(
                {
                    "Pseudo_Count": "0",
                    "Pseudos_Annotated": "None",
                    "Pseudos_IDs": "None",
                }
            )
        output_rows.append(result)
    return output_header, output_rows


def add_flanks(
    rows: list[dict[str, str]],
    chrom_genes: dict[str, list[str]],
    positions: dict[str, int],
    names: dict[str, str],
    flank_gene_count: int,
) -> None:
    for row in rows:
        ids = [
            gene_id
            for gene_id in (
                row["Members_IDs"] + "," + row["Pseudos_IDs"]
            ).split(",")
            if gene_id and gene_id != "None"
        ]
        chromosome = row["Chrom"]
        indices = sorted(
            positions[gene_id] for gene_id in ids if gene_id in positions
        )
        if chromosome not in chrom_genes or not indices:
            continue
        ordered = chrom_genes[chromosome]
        first, last = indices[0], indices[-1]
        upstream = [
            names[ordered[index]]
            for index in range(max(0, first - flank_gene_count), first)
            if ordered[index] in names
        ]
        downstream = [
            names[ordered[index]]
            for index in range(last + 1, min(len(ordered), last + flank_gene_count + 1))
            if ordered[index] in names
        ]
        row["Up_5th_ID"] = ",".join(upstream)
        row["Down_5th_ID"] = ",".join(downstream)


def strict_qc(
    rows: list[dict[str, str]],
    coordinates: dict[str, tuple[str, int, int]],
    fasta: str | Path,
    flank_gene_count: int,
) -> list[dict[str, str]]:
    fasta_index = SeqIO.index(str(fasta), "fasta")
    retained: list[dict[str, str]] = []
    try:
        for row in rows:
            upstream = [item for item in row["Up_5th_ID"].split(",") if item]
            downstream = [item for item in row["Down_5th_ID"].split(",") if item]
            if len(upstream) < flank_gene_count or len(downstream) < flank_gene_count:
                continue
            ids = [
                gene_id
                for gene_id in (
                    row["Members_IDs"] + "," + row["Pseudos_IDs"]
                ).split(",")
                if gene_id and gene_id != "None"
            ]
            valid = True
            for gene_id in ids:
                if gene_id not in coordinates:
                    valid = False
                    break
                chrom, start, end = coordinates[gene_id]
                if chrom not in fasta_index:
                    valid = False
                    break
                if "N" in str(fasta_index[chrom].seq[start:end]).upper():
                    valid = False
                    break
            if valid:
                retained.append(row)
    finally:
        fasta_index.close()
    return retained


def write_table(
    path: str | Path, header: list[str], rows: list[dict[str, str]]
) -> None:
    output = ensure_parent(path)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def command_run(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    tagged_proteins = workdir / f"{args.sample_id}.tagged_prots.faa"
    pseudo_fasta = workdir / f"{args.sample_id}.tagged_pseudos.fa"
    database = workdir / f"{args.sample_id}.pro.db"
    raw_hits = workdir / f"{args.sample_id}.blastx.raw"
    merged_hits = workdir / f"{args.sample_id}.blastx.raw.add"
    filtered_hits = workdir / f"{args.sample_id}.blastx"

    coordinates, chrom_genes, positions, names, pseudogenes = load_gene_records(
        args.gff
    )
    tag_member_proteins(args.cluster, args.proteins, tagged_proteins)
    pseudo_count = extract_nearby_pseudogenes(
        args.cluster, pseudogenes, args.fasta, pseudo_fasta, args.flank_bp
    )

    if pseudo_count:
        run_command(
            [
                args.diamond_bin,
                "makedb",
                "--in",
                str(tagged_proteins),
                "--db",
                str(database),
                "--quiet",
            ]
        )
        run_command(
            [
                args.diamond_bin,
                "blastx",
                "--query",
                str(pseudo_fasta),
                "--db",
                str(database),
                "--out",
                str(raw_hits),
                "--evalue",
                "1e-5",
                "-F",
                "15",
                "--max-hsps",
                "5",
                "--max-target-seqs",
                "5",
                "--threads",
                str(args.threads),
                "--range-culling",
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
                "qlen",
                "qcovhsp",
            ]
        )
    else:
        raw_hits.touch()

    merge_hsp_coverage(raw_hits, merged_hits)
    hit_count = filter_blastx(
        merged_hits,
        filtered_hits,
        args.min_identity,
        args.min_query_coverage,
    )
    header, rows = attach_pseudogenes(args.cluster, filtered_hits)
    add_flanks(rows, chrom_genes, positions, names, args.flank_gene_count)
    retained = strict_qc(rows, coordinates, args.fasta, args.flank_gene_count)
    write_table(args.output, header, retained)
    print(
        f"[Pan-TAG] {args.sample_id}: {pseudo_count} nearby pseudogene(s), "
        f"{hit_count} filtered hit(s), {len(retained)} valid TAG(s)."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--cluster", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--flank-bp", type=int, default=50000)
    parser.add_argument("--flank-gene-count", type=int, default=10)
    parser.add_argument("--min-identity", type=float, default=20.0)
    parser.add_argument("--min-query-coverage", type=float, default=20.0)
    parser.add_argument("--diamond-bin", default="diamond")
    parser.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return command_run(args)
    except (PseudogeneError, OSError, ValueError, KeyError) as exc:
        print(f"[Pan-TAG ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
