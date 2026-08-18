"""Regression tests for core refine and optional sequence/clustering stages."""

from __future__ import annotations

import argparse
import csv
import tempfile
import unittest
from pathlib import Path

import pysam

try:
    from pantag import fasta, refine
    from pantag.downstream import cluster, extract
except ImportError as exc:
    cluster = extract = fasta = refine = None
    IMPORT_ERROR = str(exc)
else:
    IMPORT_ERROR = ""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@unittest.skipIf(refine is None, f"downstream dependencies unavailable: {IMPORT_ERROR}")
class CoreRefineTest(unittest.TestCase):
    def test_fasta_index_can_be_stored_outside_read_only_input_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            input_dir.mkdir()
            source = input_dir / "genome.fa"
            source.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
            index = root / "work" / "fasta_indexes" / "sample.fai"
            fasta.build_index(source, index)
            self.assertTrue(index.is_file())
            self.assertFalse(Path(f"{source}.fai").exists())
            handle = pysam.FastaFile(str(source), filepath_index=str(index))
            try:
                self.assertEqual(handle.fetch("chr1", 1, 5), "CGTA")
            finally:
                handle.close()

    def test_refine_filters_span_renumbers_and_feeds_downstream(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ref_fasta = root / "ref.fa"
            sample_fasta = root / "sample.fa"
            ref_chr1 = list("A" * 200)
            ref_chr1[20:30] = list("C" * 10)
            ref_fasta.write_text(
                ">chr1\n" + "".join(ref_chr1) + "\n>chr2\n" + "A" * 200 + "\n",
                encoding="utf-8",
            )
            sample_chr1 = list("A" * 200)
            sample_chr1[20:30] = list("C" * 10)
            sample_chr1[60:70] = list("N" * 10)
            sample_fasta.write_text(">chr1\n" + "".join(sample_chr1) + "\n>chr2\n" + "A" * 200 + "\n", encoding="utf-8")
            pysam.faidx(str(ref_fasta))
            pysam.faidx(str(sample_fasta))

            def gene(seqid: str, name: str, start: int) -> str:
                return f"{seqid}\ttest\tgene\t{start}\t{start + 9}\t.\t+\t.\tID={seqid}_{name}_{start};Name={name}\n"

            ref_gff = root / "ref.gff3"
            sample_gff = root / "sample.gff3"
            ref_gff.write_text(
                "".join([
                    gene("chr1", "A", 1), gene("chr1", "B", 21),
                    gene("chr1", "C", 61), gene("chr1", "D", 81),
                    gene("chr1", "E", 101), gene("chr1", "F", 121),
                    gene("chr1", "X", 141), gene("chr1", "G", 161),
                    gene("chr2", "X", 1), gene("chr2", "G", 21),
                ]), encoding="utf-8",
            )
            sample_gff.write_text(
                "".join([
                    gene("chr1", "A", 1), gene("chr1", "B", 21),
                    gene("chr1", "C", 61), gene("chr1", "E", 101),
                    gene("chr1", "F", 121), gene("chr1", "X", 141),
                    gene("chr1", "G", 161), gene("chr2", "E", 1),
                    gene("chr2", "F", 21), gene("chr2", "X", 41),
                    gene("chr2", "G", 61),
                ]), encoding="utf-8",
            )
            samples = root / "samples.tsv"
            samples.write_text(
                "sample_id\tfasta\tgff\n"
                f"ref\t{ref_fasta}\t{ref_gff}\n"
                f"sample1\t{sample_fasta}\t{sample_gff}\n",
                encoding="utf-8",
            )
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "candidate_global_tag_id\tpan_cluster_size\treference_seqid\tsource_sample_id\toccupancy\tgene_names\tgene_ids\tleft_flank_genes\tright_flank_genes\n"
                "Candidate_cluster_0001\t2\tchr1\tref\t2\tA,B\ta,b\t\tC\n"
                "Candidate_cluster_0002\t1\tchr1\tref\t2\tC\tc\tB\t\n"
                "Candidate_cluster_0003\t1\tchr1\tref\t1\tD\td\tC\t\n"
                "Candidate_cluster_0004\t1\tchr1\tref\t2\tE\te\t\tF\n"
                "Candidate_cluster_0005\t1\tchr1\tref\t2\tG\tg\tX\t\n"
                "Candidate_cluster_0006\t1\tchr1\tref\t2\tF\tf\tZ\t\n",
                encoding="utf-8",
            )
            config = root / "config.yaml"
            config.write_text("annotation:\n  mode: provided\n", encoding="utf-8")
            outputs = {
                key: root / value for key, value in {
                    "tags": "pangenome_TAGs.tsv",
                    "matrix": "refined_TAG_count_matrix.tsv",
                    "manifest": "TAG_sequence_manifest.tsv",
                    "map": "refined_cluster_id_map.tsv",
                    "span": "span_removed_clusters.tsv",
                }.items()
            }
            refine.run(argparse.Namespace(
                tags=str(candidates), config=str(config), samples=str(samples),
                reference_sample="ref", reference_fasta=str(ref_fasta),
                reference_gff=str(ref_gff), target_gff_dir=str(root / "liftoff"),
                flank_gene_count=1, flank_overlap_threshold=0.20, threads=2,
                workdir=str(root / "refine_work"),
                output_tags=str(outputs["tags"]), output_matrix=str(outputs["matrix"]),
                output_manifest=str(outputs["manifest"]), output_id_map=str(outputs["map"]),
                output_span_audit=str(outputs["span"]),
            ))

            final_tags = read_tsv(outputs["tags"])
            matrix = read_tsv(outputs["matrix"])
            manifest = read_tsv(outputs["manifest"])
            id_map = read_tsv(outputs["map"])
            spans = read_tsv(outputs["span"])
            final_ids = [
                "Pan_cluster_0001", "Pan_cluster_0002",
                "Pan_cluster_0003", "Pan_cluster_0004",
            ]
            self.assertEqual([row["global_tag_id"] for row in final_tags], final_ids)
            self.assertEqual([row["global_tag_id"] for row in matrix], final_ids)
            self.assertEqual(list(matrix[0])[-2:], ["ref", "sample1"])
            self.assertEqual(
                sorted({row["global_tag_id"] for row in manifest}), final_ids[:3]
            )
            self.assertEqual((matrix[0]["ref"], matrix[0]["sample1"]), ("2", "2"))
            self.assertEqual((matrix[1]["ref"], matrix[1]["sample1"]), ("1", "1(gap)"))
            self.assertEqual((matrix[2]["ref"], matrix[2]["sample1"]), ("1", "NA"))
            self.assertEqual((matrix[3]["ref"], matrix[3]["sample1"]), ("NA", "NA"))
            self.assertEqual(matrix[1]["detected_sample_count"], "2")
            self.assertEqual(matrix[1]["gap_free_sample_count"], "1")
            self.assertEqual(matrix[1]["gap_free_sample_fraction"], "0.50")
            self.assertEqual(matrix[1]["gene_count_mean"], "1.00")
            self.assertEqual(matrix[1]["gene_count_std"], "0.00")
            self.assertEqual(matrix[1]["gene_count_cv"], "0.00")
            sample_columns = ["ref", "sample1"]
            self.assertFalse(any(
                "(span)" in row[sample] or row[sample] == "0"
                for row in matrix
                for sample in sample_columns
            ))
            self.assertEqual(
                [(row["candidate_global_tag_id"], row["global_tag_id"], row["status"]) for row in id_map],
                [
                    ("Candidate_cluster_0001", "Pan_cluster_0001", "retained"),
                    ("Candidate_cluster_0002", "Pan_cluster_0002", "retained"),
                    ("Candidate_cluster_0003", "Pan_cluster_0003", "retained"),
                    ("Candidate_cluster_0004", "", "removed_span"),
                    ("Candidate_cluster_0005", "", "removed_span"),
                    ("Candidate_cluster_0006", "Pan_cluster_0004", "retained"),
                ],
            )
            self.assertEqual(len([row for row in spans if row["candidate_global_tag_id"] == "Candidate_cluster_0004"]), 1)
            self.assertEqual(len([row for row in spans if row["candidate_global_tag_id"] == "Candidate_cluster_0005"]), 2)
            self.assertTrue(all(row["reason"] == "MULTIPLE_FLANK_SUPPORTED_BLOCKS" for row in spans))
            self.assertTrue(all(row["valid_block_count"] == "2" for row in spans))
            self.assertTrue(all(row["block_coordinates"] for row in spans))
            self.assertFalse(any("Candidate_cluster" in row["record_id"] for row in manifest))
            self.assertFalse(any(row["qc_status"] != "pass" for row in manifest))
            self.assertFalse(any(
                row["global_tag_id"] == "Pan_cluster_0002" and row["sample_id"] == "sample1"
                for row in manifest
            ))
            for sample_id in ("ref", "sample1"):
                for suffix in ("counts.tsv", "manifest.tsv", "span.tsv"):
                    self.assertTrue((root / "refine_work" / f"{sample_id}.{suffix}").is_file())

            fasta_dir = root / "per_TAG"
            records = root / "sequence_records.tsv"
            extract.run(argparse.Namespace(
                manifest=str(outputs["manifest"]), matrix=str(outputs["matrix"]),
                samples=str(samples), reference_sample="ref", reference_fasta=str(ref_fasta),
                output_dir=str(fasta_dir), output_records=str(records),
            ))
            self.assertEqual(sorted(path.stem for path in fasta_dir.glob("*.fa")), final_ids)

            assignments, cross, cluster_summary = root / "assignments.tsv", root / "cross.tsv", root / "cluster.tsv"
            cluster.run_one(argparse.Namespace(
                global_tag_id="Pan_cluster_0001", fasta=str(fasta_dir / "Pan_cluster_0001.fa"),
                records=str(records), assignments=str(assignments), cross_table=str(cross),
                summary=str(cluster_summary), pcoa="", k=3, ubiquitous_threshold=0.95,
                max_clusters=0, min_sequences=3, seed=1, pcoa_jitter_frac=0.0,
                pcoa_point_size=4.0, pcoa_alpha=0.25,
            ))
            result = read_tsv(cluster_summary)[0]
            self.assertEqual(result["status"], "ok")
            self.assertEqual(result["ubiquitous_threshold"], "0.95")


if __name__ == "__main__":
    unittest.main()
