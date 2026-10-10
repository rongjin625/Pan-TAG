"""Regression tests for core refine and optional sequence/clustering stages."""

from __future__ import annotations

import argparse
import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam
import yaml

try:
    from pantag import fasta, refine
    from pantag.downstream import cluster, extract, structure
except ImportError as exc:
    cluster = extract = structure = fasta = refine = None
    IMPORT_ERROR = str(exc)
else:
    IMPORT_ERROR = ""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@unittest.skipIf(refine is None, f"downstream dependencies unavailable: {IMPORT_ERROR}")
class CoreRefineTest(unittest.TestCase):
    def test_structure_selection_and_gene_id_normalization(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            matrix = root / "matrix.tsv"
            matrix.write_text(
                "global_tag_id\tgene_names\tS1\tS2\tS3\tS4\n"
                "Pan_cluster_0001\tGeneA,GeneB\t2\t0\t*\t#\n",
                encoding="utf-8",
            )
            selected = structure.read_matrix(
                str(matrix), "Pan_cluster_0001", ["S1", "S2", "S3", "S4"]
            )
            self.assertEqual(selected, [(2, "S1")])
            self.assertEqual(
                structure.normalize_gene_identifier(
                    "S1_gene:ENSMUSG0001", "S1", ("gene:", "transcript:")
                ),
                "ENSMUSG0001",
            )

    def test_cluster_summary_merge_is_cumulative_and_unique(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first.tsv"
            second = root / "second.tsv"
            output = root / "global.tsv"

            def write_summary(path: Path, tag: str, status: str) -> None:
                row = {field: "" for field in cluster.SUMMARY_FIELDS}
                row.update({"global_tag_id": tag, "status": status})
                cluster.write_tsv(str(path), cluster.SUMMARY_FIELDS, [row])

            write_summary(first, "Pan_cluster_0002", "too_few_sequences")
            write_summary(second, "Pan_cluster_0001", "ok")
            cluster.merge_summaries(argparse.Namespace(
                inputs=[str(first), str(second), str(second)], output=str(output)
            ))
            rows = read_tsv(output)
            self.assertEqual(
                [row["global_tag_id"] for row in rows],
                ["Pan_cluster_0001", "Pan_cluster_0002"],
            )

    def test_structure_pdf_uses_completed_sequence_groups(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gff = root / "sample.gff3"
            gff.write_text(
                "chr1\ttest\tgene\t1\t10\t.\t+\t.\tID=gene:LEFT;Name=Left\n"
                "chr1\ttest\tgene\t21\t30\t.\t+\t.\tID=gene:TARGET1;Name=GeneA\n"
                "chr1\ttest\tgene\t41\t50\t.\t-\t.\tID=gene:TARGET2;Name=GeneB\n"
                "chr1\ttest\tgene\t61\t70\t.\t+\t.\tID=gene:RIGHT;Name=Right\n",
                encoding="utf-8",
            )
            samples = root / "samples.tsv"
            samples.write_text(
                f"sample_id\tfasta\tgff\nS1\tunused.fa\t{gff}\n",
                encoding="utf-8",
            )
            matrix = root / "matrix.tsv"
            matrix.write_text(
                "global_tag_id\tgene_names\tS1\n"
                "Pan_cluster_0001\tGeneA,GeneB\t2\n",
                encoding="utf-8",
            )
            coordinates = root / "coordinates.tsv"
            coordinates.write_text(
                "global_tag_id\trecord_id\tsample_id\tgene_id\tgene_name\tcopy_index\tsequence_group_id\tsequence_group_color\n"
                "Pan_cluster_0001\tr1\tS1\tTARGET1\tGeneA\t1\tPan_cluster_0001-SG01\t#112233\n"
                "Pan_cluster_0001\tr2\tS1\tTARGET2\tGeneB\t2\tPan_cluster_0001-SG02\t#445566\n",
                encoding="utf-8",
            )
            config = root / "config.yaml"
            template = yaml.safe_load(
                Path(__file__).with_name("config.yaml").read_text(encoding="utf-8")
            )
            template["annotation"]["mode"] = "provided"
            config.write_text(yaml.safe_dump(template, sort_keys=False), encoding="utf-8")
            gene_index = root / "structure_genes.sqlite"
            structure.build_index(argparse.Namespace(
                samples=str(samples), config=str(config), run_dir=str(root),
                liftoff_dir=str(root / "liftoff"), output=str(gene_index), threads=2,
            ))
            output = root / "TAG_structure.pdf"
            structure.run(argparse.Namespace(
                global_tag_id="Pan_cluster_0001", matrix=str(matrix),
                coordinates=str(coordinates), gene_index=str(gene_index),
                samples=str(samples), config=str(config),
                run_dir=str(root), liftoff_dir=str(root / "liftoff"),
                output=str(output), flank=1,
            ))
            self.assertGreater(output.stat().st_size, 0)

    def test_extract_spool_accepts_sequences_above_csv_default_limit(self) -> None:
        extract.allow_large_csv_fields()
        with tempfile.TemporaryDirectory() as temporary:
            spool = Path(temporary) / "long.tsv"
            sequence = "A" * 150000
            with spool.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=extract.SPOOL_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({
                    "record_id": "long", "global_tag_id": "Pan_cluster_0001",
                    "sample_id": "sample", "gene_id": "gene", "gene_name": "gene",
                    "copy_index": "1", "sequence_length": len(sequence),
                    "sequence": sequence,
                })
            with spool.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(row["sequence"]), len(sequence))

    def test_extract_threads_are_required(self) -> None:
        action = next(
            item for item in extract.build_parser()._actions
            if item.dest == "threads"
        )
        self.assertTrue(action.required)

    def test_cluster_threads_are_required(self) -> None:
        run_parser = next(
            action for action in cluster.build_parser()._actions
            if action.dest == "command"
        ).choices["run"]
        action = next(item for item in run_parser._actions if item.dest == "threads")
        self.assertTrue(action.required)

    def test_refine_flank_threshold_is_required(self) -> None:
        action = next(
            item for item in refine.build_parser()._actions
            if item.dest == "flank_overlap_threshold"
        )
        self.assertTrue(action.required)

    def test_span_removal_requires_fraction_strictly_above_point_nine(self) -> None:
        rows = [
            {"candidate_global_tag_id": "at_0_9", "sample_id": f"S{i}"}
            for i in range(9)
        ] + [
            {"candidate_global_tag_id": "above_0_9", "sample_id": f"S{i}"}
            for i in range(10)
        ]
        self.assertEqual(refine.span_candidates_to_remove(rows, 10), {"above_0_9"})

    def test_flank_score_uses_recovered_sample_flank_as_denominator(self) -> None:
        genes = [
            refine.GeneLite(str(i), name, name, i, i, "+")
            for i, name in enumerate(("L1", "L2", "TAG", "R1", "R2"), start=1)
        ]
        hit = refine.TagHit("tag", "TAG", "TAG", "chr1", 3, 3, "+", 2)
        self.assertEqual(
            refine.flank_score(genes, [hit], frozenset({"L1", "R1"}), 2),
            0.5,
        )
        isolated_hit = refine.TagHit("tag", "TAG", "TAG", "chr1", 3, 3, "+", 0)
        self.assertEqual(
            refine.flank_score([genes[2]], [isolated_hit], frozenset({"L1"}), 2),
            -1.0,
        )

    def test_refine_reads_mouse_style_pseudogene_features(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gff = root / "mouse.gff3"
            gff.write_text(
                "chr1\ttest\tgene\t1\t10\t.\t+\t.\t"
                "ID=gene:CODING;Name=Coding;biotype=protein_coding\n"
                "chr1\ttest\tpseudogene\t21\t30\t.\t+\t.\t"
                "ID=gene:MOUSE_PSEUDO;Name=MousePseudo;"
                "biotype=processed_pseudogene\n"
                # Merely using an allowed feature is insufficient without a
                # configured pseudogene biotype.
                "chr1\ttest\tpseudogene\t41\t50\t.\t+\t.\t"
                "ID=gene:NOT_PSEUDO;Name=NotPseudo;biotype=protein_coding\n"
                # Transcript-level records must not be counted as gene copies.
                "chr1\ttest\tpseudogenic_transcript\t61\t70\t.\t+\t.\t"
                "ID=transcript:PSEUDO_TX;Parent=gene:MOUSE_PSEUDO;"
                "Name=MousePseudo-201;biotype=processed_pseudogene\n",
                encoding="utf-8",
            )
            schema = refine.load_schema(Path(__file__).with_name("config.yaml"))
            chrom_genes, sample_tags = refine.read_gff_once(
                str(gff),
                schema,
                {
                    "Coding": ("candidate",),
                    "MousePseudo": ("candidate",),
                    "NotPseudo": ("candidate",),
                },
            )

            self.assertEqual(
                [gene.gene_id for gene in chrom_genes["chr1"]],
                ["CODING", "MOUSE_PSEUDO"],
            )
            self.assertEqual(
                [hit.gene_id for hit in sample_tags["candidate"]["chr1"]],
                ["CODING", "MOUSE_PSEUDO"],
            )

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

    def test_extract_selects_tags_and_uses_external_fasta_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta_path = root / "sample.fa"
            fasta_path.write_text(">chr1\nAACCGGTTAACC\n", encoding="utf-8")
            index_path = root / "indexes" / "sample.fai"
            fasta.build_index(fasta_path, index_path)
            samples = root / "samples.tsv"
            samples.write_text(
                "sample_id\tfasta\tgff\n"
                f"sample\t{fasta_path}\tunused.gff3\n"
                f"unselected\t{root / 'missing.fa'}\tunused.gff3\n",
                encoding="utf-8",
            )
            matrix = root / "matrix.tsv"
            matrix.write_text(
                "global_tag_id\nPan_cluster_0001\nPan_cluster_0002\n",
                encoding="utf-8",
            )
            manifest = root / "manifest.tsv"
            fields = [
                "record_id", "global_tag_id", "sample_id", "gene_id", "gene_name",
                "seqid", "start", "end", "strand", "copy_index", "qc_status",
                "extraction_eligible",
            ]
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows([
                    {
                        "record_id": "r1", "global_tag_id": "Pan_cluster_0001",
                        "sample_id": "sample", "gene_id": "g1", "gene_name": "G1",
                        "seqid": "chr1", "start": 1, "end": 4, "strand": "+",
                        "copy_index": 1, "qc_status": "pass", "extraction_eligible": "true",
                    },
                    {
                        "record_id": "r2", "global_tag_id": "Pan_cluster_0002",
                        "sample_id": "sample", "gene_id": "g2", "gene_name": "G2",
                        "seqid": "chr1", "start": 5, "end": 8, "strand": "-",
                        "copy_index": 1, "qc_status": "pass", "extraction_eligible": "true",
                    },
                ])
            output_dir = root / "per_TAG"
            records = root / "sequence_records.tsv"
            args = argparse.Namespace(
                manifest=str(manifest), matrix=str(matrix), samples=str(samples),
                reference_sample="sample", reference_fasta=str(fasta_path),
                fasta_indexes=[str(index_path), str(root / "missing.fai")],
                tags=["Pan_cluster_0002", "Pan_cluster_0002"],
                threads=2, output_dir=str(output_dir), output_records=str(records),
            )
            with manifest.open(encoding="utf-8", newline="") as handle:
                selected_rows = [
                    row for row in csv.DictReader(handle, delimiter="\t")
                    if row["global_tag_id"] == "Pan_cluster_0002"
                ]
            assets = extract.sample_assets(
                args.samples, args.reference_sample, args.reference_fasta,
                args.fasta_indexes,
            )
            temporary = output_dir.with_name(output_dir.name + ".tmp")
            spool_directory = temporary / ".sample_spools"
            spool_directory.mkdir(parents=True)
            extract.write_state(
                temporary / ".extract_state.json",
                extract.extraction_state(
                    args, ["Pan_cluster_0002"], assets, {"sample": selected_rows}
                ),
            )
            extract.extract_sample(
                "sample", selected_rows, str(fasta_path), str(index_path),
                spool_directory / "000000.tsv",
            )
            with patch.object(
                extract, "extract_sample",
                side_effect=AssertionError("completed spool should be reused"),
            ) as worker:
                extract.run(args)
            worker.assert_not_called()

            self.assertEqual(
                [path.name for path in output_dir.glob("*.fa")],
                ["Pan_cluster_0002.fa"],
            )
            self.assertEqual(
                (output_dir / "Pan_cluster_0002.fa").read_text(encoding="utf-8"),
                ">r2\nAACC\n",
            )
            self.assertEqual(
                [row["global_tag_id"] for row in read_tsv(records)],
                ["Pan_cluster_0002"],
            )

    def test_extract_rejects_unknown_selected_tags(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta_path = root / "sample.fa"
            fasta_path.write_text(">chr1\nAAAA\n", encoding="utf-8")
            index_path = root / "sample.fai"
            fasta.build_index(fasta_path, index_path)
            samples = root / "samples.tsv"
            samples.write_text(
                "sample_id\tfasta\tgff\n"
                f"sample\t{fasta_path}\tunused\n", encoding="utf-8",
            )
            matrix = root / "matrix.tsv"
            matrix.write_text("global_tag_id\nPan_cluster_0001\n", encoding="utf-8")
            manifest = root / "manifest.tsv"
            manifest.write_text("record_id\n", encoding="utf-8")
            with self.assertRaisesRegex(extract.ExtractError, "Pan_cluster_9999"):
                extract.run(argparse.Namespace(
                    manifest=str(manifest), matrix=str(matrix), samples=str(samples),
                    reference_sample="sample", reference_fasta=str(fasta_path),
                    fasta_indexes=[str(index_path)], tags=["Pan_cluster_9999"],
                    threads=1, output_dir=str(root / "out"),
                    output_records=str(root / "records.tsv"),
                ))

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
            template = yaml.safe_load(
                Path(__file__).with_name("config.yaml").read_text(encoding="utf-8")
            )
            annotation = dict(template["annotation"])
            annotation["mode"] = "provided"
            config.write_text(
                yaml.safe_dump({"annotation": annotation}, sort_keys=False),
                encoding="utf-8",
            )
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
                "Pan_cluster_0001", "Pan_cluster_0002", "Pan_cluster_0003",
                "Pan_cluster_0004", "Pan_cluster_0005",
            ]
            self.assertEqual([row["global_tag_id"] for row in final_tags], final_ids)
            self.assertEqual([row["global_tag_id"] for row in matrix], final_ids)
            self.assertEqual(list(matrix[0])[-2:], ["ref", "sample1"])
            self.assertEqual(
                sorted({row["global_tag_id"] for row in manifest}), final_ids[:4]
            )
            self.assertEqual((matrix[0]["ref"], matrix[0]["sample1"]), ("2", "2"))
            self.assertEqual((matrix[1]["ref"], matrix[1]["sample1"]), ("1", "1*"))
            self.assertEqual((matrix[2]["ref"], matrix[2]["sample1"]), ("1", "NA"))
            self.assertEqual((matrix[3]["ref"], matrix[3]["sample1"]), ("1", "2#"))
            self.assertEqual((matrix[4]["ref"], matrix[4]["sample1"]), ("NA", "NA"))
            self.assertEqual(matrix[1]["detected_sample_count"], "2")
            self.assertEqual(matrix[1]["gap_free_sample_count"], "1")
            self.assertEqual(matrix[1]["gap_free_sample_fraction"], "0.50")
            self.assertEqual(matrix[1]["gene_count_mean"], "1.00")
            self.assertEqual(matrix[1]["gene_count_std"], "0.00")
            self.assertEqual(matrix[1]["gene_count_cv"], "0.00")
            sample_columns = ["ref", "sample1"]
            self.assertFalse(any(
                row[sample] == "0"
                for row in matrix
                for sample in sample_columns
            ))
            self.assertEqual(
                [(row["candidate_global_tag_id"], row["global_tag_id"], row["status"]) for row in id_map],
                [
                    ("Candidate_cluster_0001", "Pan_cluster_0001", "retained"),
                    ("Candidate_cluster_0002", "Pan_cluster_0002", "retained"),
                    ("Candidate_cluster_0003", "Pan_cluster_0003", "retained"),
                    ("Candidate_cluster_0004", "Pan_cluster_0004", "retained"),
                    ("Candidate_cluster_0005", "", "removed_span"),
                    ("Candidate_cluster_0006", "Pan_cluster_0005", "retained"),
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
            self.assertFalse(any(
                row["global_tag_id"] == "Pan_cluster_0004" and row["sample_id"] == "sample1"
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
                fasta_indexes=[str(ref_fasta) + ".fai", str(sample_fasta) + ".fai"],
                tags=None, threads=2, output_dir=str(fasta_dir), output_records=str(records),
            ))
            self.assertEqual(sorted(path.stem for path in fasta_dir.glob("*.fa")), final_ids)

            coordinates = root / "pcoa_coordinates.tsv"
            cluster_summary = root / "cluster.tsv"
            pcoa = root / "pcoa.pdf"
            with patch.object(cluster.plt, "close") as close_plot:
                cluster.run_one(argparse.Namespace(
                    global_tag_id="Pan_cluster_0001", fasta=str(fasta_dir / "Pan_cluster_0001.fa"),
                    records=str(records), summary=str(cluster_summary), pcoa=str(pcoa),
                    pcoa_coordinates=str(coordinates), k=3, ubiquitous_threshold=0.95,
                    max_clusters=0, silhouette_tolerance=0.05, min_sequences=3,
                    seed=1, pcoa_jitter_frac=0.0,
                    pcoa_point_size=4.0, pcoa_alpha=0.25, threads=2,
                ))
                figure = close_plot.call_args.args[0]
            self.assertTrue(pcoa.is_file())
            self.assertEqual(len(figure.axes), 2)
            left_points = sorted(
                tuple(point) for collection in figure.axes[0].collections
                for point in collection.get_offsets()
            )
            right_points = sorted(
                tuple(point) for collection in figure.axes[1].collections
                for point in collection.get_offsets()
            )
            self.assertEqual(left_points, right_points)
            result = read_tsv(cluster_summary)[0]
            self.assertEqual(result["status"], "ok")
            self.assertEqual(result["ubiquitous_threshold"], "0.95")
            self.assertEqual(list(result), cluster.SUMMARY_FIELDS)
            coordinate_rows = read_tsv(coordinates)
            self.assertEqual(len(coordinate_rows), 4)
            self.assertEqual(
                list(coordinate_rows[0]), cluster.PCOA_FIELDS
            )
            self.assertNotIn("cluster", coordinate_rows[0])
            self.assertTrue(all(
                row["sequence_group_id"].startswith("Pan_cluster_0001-SG")
                for row in coordinate_rows
            ))
            colors_by_group: dict[str, set[str]] = {}
            for row in coordinate_rows:
                colors_by_group.setdefault(row["sequence_group_id"], set()).add(
                    row["sequence_group_color"]
                )
            self.assertTrue(all(len(colors) == 1 for colors in colors_by_group.values()))
            self.assertTrue(all(
                next(iter(colors)).startswith("#")
                for colors in colors_by_group.values()
            ))
            plotted_group_colors = {
                collection.get_label(): cluster.to_hex(collection.get_facecolor()[0])
                for collection in figure.axes[1].collections
            }
            self.assertEqual(
                plotted_group_colors,
                {group: next(iter(colors)) for group, colors in colors_by_group.items()},
            )
            cluster.plt.close(figure)

    def test_cluster_count_prefers_the_simplest_near_best_solution(self) -> None:
        self.assertEqual(
            cluster.choose_cluster_count({2: 0.84, 3: 0.86, 7: 0.88}, 0.05),
            2,
        )
        self.assertEqual(
            cluster.choose_cluster_count({2: 0.82, 3: 0.86, 7: 0.88}, 0.05),
            3,
        )

    def test_sequence_group_ids_and_colors_are_reproducible(self) -> None:
        records = [("record_b", "AAAA"), ("record_a", "CCCC"), ("record_c", "GGGG")]
        labels = cluster.np.asarray([9, 4, 9])
        first = cluster.sequence_group_metadata("Pan_cluster_0433", records, labels)
        second = cluster.sequence_group_metadata("Pan_cluster_0433", records, labels)
        self.assertEqual(first, second)
        self.assertEqual(first[0], ["Pan_cluster_0433-SG02", "Pan_cluster_0433-SG01", "Pan_cluster_0433-SG02"])


if __name__ == "__main__":
    unittest.main()
