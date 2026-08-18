"""Regression tests for the public Pan-TAG table schemas."""

from __future__ import annotations

import argparse
import csv
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from pantag import merge, prepare, retrieve_pseudogenes, single_sample_tag


def read_one(path: Path) -> tuple[list[str], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), next(reader)


class PublicSchemaTest(unittest.TestCase):
    def test_pseudogene_feature_types_require_gene_level_feature_and_biotype(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  mode: provided\n"
                "  gene_feature_types: [gene]\n"
                "  pseudogene_feature_types: [gene, pseudogene]\n"
                "  gene_id_keys: [gene_id, ID]\n"
                "  gene_name_keys: [gene_name, Name]\n"
                "  biotype_keys: [gene_type, gene_biotype, biotype]\n"
                "  pseudogene_values: [pseudogene, processed_pseudogene]\n",
                encoding="utf-8",
            )
            gff = root / "annotation.gff3"
            gff.write_text(
                # Human-style gene feature with a pseudogene biotype.
                "chr1\ttest\tgene\t1\t10\t.\t+\t.\t"
                "ID=gene:H_PSEUDO;Name=HP;gene_type=processed_pseudogene\n"
                # Mouse-style pseudogene feature with a pseudogene biotype.
                "chr1\ttest\tpseudogene\t21\t30\t.\t+\t.\t"
                "ID=gene:M_PSEUDO;Name=MP;biotype=processed_pseudogene\n"
                "chr1\ttest\tmRNA\t21\t30\t.\t+\t.\t"
                "ID=transcript:M_TX;Parent=gene:M_PSEUDO\n"
                # An allowed feature without a pseudogene biotype is ordinary only.
                "chr1\ttest\tgene\t41\t50\t.\t+\t.\t"
                "ID=gene:CODING;Name=PC;biotype=protein_coding\n"
                "chr1\ttest\tmRNA\t41\t50\t.\t+\t.\t"
                "ID=transcript:C_TX;Parent=gene:CODING\n"
                # Transcript-level pseudogene records are not allowed by default.
                "chr1\ttest\tpseudogenic_transcript\t61\t70\t.\t+\t.\t"
                "ID=transcript:TXP;Parent=gene:M_PSEUDO;Name=MP-201;"
                "biotype=processed_pseudogene\n",
                encoding="utf-8",
            )
            schema = retrieve_pseudogenes.load_schema(config)
            coordinates, chrom_genes, _, _, pseudogenes = (
                retrieve_pseudogenes.load_gene_records(gff, schema)
            )
            self.assertEqual(
                {item["id"] for item in pseudogenes},
                {"H_PSEUDO", "M_PSEUDO"},
            )
            self.assertEqual(set(coordinates), {"H_PSEUDO", "M_PSEUDO", "CODING"})
            self.assertEqual(chrom_genes["chr1"], ["CODING"])
            self.assertNotIn("TXP", coordinates)

            proteins = root / "raw.faa"
            proteins.write_text(">M_TX\nMAAA\n>C_TX\nMAAAAA\n", encoding="utf-8")
            outputs = {
                name: root / name
                for name in ("genes.gff", "map.tsv", "proteins.faa", "qc.tsv")
            }
            prepare.command_tag_inputs(argparse.Namespace(
                sample_id="mouse",
                config=str(config),
                annotation_gff=str(gff),
                raw_proteins=str(proteins),
                output_mcscanx_gff=str(outputs["genes.gff"]),
                output_gene_map=str(outputs["map.tsv"]),
                output_proteins=str(outputs["proteins.faa"]),
                output_qc=str(outputs["qc.tsv"]),
            ))
            with outputs["map.tsv"].open(encoding="utf-8", newline="") as handle:
                ordinary_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["Gene_ID"] for row in ordinary_rows],
                ["CODING"],
            )

    def test_local_tag_qc_checks_the_complete_member_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "sample.fa"
            # Both gene bodies are N-free; the N is only in the intergenic gap.
            fasta.write_text(">chr1\nAAAAANAAAAA\n", encoding="utf-8")
            row = {
                "Members_IDs": "gene1,gene2",
                "Pseudos_IDs": "None",
                "Upstream_Genes": "up1",
                "Downstream_Genes": "down1",
            }
            coordinates = {
                "gene1": ("chr1", 0, 5),
                "gene2": ("chr1", 6, 11),
            }
            retained, statuses = retrieve_pseudogenes.strict_qc(
                [row], coordinates, fasta, 1
            )
            self.assertEqual(retained, [])
            self.assertEqual(statuses["DROP_TAG_INTERVAL_CONTAINS_N"], 1)

            bad_coordinates = {"gene1": ("missing", 0, 5), "gene2": ("missing", 6, 11)}
            retained, statuses = retrieve_pseudogenes.strict_qc(
                [row], bad_coordinates, fasta, 1
            )
            self.assertEqual(retained, [])
            self.assertEqual(statuses["ERROR_TAG_INTERVAL_EXTRACTION"], 1)

    def test_merge_same_seqid_option_blocks_cross_seqid_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            header = (
                "Cluster_ID\tChromosome\tUpstream_Genes\tDownstream_Genes\t"
                "Members_Annotated\tMembers_IDs\tPseudos_Annotated\tPseudos_IDs\n"
            )
            reference = root / "ref_TAGs.tsv"
            sample = root / "sample_TAGs.tsv"
            reference.write_text(
                header + "R1\tchr1\tL\tR\tA,B\tidA,idB\tNone\tNone\n",
                encoding="utf-8",
            )
            sample.write_text(
                header + "S1\tchr2\tL\tR\tA,B\tidA2,idB2\tNone\tNone\n",
                encoding="utf-8",
            )

            def run(require_same_seqid: bool, suffix: str) -> int:
                tags = root / f"tags_{suffix}.tsv"
                mapping = root / f"mapping_{suffix}.tsv"
                merge.command_run(argparse.Namespace(
                    reference_sample="ref",
                    local_tags=[str(reference), str(sample)],
                    member_overlap_threshold=0.1,
                    require_same_seqid=require_same_seqid,
                    output_tags=str(tags),
                    output_mapping=str(mapping),
                ))
                with tags.open(encoding="utf-8", newline="") as handle:
                    return len(list(csv.DictReader(handle, delimiter="\t")))

            self.assertEqual(run(False, "unrestricted"), 1)
            self.assertEqual(run(True, "same_seqid"), 2)

            local = merge.LocalTag(
                sample="bridge", cluster_id="B1", chrom="chr1",
                member_names={"A", "B"}, member_ids={"a", "b"},
                pseudo_names=set(), pseudo_ids=set(),
                left_flanks={"L"}, right_flanks={"R"},
            )
            globals_ = [
                merge.GlobalTag("chr1", "ref", {"A"}, {"a"}, {"L"}, {"R"}),
                merge.GlobalTag("chr2", "ref", {"B"}, {"b"}, {"L"}, {"R"}),
            ]
            unrestricted = merge.candidate_indices(local, globals_, 0.1, False)
            restricted = merge.candidate_indices(local, globals_, 0.1, True)
            self.assertEqual([item["index"] for item in unrestricted], [0, 1])
            self.assertEqual([item["index"] for item in restricted], [0])

            name_index: dict[str, set[int]] = defaultdict(set)
            flank_index: dict[str, set[int]] = defaultdict(set)
            for index, global_tag in enumerate(globals_):
                merge.index_values(name_index, global_tag.names, index)
                merge.index_values(flank_index, global_tag.flanks, index)
            indexed = merge.candidate_indices(
                local,
                globals_,
                0.1,
                False,
                merge.indexed_candidates(local, name_index, flank_index),
            )
            self.assertEqual(indexed, unrestricted)

            extension = merge.LocalTag(
                sample="bridge", cluster_id="B2", chrom="chr1",
                member_names={"A", "C"}, member_ids={"a", "c"},
                pseudo_names=set(), pseudo_ids=set(),
                left_flanks={"L", "X"}, right_flanks={"R"},
            )
            merge.add_local_tag(
                globals_[0], extension, 1.0, 0, name_index, flank_index
            )
            probe = merge.LocalTag(
                sample="probe", cluster_id="P1", chrom="chr1",
                member_names={"C"}, member_ids={"c2"},
                pseudo_names=set(), pseudo_ids=set(),
                left_flanks={"X"}, right_flanks=set(),
            )
            self.assertEqual(
                merge.indexed_candidates(probe, name_index, flank_index), [0]
            )

    def test_blastp_requires_query_and_subject_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "blastp.raw"
            filtered = root / "blastp.filtered"
            raw.write_text(
                "q1\ts1\t50\t10\t0\t0\t1\t10\t1\t10\t1e-20\t100\t30\t25\n"
                "q2\ts2\t50\t10\t0\t0\t1\t10\t1\t10\t1e-20\t100\t30\t15\n",
                encoding="utf-8",
            )
            retained = single_sample_tag.filter_blastp(raw, filtered, 20, 20, 20)
            self.assertEqual(retained, 1)
            self.assertTrue(filtered.read_text(encoding="utf-8").startswith("q1\t"))

    def test_blastx_recomputes_union_coverage_on_both_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "blastx.raw"
            merged = root / "blastx.merged"
            filtered = root / "blastx.filtered"
            raw.write_text(
                "q1\ts1\t50\t20\t0\t0\t1\t20\t1\t10\t1e-20\t100\t100\t50\t20\t20\n"
                "q1\ts1\t50\t20\t0\t0\t21\t40\t11\t20\t1e-20\t90\t100\t50\t20\t20\n",
                encoding="utf-8",
            )
            retrieve_pseudogenes.merge_hsp_coverage(raw, merged)
            fields = merged.read_text(encoding="utf-8").strip().split("\t")
            self.assertEqual(fields[14:16], ["40.00", "40.00"])
            self.assertEqual(
                retrieve_pseudogenes.filter_blastx(merged, filtered, 20, 30, 30),
                1,
            )
            self.assertEqual(
                retrieve_pseudogenes.filter_blastx(merged, filtered, 20, 30, 45),
                0,
            )

    def test_plant_gff_coding_genes_are_selected_by_translated_protein(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  mode: provided\n"
                "  gene_id_keys: [ID]\n"
                "  gene_name_keys: [Name]\n"
                "  transcript_id_keys: [ID]\n"
                "  transcript_parent_gene_keys: [Parent]\n"
                "  fallback_gene_name_to_id: true\n",
                encoding="utf-8",
            )
            gff = root / "plant.gff3"
            gff.write_text(
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene:PLANT01;Name=ABC1\n"
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=transcript:TX01;Parent=gene:PLANT01\n"
                "chr1\ttest\tgene\t201\t300\t.\t+\t.\tID=gene:PLANT02\n"
                "chr1\ttest\tmRNA\t201\t300\t.\t+\t.\tID=transcript:TX02;Parent=gene:PLANT02\n"
                "chr1\ttest\tgene\t401\t500\t.\t+\t.\tID=gene:NONCODING;Name=NC\n",
                encoding="utf-8",
            )
            proteins = root / "raw.faa"
            proteins.write_text(">TX01\nMAAA\n>TX02\nMAAAAA\n", encoding="utf-8")
            outputs = {name: root / name for name in ("genes.gff", "map.tsv", "proteins.faa", "qc.tsv")}
            prepare.command_tag_inputs(argparse.Namespace(
                sample_id="plant", config=str(config), annotation_gff=str(gff),
                raw_proteins=str(proteins), output_mcscanx_gff=str(outputs["genes.gff"]),
                output_gene_map=str(outputs["map.tsv"]), output_proteins=str(outputs["proteins.faa"]),
                output_qc=str(outputs["qc.tsv"]),
            ))
            with outputs["map.tsv"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(row["Gene_ID"], row["Gene_Name"]) for row in rows],
                [("PLANT01", "ABC1"), ("PLANT02", "PLANT02")],
            )

    def test_gffread_id_can_match_any_configured_transcript_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  mode: provided\n"
                "  gene_id_keys: [ID]\n"
                "  gene_name_keys: [Name]\n"
                "  transcript_id_keys: [transcript_id, ID]\n"
                "  transcript_parent_gene_keys: [Parent]\n"
                "  coding_filter_mode: auto\n",
                encoding="utf-8",
            )
            gff = root / "yeast.gff3"
            gff.write_text(
                "chrI\ttest\tgene\t1\t300\t.\t+\t.\t"
                "ID=gene-YAL068C;Name=PAU8;gene_biotype=protein_coding\n"
                "chrI\ttest\tmRNA\t1\t300\t.\t+\t.\t"
                "ID=rna-NM_001180043.1;Parent=gene-YAL068C;"
                "transcript_id=NM_001180043.1\n",
                encoding="utf-8",
            )
            proteins = root / "raw.faa"
            proteins.write_text(
                ">rna-NM_001180043.1\nMAAAAA\n", encoding="utf-8"
            )
            outputs = {
                name: root / name
                for name in ("genes.gff", "map.tsv", "proteins.faa", "qc.tsv")
            }
            prepare.command_tag_inputs(argparse.Namespace(
                sample_id="yeast", config=str(config), annotation_gff=str(gff),
                raw_proteins=str(proteins), output_mcscanx_gff=str(outputs["genes.gff"]),
                output_gene_map=str(outputs["map.tsv"]), output_proteins=str(outputs["proteins.faa"]),
                output_qc=str(outputs["qc.tsv"]),
            ))
            self.assertEqual(
                outputs["proteins.faa"].read_text(encoding="utf-8").count(">"), 1
            )
            self.assertIn("retained_mcscanx_genes\t1", outputs["qc.tsv"].read_text(encoding="utf-8"))

    def test_auto_coding_filter_handles_human_and_mouse_readthrough(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  mode: provided\n"
                "  coding_filter_mode: auto\n"
                "  coding_values: [protein_coding]\n"
                "  excluded_attribute_terms: [readthrough_gene, readthrough]\n",
                encoding="utf-8",
            )
            genes = (
                ("H_OK", "gene_type=protein_coding", "HumanCoding"),
                ("H_RT", "gene_type=protein_coding;tag=readthrough_gene", "HumanRT"),
                ("M_OK", "biotype=protein_coding;description=ordinary gene", "MouseCoding"),
                ("M_RT", "biotype=protein_coding;description=xxx readthrough locus", "MouseRT"),
                ("M_NC", "biotype=lncRNA", "MouseNoncoding"),
            )
            gff_lines: list[str] = []
            protein_lines: list[str] = []
            for index, (gene_id, extra, name) in enumerate(genes, start=1):
                start = index * 200
                transcript = f"TX_{gene_id}"
                gff_lines.append(
                    f"chr1\ttest\tgene\t{start}\t{start + 99}\t.\t+\t.\t"
                    f"ID={gene_id};Name={name};{extra}\n"
                )
                gff_lines.append(
                    f"chr1\ttest\tmRNA\t{start}\t{start + 99}\t.\t+\t.\t"
                    f"ID={transcript};Parent={gene_id}\n"
                )
                protein_lines.append(f">{transcript}\nMAAAAA\n")

            gff = root / "annotation.gff3"
            proteins = root / "raw.faa"
            gff.write_text("".join(gff_lines), encoding="utf-8")
            proteins.write_text("".join(protein_lines), encoding="utf-8")
            outputs = {
                name: root / name
                for name in ("genes.gff", "map.tsv", "proteins.faa", "qc.tsv")
            }
            prepare.command_tag_inputs(argparse.Namespace(
                sample_id="mixed", config=str(config), annotation_gff=str(gff),
                raw_proteins=str(proteins), output_mcscanx_gff=str(outputs["genes.gff"]),
                output_gene_map=str(outputs["map.tsv"]), output_proteins=str(outputs["proteins.faa"]),
                output_qc=str(outputs["qc.tsv"]),
            ))

            with outputs["map.tsv"].open(encoding="utf-8", newline="") as handle:
                retained = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["Gene_ID"] for row in retained],
                ["H_OK", "M_OK"],
            )
            with outputs["qc.tsv"].open(encoding="utf-8", newline="") as handle:
                qc = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(qc["biotype_filter_applied"], "true")
            self.assertEqual(qc["noncoding_biotype_genes_filtered"], "1")
            self.assertEqual(qc["excluded_attribute_genes"], "2")

    def test_single_sample_flanks_are_gene_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mcscanx_gff = root / "sample.gff"
            tandem = root / "sample.tandem"
            gene_map = root / "sample.map"
            output = root / "sample_gene.cluster"
            mcscanx_gff.write_text(
                "".join(
                    f"chr1\tgene{index}\t{index * 100}\t{index * 100 + 49}\n"
                    for index in range(1, 13)
                ),
                encoding="utf-8",
            )
            tandem.write_text("gene6,gene7\n", encoding="utf-8")
            gene_map.write_text(
                "Gene_ID\tGene_Name\n"
                + "".join(
                    f"gene{index}\tNAME{index}\n" for index in range(1, 13)
                ),
                encoding="utf-8",
            )

            self.assertEqual(
                single_sample_tag.build_clusters(
                    mcscanx_gff, tandem, gene_map, output
                ),
                1,
            )
            columns, row = read_one(output)
            self.assertIn("Upstream_5th_Gene", columns)
            self.assertIn("Downstream_5th_Gene", columns)
            self.assertEqual(row["Upstream_5th_Gene"], "NAME1")
            self.assertEqual(row["Downstream_5th_Gene"], "NAME12")
            self.assertEqual(row["Length"], "150")

    def test_final_metrics_include_pseudogenes_and_exclude_them_from_density(self) -> None:
        coordinates = {
            "left_pseudo": ("chr1", 0, 10),
            "gene1": ("chr1", 10, 20),
            "internal_pseudo": ("chr1", 12, 18),
            "gene2": ("chr1", 30, 40),
            "right_pseudo": ("chr1", 40, 50),
        }
        chrom_genes = {"chr1": list(coordinates)}
        pseudo_ids = {"left_pseudo", "internal_pseudo", "right_pseudo"}
        cases = (
            ("left_pseudo", "1", "40", "40"),
            ("internal_pseudo", "11", "40", "30"),
            ("right_pseudo", "11", "50", "40"),
        )
        for pseudo, start, end, length in cases:
            with self.subTest(pseudo=pseudo):
                row = {
                    "Members_IDs": "gene1,gene2,gene1,None",
                    "Pseudos_IDs": pseudo,
                }
                retrieve_pseudogenes.finalize_tag_metrics(
                    [row], coordinates, chrom_genes, pseudo_ids
                )
                self.assertEqual(row["Size"], "3")
                self.assertEqual(row["Start"], start)
                self.assertEqual(row["End"], end)
                self.assertEqual(row["Length"], length)
                self.assertEqual(
                    row["Gene_Density_per_kb"], f"{2 / int(length) * 1000:.2f}"
                )

    def test_global_table_uses_compact_public_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs: list[str] = []
            header = (
                "Cluster_ID\tChromosome\tUpstream_Genes\tDownstream_Genes\t"
                "Members_Annotated\tMembers_IDs\tPseudos_Annotated\tPseudos_IDs\n"
            )
            for sample, member_ids, left, right in (
                ("ref", "idA1,idB1", "L1,L2", "R1,R2"),
                ("sample1", "idA2,idB2", "L2,L3", "R2,R3"),
            ):
                path = root / f"{sample}_TAGs.tsv"
                path.write_text(
                    header
                    + f"Cluster_1\tchr1\t{left}\t{right}\tA,B\t{member_ids}\tP\tpseudo_{sample}\n",
                    encoding="utf-8",
                )
                inputs.append(str(path))

            outputs = {
                name: root / f"{name}.tsv"
                for name in ("tags", "mapping")
            }
            merge.command_run(argparse.Namespace(
                reference_sample="ref",
                local_tags=inputs,
                member_overlap_threshold=0.1,
                output_tags=str(outputs["tags"]),
                output_mapping=str(outputs["mapping"]),
            ))

            columns, row = read_one(outputs["tags"])
            self.assertEqual(row["candidate_global_tag_id"], "Candidate_cluster_0001")
            self.assertEqual(row["pan_cluster_size"], "3")
            self.assertEqual(set(row["gene_names"].split(",")), {"A", "B", "P"})
            self.assertEqual(
                set(row["gene_ids"].split(",")),
                {"idA1", "idB1", "idA2", "idB2", "pseudo_ref", "pseudo_sample1"},
            )
            self.assertEqual(set(row["left_flank_genes"].split(",")), {"L1", "L2", "L3"})
            self.assertEqual(set(row["right_flank_genes"].split(",")), {"R1", "R2", "R3"})
            self.assertNotIn("local_tag_count", columns)
            self.assertNotIn("sample_ids", columns)
            self.assertNotIn("member_ids", columns)


if __name__ == "__main__":
    unittest.main()
