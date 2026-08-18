from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path

from pantag import prepare
from pantag.gff import load_schema, parse_attributes, transcript_ids


class GffLiftoffV2Test(unittest.TestCase):
    def test_all_transcript_aliases_are_returned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary) / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  transcript_id_keys: [transcript_id, ID]\n",
                encoding="utf-8",
            )
            schema = load_schema(config)
            attributes = parse_attributes(
                "ID=rna-NM_001180043.1;transcript_id=NM_001180043.1"
            )
            self.assertEqual(
                transcript_ids(attributes, schema),
                ("NM_001180043.1", "rna-NM_001180043.1"),
            )

    def test_gffread_internal_id_maps_to_parent_gene(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.yaml"
            config.write_text(
                "annotation:\n"
                "  mode: provided\n"
                "  gene_id_keys: [ID]\n"
                "  gene_name_keys: [Name]\n"
                "  transcript_id_keys: [transcript_id, ID]\n"
                "  transcript_parent_gene_keys: [Parent]\n",
                encoding="utf-8",
            )
            gff = root / "reference.gff3"
            gff.write_text(
                "chrI\ttest\tgene\t1\t300\t.\t+\t.\t"
                "ID=gene-YAL068C;Name=PAU8;gene_biotype=protein_coding\n"
                "chrI\ttest\tmRNA\t1\t300\t.\t+\t.\t"
                "ID=rna-NM_001180043.1;Parent=gene-YAL068C;"
                "transcript_id=NM_001180043.1\n",
                encoding="utf-8",
            )
            proteins = root / "raw.faa"
            proteins.write_text(">rna-NM_001180043.1\nMAAAA\n", encoding="utf-8")
            outputs = {name: root / name for name in ("genes.gff", "map.tsv", "proteins.faa", "qc.tsv")}
            prepare.command_tag_inputs(argparse.Namespace(
                sample_id="SGDref", config=str(config), annotation_gff=str(gff),
                raw_proteins=str(proteins), output_mcscanx_gff=str(outputs["genes.gff"]),
                output_gene_map=str(outputs["map.tsv"]), output_proteins=str(outputs["proteins.faa"]),
                output_qc=str(outputs["qc.tsv"]),
            ))
            self.assertIn(">gene-YAL068C", outputs["proteins.faa"].read_text(encoding="utf-8"))
            self.assertIn("unmatched_protein_records\t0", outputs["qc.tsv"].read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
