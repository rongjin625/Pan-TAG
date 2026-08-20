from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path

import yaml

from pantag import prepare
from pantag.gff import load_schema, parse_attributes, transcript_ids


def write_annotation_config(path: Path, **overrides: object) -> None:
    template = yaml.safe_load(
        Path(__file__).with_name("config.yaml").read_text(encoding="utf-8")
    )
    annotation = dict(template["annotation"])
    annotation.update(overrides)
    path.write_text(
        yaml.safe_dump({"annotation": annotation}, sort_keys=False),
        encoding="utf-8",
    )


class GffLiftoffV2Test(unittest.TestCase):
    def test_all_transcript_aliases_are_returned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary) / "config.yaml"
            write_annotation_config(
                config, transcript_id_keys=["transcript_id", "ID"]
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
            write_annotation_config(
                config,
                mode="provided",
                gene_id_keys=["ID"],
                gene_name_keys=["Name"],
                transcript_id_keys=["transcript_id", "ID"],
                transcript_parent_gene_keys=["Parent"],
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
