from __future__ import annotations

import argparse
import csv
import tempfile
import unittest
from pathlib import Path

import yaml

from pantag.cli import ConfigureError, build_parser, configure


class ConfigureCliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        for name in ("ref.fa", "ref.gff3", "sample.fa", "sample.gff3"):
            (self.root / name).write_text("test\n", encoding="utf-8")
        self.features = self.root / "features.txt"
        self.features.write_text("gene\n", encoding="utf-8")
        self.template = self.root / "template.yaml"
        template_config = yaml.safe_load(
            Path(__file__).with_name("config.yaml").read_text(encoding="utf-8")
        )
        template_config["reference"] = {
            "sample_id": "old", "fasta": "old", "gff": "old"
        }
        template_config["samples"] = "old.tsv"
        template_config["annotation"]["gene_id_keys"] = ["ID"]
        template_config["liftoff"]["feature_types"] = str(self.features)
        self.template.write_text(
            yaml.safe_dump(template_config, sort_keys=False),
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_samples(self, target_gff: bool) -> Path:
        path = self.root / "samples.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "fasta", "gff"])
            writer.writerow(
                ["REF", self.root / "ref.fa", self.root / "ref.gff3"]
            )
            writer.writerow(
                [
                    "S1",
                    self.root / "sample.fa",
                    self.root / "sample.gff3" if target_gff else "",
                ]
            )
        return path

    def args(self, samples: Path, mode: str) -> argparse.Namespace:
        return argparse.Namespace(
            samples=str(samples),
            annotation_mode=mode,
            pseudogene_flank_bp=25000,
            flank_gene_count=6,
            template=str(self.template),
            output=str(self.root / "config.yaml"),
        )

    def test_configure_provided_generates_workflow_schema(self) -> None:
        samples = self.write_samples(target_gff=True)
        output = configure(self.args(samples, "provided"))
        config = yaml.safe_load(output.read_text(encoding="utf-8"))

        self.assertEqual(config["reference"]["sample_id"], "REF")
        self.assertEqual(config["reference"]["fasta"], str(self.root / "ref.fa"))
        self.assertEqual(config["samples"], str(samples))
        self.assertEqual(config["annotation"]["mode"], "provided")
        self.assertEqual(config["annotation"]["gene_id_keys"], ["ID"])
        self.assertEqual(config["pseudogene"]["search_flank_bp"], 25000)
        self.assertEqual(config["flanks"]["gene_count"], 6)
        self.assertEqual(config["merge"]["member_overlap_threshold"], 0.3)

    def test_pseudogene_feature_types_cli_override(self) -> None:
        samples = self.write_samples(target_gff=True)
        args = build_parser().parse_args([
            "--samples", str(samples),
            "--annotation-mode", "provided",
            "--pseudogene-enabled",
            "--pseudogene-feature-types", "gene,pseudogene",
            "--template", str(self.template),
            "--output", str(self.root / "pseudo-features.yaml"),
        ])
        config = yaml.safe_load(configure(args).read_text(encoding="utf-8"))
        self.assertEqual(
            config["annotation"]["pseudogene_feature_types"],
            ["gene", "pseudogene"],
        )

    def test_incomplete_template_does_not_receive_hidden_defaults(self) -> None:
        old_template = self.root / "old-template.yaml"
        config = yaml.safe_load(self.template.read_text(encoding="utf-8"))
        config["annotation"].pop("pseudogene_feature_types")
        old_template.write_text(
            yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
        )
        args = build_parser().parse_args([
            "--samples", str(self.write_samples(target_gff=True)),
            "--annotation-mode", "provided",
            "--pseudogene-enabled",
            "--template", str(old_template),
            "--output", str(self.root / "upgraded.yaml"),
        ])
        with self.assertRaisesRegex(ConfigureError, "pseudogene_feature_types"):
            configure(args)

    def test_configure_liftoff_accepts_empty_target_gff(self) -> None:
        output = configure(self.args(self.write_samples(target_gff=False), "liftoff"))
        config = yaml.safe_load(output.read_text(encoding="utf-8"))
        self.assertEqual(config["annotation"]["mode"], "liftoff")

    def test_configure_provided_rejects_missing_target_gff(self) -> None:
        with self.assertRaisesRegex(ConfigureError, "requires target GFFs"):
            configure(self.args(self.write_samples(target_gff=False), "provided"))

    def test_grouped_cli_overrides_without_resource_options(self) -> None:
        samples = self.write_samples(target_gff=True)
        output = self.root / "full-config.yaml"
        args = build_parser().parse_args([
            "--samples", str(samples),
            "--annotation-mode", "provided",
            "--template", str(self.template),
            "--output", str(output),
            "--reference-sample", "REF",
            "--require-same-seqid",
            "--blastp-subject-coverage", "25",
            "--flank-overlap-threshold", "0.35",
            "--no-pseudogene-enabled",
            "--results-dir", "custom-results",
        ])
        configure(args)
        config = yaml.safe_load(output.read_text(encoding="utf-8"))
        self.assertTrue(config["merge"]["require_same_seqid"])
        self.assertEqual(config["tag"]["blastp_subject_coverage"], 25.0)
        self.assertEqual(config["refine"]["flank_overlap_threshold"], 0.35)
        self.assertFalse(config["pseudogene"]["enabled"])
        self.assertEqual(config["results_dir"], "custom-results")

    def test_help_shows_config_values_without_none_defaults(self) -> None:
        parser = build_parser()
        help_text = parser.format_help()
        self.assertNotIn("(default: None)", help_text)
        self.assertIn("(config default:", help_text)
        self.assertNotIn("Resource options", help_text)

    def test_custom_template_values_are_not_overwritten(self) -> None:
        samples = self.write_samples(target_gff=True)
        args = build_parser().parse_args([
            "--samples", str(samples),
            "--annotation-mode", "provided",
            "--template", str(self.template),
            "--output", str(self.root / "custom-config.yaml"),
            "--pseudogene-enabled",
        ])
        output = configure(args)
        config = yaml.safe_load(output.read_text(encoding="utf-8"))
        self.assertEqual(config["pseudogene"]["search_flank_bp"], 50000)
        self.assertEqual(config["flanks"]["gene_count"], 10)

    def test_invalid_refine_flank_overlap_threshold_is_rejected(self) -> None:
        config = yaml.safe_load(self.template.read_text(encoding="utf-8"))
        config["refine"] = {"flank_overlap_threshold": 1.01}
        self.template.write_text(
            yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
        )
        with self.assertRaisesRegex(ConfigureError, "must be between 0 and 1"):
            configure(self.args(self.write_samples(target_gff=True), "provided"))


if __name__ == "__main__":
    unittest.main()
