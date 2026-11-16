"""Check index temporary-file isolation without requiring a local pysam build."""

import os
import importlib.util
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pantag.fasta import build_index


class FastaTemporaryDirectoryTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam is not installed")
    def test_real_index_with_unusable_inherited_temp_directory(self):
        import pysam

        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "genome.fa"
            source.write_bytes(b">chr1\nACGTACGT\n")
            output = Path(root) / "work" / "genome.fa.fai"
            unavailable = str(Path(root) / "missing" / "temp")
            with patch.dict(os.environ, {
                "TMPDIR": unavailable, "TEMP": unavailable, "TMP": unavailable,
            }):
                build_index(source, output)
                self.assertEqual(os.environ["TMPDIR"], unavailable)
            with pysam.FastaFile(str(source), filepath_index=str(output)) as handle:
                self.assertEqual(handle.fetch("chr1", 1, 5), "CGTA")
            self.assertFalse(Path(str(source) + ".fai").exists())
            self.assertEqual(list(output.parent.iterdir()), [output])

    def test_uses_private_directory_and_preserves_parent_environment(self):
        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "genome.fa"
            source.write_text(">chr1\nACGT\n")
            output = Path(root) / "work" / "genome.fa.fai"

            def run_index(command, *, env, check):
                scratch = Path(env["TMPDIR"])
                self.assertEqual(scratch.parent, output.parent)
                self.assertEqual(env["TMP"], str(scratch))
                self.assertEqual(env["TEMP"], str(scratch))
                self.assertTrue(check)
                self.assertEqual(Path(command[-1]).parent, scratch)
                Path(command[-1] + ".fai").write_text("chr1\t4\t6\t4\t5\n")

            with patch.dict(os.environ, {"TMPDIR": "/unavailable/login/tmp"}):
                with patch("pantag.fasta.os.symlink", side_effect=shutil.copyfile), \
                     patch("pantag.fasta.subprocess.run", side_effect=run_index):
                    build_index(source, output)
                self.assertEqual(os.environ["TMPDIR"], "/unavailable/login/tmp")
            self.assertTrue(output.is_file())
            self.assertEqual(list(output.parent.iterdir()), [output])
            self.assertFalse(Path(str(source) + ".fai").exists())

    def test_failure_cleans_scratch_and_preserves_existing_index(self):
        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "genome.fa"
            source.write_text(">chr1\nACGT\n")
            output = Path(root) / "genome.fa.fai"
            output.write_text("existing index")
            with patch("pantag.fasta.os.symlink", side_effect=shutil.copyfile), \
                 patch("pantag.fasta.subprocess.run",
                       side_effect=subprocess.CalledProcessError(1, "faidx")):
                with self.assertRaises(subprocess.CalledProcessError):
                    build_index(source, output)
            self.assertEqual(output.read_text(), "existing index")
            self.assertEqual(list(Path(root).glob(".faidx-*")), [])


if __name__ == "__main__":
    unittest.main()
