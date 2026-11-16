#!/usr/bin/env python3
"""Create FASTA indexes required by Pan-TAG random-access stages."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def build_index(fasta: str | Path, index: str | Path) -> None:
    source = Path(fasta).resolve()
    output = Path(index).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    # Keep pysam's output-capture files in private workspace storage. A fresh
    # interpreter also avoids Python's cached tempfile directory and changing
    # process-global settings when this function is called from other code.
    with tempfile.TemporaryDirectory(prefix=".faidx-", dir=output.parent) as scratch:
        link = Path(scratch) / ("source" + "".join(source.suffixes))
        os.symlink(source, link)
        environment = os.environ.copy()
        environment.update(TMPDIR=scratch, TEMP=scratch, TMP=scratch)
        subprocess.run(
            [sys.executable, "-c",
             "import sys, pysam; pysam.faidx(sys.argv[1])", str(link)],
            env=environment, check=True,
        )
        os.replace(Path(f"{link}.fai"), output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", nargs="+", required=True)
    parser.add_argument("--index", nargs="+", required=True)
    args = parser.parse_args()
    if len(args.fasta) != len(args.index):
        parser.error("--fasta and --index must contain the same number of paths")
    for fasta, index in zip(args.fasta, args.index):
        path = Path(fasta)
        if not path.is_file():
            parser.error(f"FASTA does not exist: {path}")
        build_index(path, index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
