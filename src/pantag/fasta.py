#!/usr/bin/env python3
"""Create FASTA indexes required by Pan-TAG random-access stages."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pysam


def build_index(fasta: str | Path, index: str | Path) -> None:
    source = Path(fasta).resolve()
    output = Path(index)
    output.parent.mkdir(parents=True, exist_ok=True)
    link = output.with_name(output.name + ".source" + "".join(source.suffixes))
    generated = Path(f"{link}.fai")
    temporary = output.with_suffix(output.suffix + ".tmp")
    for path in (link, generated, temporary):
        path.unlink(missing_ok=True)
    try:
        os.symlink(source, link)
        pysam.faidx(str(link))
        os.replace(generated, temporary)
        os.replace(temporary, output)
    finally:
        link.unlink(missing_ok=True)
        generated.unlink(missing_ok=True)
        temporary.unlink(missing_ok=True)


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
