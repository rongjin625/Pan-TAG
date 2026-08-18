#!/usr/bin/env python3
"""Build one reusable gffutils database for all Liftoff sample jobs."""

from __future__ import annotations

import argparse
from pathlib import Path

import gffutils


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    gffutils.constants.ignore_url_escape_characters = True
    gffutils.create_db(
        args.gff,
        str(output),
        force=True,
        merge_strategy="create_unique",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        verbose=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
