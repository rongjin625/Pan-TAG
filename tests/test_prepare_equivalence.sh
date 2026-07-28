#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLE="${SAMPLE:-chm13}"
RAW_GFF="${RAW_GFF:-${BASE}/tests/data/reference/chm13.gff3}"
RAW_FA="${RAW_FA:-${BASE}/tests/data/reference/chm13.fa}"
OUT="${BASE}/comparison"

mkdir -p "${OUT}/old" "${OUT}/new" "${OUT}/shared"

python -m pantag.prepare validate-inputs \
    --config "${BASE}/config/config.yaml" \
    --samples "${BASE}/config/samples.tsv" \
    --output "${OUT}/input_validation.json"

gffread "${RAW_GFF}" \
    -g "${RAW_FA}" \
    -y "${OUT}/shared/${SAMPLE}.raw_protein.faa"

python "${BASE}/tests/legacy/pipeline1_get_clean_gff.py" \
    "${RAW_GFF}" \
    "${OUT}/old/${SAMPLE}.gff" \
    "${OUT}/old/${SAMPLE}.map"

python "${BASE}/tests/legacy/pipeline2_filter_longest_and_rename.py" \
    "${OUT}/shared/${SAMPLE}.raw_protein.faa" \
    "${OUT}/old/${SAMPLE}.gff" \
    "${RAW_GFF}" \
    "${OUT}/old/${SAMPLE}.proteins.faa"

python -m pantag.prepare tag-inputs \
    --sample-id "${SAMPLE}" \
    --annotation-gff "${RAW_GFF}" \
    --raw-proteins "${OUT}/shared/${SAMPLE}.raw_protein.faa" \
    --output-mcscanx-gff "${OUT}/new/${SAMPLE}.gff" \
    --output-gene-map "${OUT}/new/${SAMPLE}.map" \
    --output-proteins "${OUT}/new/${SAMPLE}.proteins.faa" \
    --output-qc "${OUT}/new/${SAMPLE}.prepare_qc.tsv"

cmp "${OUT}/old/${SAMPLE}.gff" "${OUT}/new/${SAMPLE}.gff"
cmp "${OUT}/old/${SAMPLE}.map" "${OUT}/new/${SAMPLE}.map"
cmp "${OUT}/old/${SAMPLE}.proteins.faa" "${OUT}/new/${SAMPLE}.proteins.faa"

echo "[RESULT] PASS: prepare outputs are byte-identical to pipeline1-2"
