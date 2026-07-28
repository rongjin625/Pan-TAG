#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLE="${SAMPLE:-chm13}"
INPUT="${BASE}/comparison/new"
WORK_DIR="${BASE}/work/mcscanx_test/${SAMPLE}"
RESULT_DIR="${BASE}/results/local_tags/${SAMPLE}"

mkdir -p "${WORK_DIR}" "${RESULT_DIR}"

if [[ ! -x "${BASE}/work/software/MCScanX_pantag" ]]; then
    bash "${BASE}/src/pantag/mcscanx/build.sh"
fi

python -m pantag.single_sample_tag \
    --sample-id "${SAMPLE}" \
    --workdir "${WORK_DIR}" \
    --mcscanx-bin "${BASE}/work/software/MCScanX_pantag" \
    --mcscanx-gff "${INPUT}/${SAMPLE}.gff" \
    --gene-map "${INPUT}/${SAMPLE}.map" \
    --proteins "${INPUT}/${SAMPLE}.proteins.faa" \
    --threads "${THREADS:-4}" \
    --min-identity 20 \
    --min-query-coverage 20 \
    --output-blast-raw "${WORK_DIR}/${SAMPLE}.blast.raw" \
    --output-blast "${WORK_DIR}/${SAMPLE}.blast" \
    --output-tandem "${WORK_DIR}/${SAMPLE}.tandem" \
    --output-cluster "${RESULT_DIR}/${SAMPLE}_gene.cluster"

test -f "${WORK_DIR}/${SAMPLE}.tandem"
test -s "${RESULT_DIR}/${SAMPLE}_gene.cluster"

echo "[RESULT] PASS: single-sample TAG outputs generated"
