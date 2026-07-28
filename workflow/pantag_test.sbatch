#!/usr/bin/env bash
#SBATCH --job-name=pantag_test
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/pantag_test.%j.out
#SBATCH --error=logs/pantag_test.%j.err

set -euo pipefail

if [[ -n "${PANTAG_ROOT:-}" ]]; then
    PROJECT_ROOT="${PANTAG_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
else
    echo "[ERROR] Cannot determine the Pan-TAG project directory." >&2
    echo "Set PANTAG_ROOT to the absolute Pan_TAG_v0.1 path." >&2
    exit 2
fi

if [[ ! -f "${PROJECT_ROOT}/pixi.toml" ]]; then
    echo "[ERROR] pixi.toml not found under PROJECT_ROOT=${PROJECT_ROOT}" >&2
    echo "Submit from the Pan_TAG_v0.1 root, or set PANTAG_ROOT explicitly." >&2
    exit 2
fi

cd "${PROJECT_ROOT}"
PIXI_MANIFEST="${PROJECT_ROOT}/pixi.toml"

if [[ -z "${PIXI_BIN:-}" ]]; then
    PIXI_BIN="$(command -v pixi || true)"
fi
if [[ -z "${PIXI_BIN}" || ! -x "${PIXI_BIN}" ]]; then
    echo "[ERROR] pixi executable not found." >&2
    echo "Submit with: sbatch --export=ALL,PIXI_BIN=/absolute/path/to/pixi workflow/pantag_test.sbatch" >&2
    exit 2
fi

export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

echo "Start:       $(date)"
echo "Host:        $(hostname)"
echo "Project:     ${PROJECT_ROOT}"
echo "Pixi:        ${PIXI_BIN}"
echo "Manifest:    ${PIXI_MANIFEST}"
echo "SLURM job:   ${SLURM_JOB_ID:-NA}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK:-NA}"

echo "[1/2] Environment preflight"
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" python --version
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" snakemake --version
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" diamond version
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" gffread --version
"${PIXI_BIN}" run \
    --manifest-path "${PIXI_MANIFEST}" \
    --environment liftoff \
    liftoff -h >/dev/null
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" g++ --version | head -1

echo "[2/2] Running Pan-TAG"
"${PIXI_BIN}" run --manifest-path "${PIXI_MANIFEST}" snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --cores "${SLURM_CPUS_PER_TASK:-8}" \
    --resources mem_mb=64000 \
    --printshellcmds \
    --show-failed-logs \
    --rerun-incomplete \
    --latency-wait 60

echo "End:         $(date)"
echo "[RESULT] Pan-TAG test completed"
