#!/usr/bin/env bash
#SBATCH --job-name=pantag_test
#SBATCH --partition=cpu64,cpu128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/pantag_test.%j.out
#SBATCH --error=logs/pantag_test.%j.err

set -euo pipefail

CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_ROOT="${SLURM_SUBMIT_DIR:-${PWD}}"
CONFIG_PATH="${PANTAG_CONFIG:-${RUN_ROOT}/config/config.yaml}"
if [[ "${CONFIG_PATH}" != /* ]]; then
    CONFIG_PATH="${RUN_ROOT}/${CONFIG_PATH}"
fi
SNAKEMAKE_BIN="${CODE_ROOT}/.pixi/envs/default/bin/snakemake"
test -x "${SNAKEMAKE_BIN}"

PYTHONPATH="${CODE_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
"${SNAKEMAKE_BIN}" "${PANTAG_TARGET:-all}" \
  --directory "${RUN_ROOT}" \
  --snakefile "${CODE_ROOT}/workflow/Snakefile" \
  --configfile "${CONFIG_PATH}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}"
