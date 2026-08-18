#!/usr/bin/env bash
set -euo pipefail

# Build Pan-TAG's modified MCScanX from the bundled source archive.
#
# Usage:
#   build.sh OUTPUT_BIN BUILD_INFO
#
# Snakefile supplies both paths from config.work_dir. Optional defaults keep
# the script convenient for a manual build in the source checkout.
#
# The temporary extracted source directory is deleted after a successful build.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

SOURCE_ARCHIVE="${SCRIPT_DIR}/MCScanX-source.tar.gz"
OUTPUT_BIN="${1:-${PROJECT_ROOT}/work/software/MCScanX_pantag}"
BUILD_INFO="${2:-${PROJECT_ROOT}/work/software/MCScanX_pantag.build_info.txt}"
SOFTWARE_DIR="$(dirname "${OUTPUT_BIN}")"
BUILD_DIR="${SOFTWARE_DIR}/MCScanX-build"

if [[ "$(dirname "${BUILD_INFO}")" != "${SOFTWARE_DIR}" ]]; then
    echo "[ERROR] OUTPUT_BIN and BUILD_INFO must use the same directory." >&2
    exit 2
fi

CXX_BIN="${CXX:-g++}"

echo "[1/6] Checking files and compiler"

if [[ ! -f "${SOURCE_ARCHIVE}" ]]; then
    echo "[ERROR] Missing source archive: ${SOURCE_ARCHIVE}" >&2
    exit 2
fi

if ! command -v "${CXX_BIN}" >/dev/null 2>&1; then
    echo "[ERROR] Cannot find compiler: ${CXX_BIN}" >&2
    echo "Load a GCC module first, then rerun this script." >&2
    exit 2
fi

"${CXX_BIN}" --version | head -1

echo "[2/6] Preparing temporary build directory"
mkdir -p "${SOFTWARE_DIR}"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

tar -xzf "${SOURCE_ARCHIVE}" -C "${BUILD_DIR}"

# The archive may contain files directly, or one top-level directory.
if [[ -f "${BUILD_DIR}/msa.cc" ]]; then
    SOURCE_DIR="${BUILD_DIR}"
else
    MSA_MATCHES=()
    while IFS= read -r -d '' path; do
        MSA_MATCHES+=("$path")
    done < <(find "${BUILD_DIR}" -type f -name msa.cc -print0)

    if [[ "${#MSA_MATCHES[@]}" -ne 1 ]]; then
        echo "[ERROR] Expected exactly one msa.cc after extraction; found ${#MSA_MATCHES[@]}." >&2
        exit 2
    fi

    SOURCE_DIR="$(dirname "${MSA_MATCHES[0]}")"
fi

echo "[3/6] Checking required MCScanX source files"
REQUIRED_FILES=(
    struct.cc struct.h
    mcscan.cc mcscan.h
    read_data.cc read_data.h
    out_utils.cc out_utils.h
    dagchainer.cc dagchainer.h
    msa.cc msa.h
    permutation.cc permutation.h
)

for name in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "${SOURCE_DIR}/${name}" ]]; then
        echo "[ERROR] Missing source file: ${SOURCE_DIR}/${name}" >&2
        exit 2
    fi
done

echo "[4/6] Changing tandem gene-order distance from 1 to 15"
python - "${SOURCE_DIR}/msa.cc" <<'PY'
from pathlib import Path
import sys

msa = Path(sys.argv[1])
text = msa.read_text(encoding="utf-8")

old = (
    "if (fabs(it8->second.gene_id-it9->second.gene_id)==1"
    "&&it8->second.mol==it9->second.mol)"
)
new = (
    "if (fabs(it8->second.gene_id-it9->second.gene_id)<=15"
    "&&it8->second.mol==it9->second.mol)"
)

count = text.count(old)
if count != 1:
    raise SystemExit(
        f"[ERROR] Expected the original tandem condition exactly once; found {count}."
    )

msa.write_text(text.replace(old, new, 1), encoding="utf-8")

if msa.read_text(encoding="utf-8").count(new) != 1:
    raise SystemExit("[ERROR] Failed to install the <=15 condition.")

print("[PASS] msa.cc now uses tandem gene-order distance <=15")
PY

echo "[5/6] Compiling MCScanX_pantag"
(
    cd "${SOURCE_DIR}"
    "${CXX_BIN}" -O2 \
        struct.cc \
        mcscan.cc \
        read_data.cc \
        out_utils.cc \
        dagchainer.cc \
        msa.cc \
        permutation.cc \
        -o MCScanX_pantag
)

if [[ ! -x "${SOURCE_DIR}/MCScanX_pantag" ]]; then
    echo "[ERROR] Compilation failed: MCScanX_pantag was not created." >&2
    exit 2
fi

cp "${SOURCE_DIR}/MCScanX_pantag" "${OUTPUT_BIN}"
chmod +x "${OUTPUT_BIN}"

echo "[6/6] Recording build information and cleaning temporary source"
{
    echo "program=MCScanX_pantag"
    echo "source_archive=${SOURCE_ARCHIVE}"
    echo "modification=tandem_gene_order_distance_1_to_15"
    echo "compiler=${CXX_BIN}"
    echo "compiler_path=$(command -v "${CXX_BIN}")"
    echo "compiler_version=$("${CXX_BIN}" --version | head -1)"
    echo "binary_sha256=$(sha256sum "${OUTPUT_BIN}" | awk '{print $1}')"
} > "${BUILD_INFO}"

rm -rf "${BUILD_DIR}"

echo
echo "[RESULT] PASS"
echo "Binary:     ${OUTPUT_BIN}"
echo "Build info: ${BUILD_INFO}"
echo "Temporary source removed: ${BUILD_DIR}"
