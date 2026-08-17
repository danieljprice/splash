#!/usr/bin/env bash
# Run splash CLI render matrix for regression tests.
# Stages config prefix files into WORK_DIR, then writes PNG outputs there.
#
# Usage:
#   TESTDATA_DIR=/path/to/splash-testdata WORK_DIR=/path/to/work \
#     ./scripts/run_render_regression.sh
#
# Requires: splash on PATH (and giza libs as needed)

set -euo pipefail

TESTDATA_DIR="${TESTDATA_DIR:?Set TESTDATA_DIR to a splash-testdata checkout}"
WORK_DIR="${WORK_DIR:-.}"

DUMP="${DUMP_NAME:-binary_00000}"
PREFIX=render_testing
CONFIG_SRC="${TESTDATA_DIR}/config"
DUMP_SRC="${TESTDATA_DIR}/datafiles/${DUMP}"
# Phantom binary dump (override with SPLASH_FORMAT if needed)
FORMAT="${SPLASH_FORMAT:-phantom}"

if [[ ! -d "${CONFIG_SRC}" ]]; then
  echo "error: missing config directory: ${CONFIG_SRC}" >&2
  exit 1
fi

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

if [[ ! -f "${DUMP}" ]]; then
  if [[ -f "${DUMP_SRC}" ]]; then
    cp -f "${DUMP_SRC}" "./${DUMP}"
  else
    echo "error: dump file '${DUMP}' not found in ${WORK_DIR} or ${DUMP_SRC}" >&2
    exit 1
  fi
fi

if [[ ! -s "${DUMP}" ]]; then
  echo "error: dump file '${DUMP}' is empty" >&2
  exit 1
fi

# splash -p PREFIX looks for PREFIX.defaults/.limits/.units in the cwd
for ext in defaults limits units; do
  src="${CONFIG_SRC}/${PREFIX}.${ext}"
  if [[ ! -f "${src}" ]]; then
    echo "error: missing ${src}" >&2
    exit 1
  fi
  cp -f "${src}" "./${PREFIX}.${ext}"
done

if ! command -v splash >/dev/null 2>&1; then
  echo "error: splash not found on PATH" >&2
  exit 1
fi

run_one() {
  local out="$1"
  shift
  echo "+ splash $* -dev ${out}"
  splash "$@" -dev "${out}"
  if [[ ! -s "${out}" ]]; then
    echo "error: splash did not write ${out}" >&2
    exit 1
  fi
}

# CLI matrix (same coverage as PR #122)
run_one log_rho_render.png      -f "${FORMAT}" "${DUMP}" -r 6  -p "${PREFIX}"
run_one log_u_render.png        -f "${FORMAT}" "${DUMP}" -r 10 -p "${PREFIX}"
run_one log_rho_v_render.png    -f "${FORMAT}" "${DUMP}" -r 6  -p "${PREFIX}" -vec 7
run_one log_rho_sink0_render.png -f "${FORMAT}" "${DUMP}" -r 6  -p "${PREFIX}" --sink=0
run_one log_rho_sink1_render.png -f "${FORMAT}" "${DUMP}" -r 6  -p "${PREFIX}" --sink=1
run_one log_rho_sink2_render.png -f "${FORMAT}" "${DUMP}" -r 6  -p "${PREFIX}" --sink=2

echo "render regression: wrote PNGs in ${WORK_DIR}"
