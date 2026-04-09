#!/usr/bin/env bash

set -euo pipefail

# input Argument (from Condor submit file)
PARTICLE=$1
ENERGY=$2   # input energy ex. lgE_16.0
SIN_VAL=$3  # input zenith ex. 0.1

SIN="sin2_${SIN_VAL}"

OUTDIR="Xmax"
mkdir -p "$OUTDIR"

BASE_LONG="/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/${PARTICLE}/${ENERGY}"
BASE_ZEN="./ZenithAngle"

echo "--- Current Task: ${PARTICLE} | ${ENERGY} | ${SIN} ---"

echo "--- Checking local files ---"
ls -F
echo "----------------------------"

LONGBASENAME="${BASE_LONG}/${SIN}"
ZENFILE="${BASE_ZEN}/${PARTICLE}_${ENERGY}_${SIN_VAL}.npz"
OUTFILE="${OUTDIR}/${PARTICLE}_${ENERGY}_${SIN}.dat"

if [[ ! -f "${ZENFILE}" ]]; then
    echo "Skipping: Zenith file missing"
    echo "${ZENFILE}"
    exit 0
fi

if [[ ! -d "${LONGBASENAME}" ]]; then
    echo "Skipping: Longbase directory missing"
    exit 0
fi

mapfile -t long_files < <(find "${LONGBASENAME}" -name "DAT*.long" -size +0c 2>/dev/null)

if [[ ${#long_files[@]} -eq 0 ]]; then
    echo "Skipping: No .long files found"
    exit 0
fi

echo "Processing ${PARTICLE} ${ENERGY} ${SIN} (${#long_files[@]} files)"

python3 submit_Xmax_ParseAndFitLongitudinalProfile.py \
    --longbase "${LONGBASENAME}" \
    --zenfile "${ZENFILE}" \
    --outfile "${OUTFILE}" \
    --removeFinal20gcm2



echo "Task completed!"