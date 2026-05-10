#!/bin/bash
#SBATCH -J nnbar-cry-cosmic
#SBATCH -A lu2026-2-51
#SBATCH -p lu48
#SBATCH --array=0-29
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o /projects/hep/fs10/shared/nnbar/billy/NNBAR_Detector_sim/slurm/cosmic-%A_%a.out
#SBATCH -e /projects/hep/fs10/shared/nnbar/billy/NNBAR_Detector_sim/slurm/cosmic-%A_%a.err

set -euo pipefail

NNBAR_ROOT=${NNBAR_ROOT:-/projects/hep/fs10/shared/nnbar/billy/NNBAR_Detector_sim}
BUILD_DIR=${BUILD_DIR:-${NNBAR_ROOT}/build_lunarc}
CRY_DIR=${CRY_DIR:-/projects/hep/fs10/shared/nnbar/billy/packages/cry_v1.7}
HIBENV=${HIBENV:-/projects/hep/fs10/shared/nnbar/billy/packages/hibeam_env}
EVENTS=${EVENTS:-1000000}  # default macro command: /run/beamOn 1000000

activate_hibeam_env() {
    export CONDA_PREFIX="${HIBENV}"
    export CONDA_DEFAULT_ENV="hibeam_env"
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    set +u
    for f in "${CONDA_PREFIX}"/etc/conda/activate.d/*.sh; do
        [ -f "${f}" ] && source "${f}"
    done
    set -u
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
    export Geant4_DIR="${CONDA_PREFIX}/lib/cmake/Geant4"
    export CRY_DIR="${CRY_DIR}"
}

echo "[cosmic] Activating conda-forge Geant4 runtime..."
activate_hibeam_env

PARTICLES=(mu- gamma e- neutron proton)
EBINS=(0 1 2 3 4 5)
EMIN=(0 0.5 1 5 10 50)
EMAX=(0.5 1 5 10 50 200)

particle_idx=$((SLURM_ARRAY_TASK_ID / 6))
ebin=$((SLURM_ARRAY_TASK_ID % 6))
particle=${PARTICLES[$particle_idx]}

# Zero N_ij bins from thesis Table 6.1: gamma:4 gamma:5 e-:4 e-:5
case "${particle}:${ebin}" in
  gamma:4|gamma:5|e-:4|e-:5)
    echo "Skipping zero-weight bin ${particle}:${ebin}"
    exit 0
    ;;
esac

out_dir="output/cosmic_${particle}_${ebin}"
macro_file="${BUILD_DIR}/cosmic_${particle}_${ebin}.mac"
mkdir -p "${BUILD_DIR}/output" "${BUILD_DIR}/slurm" "${BUILD_DIR}/${out_dir}"

cd "${BUILD_DIR}"
rm -rf data config macro
ln -s "${NNBAR_ROOT}/data" data
ln -s "${NNBAR_ROOT}/config" config
ln -s "${NNBAR_ROOT}/macro" macro

cat > "${macro_file}" <<MACRO
/run/initialize
/particle_generator/set_folder_name cosmic_${particle}_${ebin}
/cosmic/mode true
/cosmic/particle ${particle}
/cosmic/energyMin ${EMIN[$ebin]} GeV
/cosmic/energyMax ${EMAX[$ebin]} GeV
/cosmic/energyBin ${ebin}
/cosmic/particleIdx ${particle_idx}
/cosmic/dataPath ${CRY_DIR}/data
/run/beamOn ${EVENTS}
MACRO

./nnbar-detector-simulation -m "${macro_file}" -t 4 -g
