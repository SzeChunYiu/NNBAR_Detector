from pathlib import Path
import re


def _text(path: str) -> str:
    return Path(path).read_text(encoding="utf-8")


def test_cmake_exposes_and_links_cry_option() -> None:
    cmake = _text("CMakeLists.txt")

    assert re.search(r'option\(\s*WITH_CRY\s+"[^"]*"\s+OFF\s*\)', cmake)
    assert "find_path(CRY_INCLUDE_DIR CRYSetup.h" in cmake
    assert "$ENV{CRY_DIR}/src" in cmake or "${CRY_DIR}/src" in cmake
    assert "find_library(CRY_LIBRARY NAMES CRY" in cmake
    assert "$ENV{CRY_DIR}/lib" in cmake or "${CRY_DIR}/lib" in cmake
    assert "WITH_CRY_VALUE" in cmake
    assert "add_compile_definitions(WITH_CRY=1)" in cmake
    assert "${CRY_INCLUDE_DIR}" in cmake
    assert "target_link_libraries(nnbar-detector-simulation PRIVATE ${CRY_LIBRARY})" in cmake
    assert "#define WITH_CRY" in _text("config.h.in")


def test_cry_primary_generator_files_implement_thesis_parameters() -> None:
    header = _text("include/core/CRYPrimaryGenerator.hh")
    source = _text("src/core/CRYPrimaryGenerator.cc")

    assert "#ifdef WITH_CRY" in header
    assert "class CRYSetup;" in header
    assert "class CRYGenerator;" in header
    assert "void GenerateCRYPrimaries(G4Event* anEvent);" in header
    for method in [
        "SetParticleType",
        "SetEnergyMin",
        "SetEnergyMax",
        "SetEnergyBinIndex",
        "SetParticleIndex",
        "SetCRYDataPath",
    ]:
        assert method in header

    assert "date 1-1-2024" in source
    assert "latitude 55.71" in source
    assert "altitude 0" in source
    # CRY v1.7 reports particle positions in meters, so the 24 m x 24 m
    # thesis plane must be configured as a 24 m CRY subbox and converted to
    # Geant4 length units with CLHEP::m.  Treating CRY's coordinates as cm
    # shrinks the generated plane by a factor of 100 in the Parquet output.
    assert re.search(r"subboxLength\s+24(?!\d)", source)
    assert "selected->x() * CLHEP::m" in source
    assert "selected->y() * CLHEP::m" in source
    assert "selected->x() * CLHEP::cm" not in source
    assert "selected->y() * CLHEP::cm" not in source
    assert "500.0 * CLHEP::cm" in source
    assert "G4UniformRand() * (fEmax - fEmin)" in source
    assert "N_ij[6][5]" in source
    assert "1.69e11" in source and "2.30e12" in source and "5.00e6" in source
    assert re.search(r"for\s*\(\s*int\s+i\s*=\s*0\s*;\s*i\s*<\s*6\s*;\s*\+\+i\s*\)", source)
    assert "return (n_ij / S) * (n_ij / sum_i);" in source
    assert "ParquetOutputManager::Instance().WriteParticle" in source
    assert "rec.weight = fWeight" in source


def test_cry_setup_is_lazy_until_macro_configuration_is_applied() -> None:
    header = _text("include/core/CRYPrimaryGenerator.hh")
    source = _text("src/core/CRYPrimaryGenerator.cc")

    # /cosmic/dataPath is processed after /run/initialize, so construction must
    # not open CRY's data files with the default ./data path.
    constructor_match = re.search(
        r"CRYPrimaryGenerator::CRYPrimaryGenerator\(\).*?\{(?P<body>.*?)\n\}",
        source,
        re.DOTALL,
    )
    assert constructor_match is not None
    assert "UpdateCRYSetup()" not in constructor_match.group("body")

    assert "bool fSetupDirty = true;" in header
    assert "fSetupDirty = true;" in source
    assert "if (!fGenerator || fSetupDirty)" in source
    assert "fSetupDirty = false;" in source


def test_cry_particle_matching_uses_non_const_cry_api_signature() -> None:
    header = _text("include/core/CRYPrimaryGenerator.hh")
    source = _text("src/core/CRYPrimaryGenerator.cc")

    # CRY v1.7 exposes CRYParticle::PDGid() as a non-const method, so the
    # matcher must not take a const CRYParticle reference.
    assert "bool MatchesRequestedParticle(CRYParticle& particle) const;" in header
    assert "bool CRYPrimaryGenerator::MatchesRequestedParticle(CRYParticle& particle) const" in source


def test_cry_setup_is_lazy_until_macro_configuration_is_applied() -> None:
    header = _text("include/core/CRYPrimaryGenerator.hh")
    source = _text("src/core/CRYPrimaryGenerator.cc")

    constructor = re.search(
        r"CRYPrimaryGenerator::CRYPrimaryGenerator\(\)[^{]*\{(?P<body>.*?)\n\}",
        source,
        re.S,
    )
    assert constructor is not None
    assert "UpdateCRYSetup();" not in constructor.group("body")
    assert "bool fSetupDirty = true;" in header
    assert "fSetupDirty = true;" in source
    assert "if (!fGenerator || fSetupDirty)" in source
    assert "fSetupDirty = false;" in source


def test_cry_generator_uses_preallocated_event_vector_api() -> None:
    source = _text("src/core/CRYPrimaryGenerator.cc")

    assert "std::vector<CRYParticle*> particles;" in source
    assert "fGenerator->genEvent(&particles);" in source
    assert "fGenerator->genEvent()" not in source


def test_primary_generator_action_has_cosmic_messenger_and_dispatch() -> None:
    header = _text("include/core/PrimaryGeneratorAction.hh")
    source = _text("src/core/PrimaryGeneratorAction.cc")

    assert "CRYPrimaryGenerator" in header
    assert "fCRYGenerator" in header
    assert "sCRYMode" in header
    assert "sCRYParticle" in header
    assert "sCRYEnergyMin" in header
    assert "sCRYEnergyMax" in header

    assert 'new G4GenericMessenger(this, "/cosmic/"' in source
    for command in [
        'DeclareProperty("mode", sCRYMode',
        'DeclareProperty("particle", sCRYParticle',
        'DeclarePropertyWithUnit("energyMin", "GeV", sCRYEnergyMin',
        'DeclarePropertyWithUnit("energyMax", "GeV", sCRYEnergyMax',
        'DeclareProperty("energyBin", sCRYEnergyBinIdx',
        'DeclareProperty("particleIdx", sCRYParticleIdx',
        'DeclareProperty("dataPath", sCRYDataPath',
    ]:
        assert command in source

    assert "if (sCRYMode && fCRYGenerator)" in source
    assert "fCRYGenerator->GenerateCRYPrimaries(anEvent);" in source


def test_cosmic_slurm_array_maps_30_bins_and_skips_zero_weight_bins() -> None:
    slurm = _text("slurm/run_cosmic.slurm")

    assert "#SBATCH --array=0-29" in slurm
    assert "#SBATCH -p lu48" in slurm
    assert "#SBATCH --cpus-per-task=4" in slurm
    assert "#SBATCH --mem=16G" in slurm
    assert "#SBATCH -t 08:00:00" in slurm
    assert "#SBATCH -o /projects/hep/fs10/shared/nnbar/billy/NNBAR_Detector_sim/slurm/cosmic-%A_%a.out" in slurm
    assert "#SBATCH -e /projects/hep/fs10/shared/nnbar/billy/NNBAR_Detector_sim/slurm/cosmic-%A_%a.err" in slurm
    assert "activate_hibeam_env()" in slurm
    assert 'for f in "${CONDA_PREFIX}"/etc/conda/activate.d/*.sh; do' in slurm
    assert "activate_hibeam_env" in slurm
    assert "export CRY_DIR" in slurm
    assert "PARTICLES=(mu- gamma e- neutron proton)" in slurm
    assert "EMIN=(0 0.5 1 5 10 50)" in slurm
    assert "EMAX=(0.5 1 5 10 50 200)" in slurm
    for skipped in ["gamma:4", "gamma:5", "e-:4", "e-:5"]:
        assert skipped in slurm
    assert "/generator/use_particle_gun" not in slurm
    assert "/cosmic/mode true" in slurm
    assert "/cosmic/particleIdx" in slurm
    assert "/run/beamOn 1000000" in slurm
    assert "/generator/use_particle_gun true" not in slurm
    assert './nnbar-detector-simulation -m "${macro_file}" -t 4 -g' in slurm


def test_runtime_wrapper_sources_conda_geant4_data_environment() -> None:
    wrapper = _text("scripts/wrapper_template.sh.in")

    assert "CONDA_PREFIX=\"@GEANT4_INSTALL@\"" in wrapper
    assert "etc/conda/activate.d" in wrapper
    assert 'for f in "${CONDA_PREFIX}"/etc/conda/activate.d/*.sh; do' in wrapper
    assert 'source "${f}"' in wrapper


def test_cosmic_slurm_activates_geant4_runtime_data() -> None:
    slurm = _text("slurm/run_cosmic.slurm")

    assert "activate_hibeam_env()" in slurm
    assert 'for f in "${CONDA_PREFIX}"/etc/conda/activate.d/*.sh' in slurm
    assert 'source "${f}"' in slurm
    assert 'export Geant4_DIR="${CONDA_PREFIX}/lib/cmake/Geant4"' in slurm
    assert 'export CRY_DIR="${CRY_DIR}"' in slurm
    assert 'ln -s "${NNBAR_ROOT}/data" data' in slurm
    assert 'ln -s "${NNBAR_ROOT}/config" config' in slurm
    assert 'ln -s "${NNBAR_ROOT}/macro" macro' in slurm
