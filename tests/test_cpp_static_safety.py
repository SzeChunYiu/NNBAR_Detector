from pathlib import Path
import re


def _constructor_body(source: str, class_name: str) -> str:
    match = re.search(
        rf"{class_name}::{class_name}\(\)\s*:.*?\{{(?P<body>.*?)\n\}}",
        source,
        flags=re.DOTALL,
    )
    assert match is not None, f"{class_name} default constructor not found"
    return match.group("body")


def _function_body(source: str, qualified_name: str) -> str:
    match = re.search(
        rf"void\s+{re.escape(qualified_name)}\(\)\s*\{{(?P<body>.*?)\n\}}",
        source,
        flags=re.DOTALL,
    )
    assert match is not None, f"{qualified_name} function not found"
    return match.group("body")


def test_nnbarhit_default_constructor_initializes_all_output_fields() -> None:
    body = _constructor_body(Path("src/hits/NNbarHit.cc").read_text(encoding="utf-8"), "NNbarHit")

    for field in [
        "posX",
        "posY",
        "posZ",
        "px",
        "py",
        "pz",
        "posX_particle",
        "posY_particle",
        "posZ_particle",
        "posX_local",
        "posY_local",
        "posZ_local",
        "TrackLength",
        "step_info",
    ]:
        assert re.search(rf"\b{field}\s*=", body), f"{field} is not initialized"


def test_main_argument_parser_checks_missing_option_values_before_argv_access() -> None:
    source = Path("src/main.cc").read_text(encoding="utf-8")

    assert "RequireOptionValue" in source
    assert "index + 1 < argc" in source
    assert "return 1;" in source


def test_gpu_shutdown_does_not_duplicate_run_end_statistics() -> None:
    for path, qualified_name in [
        ("src/gpu/GPUManager.cc", "GPUManager::Shutdown"),
        ("src/gpu/OpticalPhotonGPU.cc", "OpticalPhotonGPU::Shutdown"),
        ("src/gpu/TPCDriftGPU.cc", "TPCDriftGPU::Shutdown"),
    ]:
        body = _function_body(Path(path).read_text(encoding="utf-8"), qualified_name)
        assert "PrintStatistics()" not in body, f"{qualified_name} duplicates run-end statistics"

    run_action = Path("src/core/RunAction.cc").read_text(encoding="utf-8")
    assert "gpuMgr.PrintStatistics()" in run_action
    assert "opticalGPU.PrintStatistics()" in run_action
    assert "tpcGPU.PrintStatistics()" in run_action


def test_primary_generator_setters_update_shared_messenger_state() -> None:
    header = Path("include/core/PrimaryGeneratorAction.hh").read_text(encoding="utf-8")

    for setter, static_field in [
        ("SetCalibrationMode", "sCalibrationMode"),
        ("SetCalibrationEnergy", "sCalibEnergy"),
        ("SetCalibrationEnergyMin", "sCalibEnergyMin"),
        ("SetCalibrationEnergyMax", "sCalibEnergyMax"),
        ("SetCalibrationParticle", "sCalibParticle"),
        ("SetCalibrationSurface", "sCalibSurface"),
        ("SetCalibrationSpread", "sCalibSpread"),
    ]:
        assert re.search(rf"{setter}.*\{{\s*{static_field}\s*=", header), f"{setter} does not update {static_field}"

    assert "fCalibrationMode" not in header
    assert "fCalibEnergy" not in header
    assert "fCalibParticle" not in header


def test_primary_generator_signal_mode_can_emit_multiple_primaries_per_event() -> None:
    header = Path("include/core/PrimaryGeneratorAction.hh").read_text(encoding="utf-8")
    source = Path("src/core/PrimaryGeneratorAction.cc").read_text(encoding="utf-8")

    assert "sSignalParticles" in header
    assert 'DeclareProperty("signal_particles", sSignalParticles' in source
    assert "SplitSignalParticles(sSignalParticles)" in source
    assert "for (const G4String& signalParticle : signalParticles)" in source
    assert "event_number_global++" in source
    assert source.count("event_number_global++") == 3


def test_primary_generator_messenger_owns_only_initialized_commands() -> None:
    header = Path("include/core/PrimaryGeneratorMessenger.hh").read_text(encoding="utf-8")
    source = Path("src/core/PrimaryGeneratorMessenger.cc").read_text(encoding="utf-8")

    assert "PrimaryGeneratorAction*      Action" not in header
    assert "MessInput" not in header
    for member in ["CRYDir", "FileCmd", "InputCmd", "UpdateCmd"]:
        assert re.search(rf"\b{member}\s*=\s*nullptr", header), f"{member} is not initialized"
        assert f"delete {member};" in source


def test_cuda_only_gpu_buffers_are_not_compiled_into_cpu_fallback() -> None:
    optical_header = Path("include/gpu/OpticalPhotonGPU.hh").read_text(encoding="utf-8")
    tpc_header = Path("include/gpu/TPCDriftGPU.hh").read_text(encoding="utf-8")
    gpu_manager = Path("src/gpu/GPUManager.cc").read_text(encoding="utf-8")

    assert "#if HAS_CUDA\n    void* d_scintPoints_" in optical_header
    assert "#if HAS_CUDA\n    void* d_clusters_" in tpc_header
    assert "(void)bytes;" in gpu_manager
    assert "(void)ptr;" in gpu_manager


def test_leadglass_position_import_fails_closed_and_resets_rows() -> None:
    source = Path("src/detector/LeadGlass_geometry.cc").read_text(encoding="utf-8")

    assert '#include "G4Exception.hh"' in source
    assert "data.clear();" in source
    assert "G4ExceptionDescription" in source
    assert "LeadGlassMissingPositionFile" in source
    assert "LeadGlassEmptyPositionFile" in source
    assert "FatalException" in source
