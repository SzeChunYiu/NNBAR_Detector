from pathlib import Path
import re


def _cmake_option_default(cmake_text: str, option_name: str) -> str:
    match = re.search(
        rf"^\s*option\(\s*{re.escape(option_name)}\s+\"[^\"]*\"\s+(ON|OFF)\s*\)",
        cmake_text,
        flags=re.MULTILINE,
    )
    assert match is not None, f"{option_name} option is missing"
    return match.group(1)


def test_gpu_backends_are_explicit_opt_in_by_default() -> None:
    cmake_text = Path("CMakeLists.txt").read_text(encoding="utf-8")

    assert _cmake_option_default(cmake_text, "WITH_CELERITAS") == "OFF"
    assert _cmake_option_default(cmake_text, "WITH_OPTICKS") == "OFF"


def test_build_script_resets_optional_physics_flags_by_default() -> None:
    build_script = Path("scripts/build.sh").read_text(encoding="utf-8")

    for option in [
        "WITH_SCINTILLATION",
        "WITH_GARFIELD",
        "WITH_GARFIELD_GPU",
        "WITH_CELERITAS",
        "WITH_OPTICKS",
        "WITH_DASHBOARD",
    ]:
        assert f"-D{option}=OFF" in build_script, (
            f"scripts/build.sh must set {option}=OFF by default so stale "
            "CMakeCache opt-ins do not leak into normal builds"
        )


def test_unavailable_opticks_disables_feature_instead_of_stub_linkage() -> None:
    cmake_text = Path("CMakeLists.txt").read_text(encoding="utf-8")

    assert "Opticks ${YELLOW}Not Found${RESET} - using CPU optical simulation" in cmake_text
    assert 'set(WITH_OPTICKS OFF CACHE BOOL "" FORCE)' in cmake_text
    assert "Stub Mode" not in cmake_text
    assert "Enabled (Stub - CPU)" not in cmake_text
