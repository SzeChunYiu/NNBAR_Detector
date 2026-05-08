import json
from pathlib import Path

import pandas as pd

from nnbar_reconstruction.calibration import scan_charged_pid_thresholds
from nnbar_reconstruction.cli import main


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def _pid_scan_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "Name": "pi+", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 1.0, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 2, "Name": "pi-", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 5.0, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 3, "Name": "proton", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 9.0, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 4, "Name": "proton", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 5.0, "trackl": 1.0},
        ]
    )
    scint = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 1, "eDep": 1.0, "particle_x": 30.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 2, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 2, "eDep": 1.0, "particle_x": 30.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 3, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 3, "eDep": 1.0, "particle_x": 30.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 4, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 4, "eDep": 1.0, "particle_x": 5.0, "particle_y": 0.0, "particle_z": 0.0},
        ]
    )
    return tpc, scint


def test_pid_threshold_scan_ranks_truth_separating_config_first() -> None:
    tpc, scint = _pid_scan_fixture()

    scan = scan_charged_pid_thresholds(
        tpc,
        scint,
        proton_dedx_values=[4.0, 8.0],
        short_range_values=[3.0, 10.0],
        short_range_proton_dedx_values=[4.0, 6.0],
    )

    best = scan.iloc[0]
    assert best["proton_dedx_min"] == 8.0
    assert best["short_range_cm"] == 10.0
    assert best["short_range_proton_dedx_min"] == 4.0
    assert best["accuracy"] == 1.0
    assert best["proton_recall"] == 1.0
    assert best["pion_recall"] == 1.0
    assert bool(best["has_both_classes"])


def test_cli_scan_pid_reports_best_thresholds(tmp_path: Path, capsys) -> None:
    tpc, scint = _pid_scan_fixture()
    _write(tpc, tmp_path / "TPC_output_0.parquet")
    _write(scint, tmp_path / "Scintillator_output_0.parquet")

    exit_code = main(
        [
            "scan-pid",
            str(tmp_path),
            "--run",
            "0",
            "--proton-dedx",
            "4,8",
            "--short-range",
            "3,10",
            "--short-range-dedx",
            "4,6",
            "--top",
            "1",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["run"] == 0
    assert payload["best"]["proton_dedx_min"] == 8.0
    assert payload["best"]["short_range_cm"] == 10.0
    assert payload["best"]["short_range_proton_dedx_min"] == 4.0
    assert payload["best"]["accuracy"] == 1.0
    assert payload["calibration_usable"]


def test_cli_scan_pid_can_combine_multiple_runs(tmp_path: Path, capsys) -> None:
    tpc, scint = _pid_scan_fixture()
    _write(tpc[tpc["Track_ID"].isin([1, 2])], tmp_path / "TPC_output_0.parquet")
    _write(scint[scint["Track_ID"].isin([1, 2])], tmp_path / "Scintillator_output_0.parquet")
    _write(tpc[tpc["Track_ID"].isin([3, 4])], tmp_path / "TPC_output_1.parquet")
    _write(scint[scint["Track_ID"].isin([3, 4])], tmp_path / "Scintillator_output_1.parquet")

    exit_code = main(
        [
            "scan-pid",
            str(tmp_path),
            "--runs",
            "0,1",
            "--proton-dedx",
            "4,8",
            "--short-range",
            "3,10",
            "--short-range-dedx",
            "4,6",
            "--top",
            "1",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["runs"] == [0, 1]
    assert payload["labeled_tracks"] == 4
    assert payload["calibration_usable"]
    assert payload["best"]["true_pion"] == 2
    assert payload["best"]["true_proton"] == 2


def test_cli_scan_pid_can_discover_all_runs(tmp_path: Path, capsys) -> None:
    tpc, scint = _pid_scan_fixture()
    _write(tpc[tpc["Track_ID"].isin([1, 2])], tmp_path / "TPC_output_0.parquet")
    _write(scint[scint["Track_ID"].isin([1, 2])], tmp_path / "Scintillator_output_0.parquet")
    _write(tpc[tpc["Track_ID"].isin([3, 4])], tmp_path / "TPC_output_1.parquet")
    _write(scint[scint["Track_ID"].isin([3, 4])], tmp_path / "Scintillator_output_1.parquet")
    (tmp_path / "._TPC_output_99.parquet").write_text("sidecar", encoding="utf-8")

    exit_code = main(
        [
            "scan-pid",
            str(tmp_path),
            "--all-runs",
            "--proton-dedx",
            "4,8",
            "--short-range",
            "3,10",
            "--short-range-dedx",
            "4,6",
            "--top",
            "1",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["runs"] == [0, 1]
    assert payload["calibration_usable"]


def test_pid_threshold_scan_marks_single_class_samples_unusable() -> None:
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "Name": "pi+", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 1.0, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 2, "Name": "pi-", "px": 1.0, "py": 0.0, "pz": 0.0, "eDep": 2.0, "trackl": 1.0},
        ]
    )

    scan = scan_charged_pid_thresholds(
        tpc,
        proton_dedx_values=[4.0],
        short_range_values=[10.0],
        short_range_proton_dedx_values=[4.0],
    )

    assert scan.iloc[0]["true_proton"] == 0
    assert scan.iloc[0]["true_pion"] == 2
    assert not bool(scan.iloc[0]["has_both_classes"])
