import json
from pathlib import Path

import pandas as pd

from nnbar_reconstruction.cli import main
from nnbar_reconstruction.reconstruction import reconstruct_run
from nnbar_reconstruction.validation import assess_validation_readiness, evaluate_reconstruction_truth


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def _validation_fixture(output_dir: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 10,
                    "Name": "pi+",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 1.0,
                    "eDep": 0.5,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 10,
                    "Name": "pi+",
                    "Origin": "CarbonPV",
                    "x": 50.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 2.0,
                    "eDep": 0.5,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 20,
                    "Name": "proton",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 1.0,
                    "eDep": 9.0,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 20,
                    "Name": "proton",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 50.0,
                    "z": 0.0,
                    "t": 2.0,
                    "eDep": 9.0,
                    "trackl": 1.0,
                },
            ]
        ),
        output_dir / "TPC_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 101, "Name": "pi+", "x": 100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
                {"Event_ID": 1, "Track_ID": 102, "Name": "gamma", "x": -100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
            ]
        ),
        output_dir / "LeadGlass_output_0.parquet",
    )


def _single_class_validation_fixture(output_dir: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 10,
                    "Name": "pi+",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 1.0,
                    "eDep": 0.5,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 10,
                    "Name": "pi+",
                    "Origin": "CarbonPV",
                    "x": 50.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 2.0,
                    "eDep": 0.5,
                    "trackl": 1.0,
                },
            ]
        ),
        output_dir / "TPC_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 101, "Name": "pi+", "x": 100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
            ]
        ),
        output_dir / "LeadGlass_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 20,
                    "Name": "proton",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 0.0,
                    "z": 0.0,
                    "t": 1.0,
                    "eDep": 9.0,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 20,
                    "Name": "proton",
                    "Origin": "CarbonPV",
                    "x": 0.0,
                    "y": 50.0,
                    "z": 0.0,
                    "t": 2.0,
                    "eDep": 9.0,
                    "trackl": 1.0,
                },
            ]
        ),
        output_dir / "TPC_output_1.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 102, "Name": "gamma", "x": -100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
            ]
        ),
        output_dir / "LeadGlass_output_1.parquet",
    )


def test_reconstruction_truth_validation_reports_supported_metrics(tmp_path: Path) -> None:
    _validation_fixture(tmp_path)

    result = reconstruct_run(tmp_path, 0)
    report = evaluate_reconstruction_truth(result)

    assert report["charged_pid"]["usable"]
    assert report["charged_pid"]["accuracy"] == 1.0
    assert report["charged_pid"]["true_proton"] == 1
    assert report["charged_pid"]["true_pion"] == 1
    assert report["photon_charged_match"]["usable"]
    assert report["photon_charged_match"]["accuracy"] == 1.0
    assert report["photon_charged_match"]["true_charged"] == 1
    assert report["photon_charged_match"]["true_neutral"] == 1
    assert report["overall_usable"]


def test_cli_validate_reco_emits_json_report(tmp_path: Path, capsys) -> None:
    _validation_fixture(tmp_path)

    exit_code = main(["validate-reco", str(tmp_path), "--run", "0"])

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["run"] == 0
    assert payload["runs"] == [0]
    assert payload["overall_usable"]
    assert payload["charged_pid"]["accuracy"] == 1.0
    assert payload["photon_charged_match"]["accuracy"] == 1.0


def test_cli_validate_reco_aggregates_class_support_across_runs(tmp_path: Path, capsys) -> None:
    _single_class_validation_fixture(tmp_path)

    exit_code = main(["validate-reco", str(tmp_path), "--runs", "0,1"])

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["run"] is None
    assert payload["runs"] == [0, 1]
    assert not payload["run_reports"][0]["overall_usable"]
    assert not payload["run_reports"][1]["overall_usable"]
    assert payload["aggregate"]["overall_usable"]
    assert payload["aggregate"]["charged_pid"]["true_pion"] == 1
    assert payload["aggregate"]["charged_pid"]["true_proton"] == 1
    assert payload["aggregate"]["photon_charged_match"]["true_charged"] == 1
    assert payload["aggregate"]["photon_charged_match"]["true_neutral"] == 1


def test_cli_validate_reco_can_discover_all_runs(tmp_path: Path, capsys) -> None:
    _single_class_validation_fixture(tmp_path)
    (tmp_path / "._LeadGlass_output_99.parquet").write_text("sidecar", encoding="utf-8")

    exit_code = main(["validate-reco", str(tmp_path), "--all-runs"])

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["runs"] == [0, 1]
    assert payload["aggregate"]["overall_usable"]


def test_validation_readiness_gate_reports_missing_class_counts(tmp_path: Path) -> None:
    _single_class_validation_fixture(tmp_path)

    report = evaluate_reconstruction_truth(reconstruct_run(tmp_path, 0))
    readiness = assess_validation_readiness(
        report,
        min_class_count=2,
        min_accuracy=0.9,
        min_balanced_f1=0.9,
    )

    assert not readiness["passed"]
    assert "charged_pid.true_proton 0 < 2" in readiness["failed_requirements"]
    assert "photon_charged_match.true_neutral 0 < 2" in readiness["failed_requirements"]


def test_cli_validate_reco_emits_thresholded_readiness(tmp_path: Path, capsys) -> None:
    _validation_fixture(tmp_path)

    exit_code = main(
        [
            "validate-reco",
            str(tmp_path),
            "--run",
            "0",
            "--min-class-count",
            "1",
            "--min-accuracy",
            "0.9",
            "--min-balanced-f1",
            "0.9",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["readiness"]["passed"]
    assert payload["readiness"]["failed_requirements"] == []
