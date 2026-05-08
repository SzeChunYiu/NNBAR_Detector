import json
from pathlib import Path

import pandas as pd

from nnbar_reconstruction.cli import main
from nnbar_reconstruction.reconstruction import reconstruct_run
from nnbar_reconstruction.validation import evaluate_reconstruction_truth


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
