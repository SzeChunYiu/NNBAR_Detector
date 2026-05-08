from pathlib import Path

import pandas as pd

from nnbar_reconstruction.io import discover_runs, load_run
from nnbar_reconstruction.reconstruction import (
    reconstruct_electron_pair_objects,
    find_pi0_candidates,
    reconstruct_charged_objects,
    reconstruct_photon_objects,
    reconstruct_run,
)


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def test_discover_runs_ignores_macos_sidecars(tmp_path: Path) -> None:
    _write(pd.DataFrame({"Event_ID": [1]}), tmp_path / "TPC_output_0.parquet")
    _write(pd.DataFrame({"Event_ID": [1]}), tmp_path / "PMT_output_0.parquet")
    (tmp_path / "._TPC_output_0.parquet").write_text("not parquet", encoding="utf-8")

    assert discover_runs(tmp_path)["TPC"] == [0]
    assert discover_runs(tmp_path)["PMT"] == [0]
    assert len(load_run(tmp_path, 0)["TPC"]) == 1


def test_charged_pid_uses_dedx_and_range() -> None:
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 10, "Name": "pi+", "px": 1, "py": 0, "pz": 0, "eDep": 0.2, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 10, "Name": "pi+", "px": 1, "py": 0, "pz": 0, "eDep": 0.3, "trackl": 1.0},
            {"Event_ID": 1, "Track_ID": 20, "Name": "proton", "px": 0, "py": 1, "pz": 0, "eDep": 10.0, "trackl": 1.0},
        ]
    )
    scint = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 10, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 10, "eDep": 1.0, "particle_x": 30.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 20, "eDep": 1.0, "particle_x": 0.0, "particle_y": 0.0, "particle_z": 0.0},
            {"Event_ID": 1, "Track_ID": 20, "eDep": 1.0, "particle_x": 5.0, "particle_y": 0.0, "particle_z": 0.0},
        ]
    )

    charged = reconstruct_charged_objects(tpc, scint)
    pid_by_track = dict(zip(charged["track_id"], charged["pid_guess"]))

    assert pid_by_track[10] == "charged_pion"
    assert pid_by_track[20] == "proton"


def test_pi0_candidate_uses_thesis_selection_cuts() -> None:
    lead = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 101, "Name": "gamma", "x": 100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
            {"Event_ID": 1, "Track_ID": 102, "Name": "gamma", "x": 0.0, "y": 100.0, "z": 0.0, "eDep": 100.0},
        ]
    )

    photons = reconstruct_photon_objects(lead, pd.DataFrame(), pd.DataFrame())
    pi0 = find_pi0_candidates(photons)

    assert len(photons) == 2
    assert len(pi0) == 1
    assert 140.0 < pi0.iloc[0]["mass"] < 142.0
    assert bool(pi0.iloc[0]["passes_selection"])


def test_electron_pair_candidates_use_tpc_entry_point_separation() -> None:
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "Name": "e-", "x": 0.0, "y": 0.0, "z": 0.0, "t": 1.0, "eDep": 0.2},
            {"Event_ID": 1, "Track_ID": 1, "Name": "e-", "x": 20.0, "y": 0.0, "z": 0.0, "t": 2.0, "eDep": 0.3},
            {"Event_ID": 1, "Track_ID": 2, "Name": "e+", "x": 3.0, "y": 4.0, "z": 0.0, "t": 1.2, "eDep": 0.4},
            {"Event_ID": 1, "Track_ID": 3, "Name": "pi+", "x": 20.0, "y": 0.0, "z": 0.0, "t": 1.1, "eDep": 1.0},
        ]
    )

    pairs = reconstruct_electron_pair_objects(tpc)

    assert len(pairs) == 1
    pair = pairs.iloc[0]
    assert pair["track1_id"] == 1
    assert pair["track2_id"] == 2
    assert pair["entry_separation"] == 5.0
    assert bool(pair["has_opposite_charge_truth"])


def test_reconstruct_run_writes_expected_tables(tmp_path: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "Name": "pi+",
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 0.0,
                    "eDep": 0.5,
                    "trackl": 2.0,
                }
            ]
        ),
        tmp_path / "TPC_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 11, "Name": "gamma", "x": 100.0, "y": 0.0, "z": 0.0, "eDep": 100.0},
                {"Event_ID": 1, "Track_ID": 12, "Name": "gamma", "x": 0.0, "y": 100.0, "z": 0.0, "eDep": 100.0},
            ]
        ),
        tmp_path / "LeadGlass_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Module_ID": 4, "Photon_Index": 0, "photons": 2, "KE": 0.001, "t": 12.0},
                {"Event_ID": 1, "Module_ID": 4, "Photon_Index": 1, "photons": 2, "KE": 0.001, "t": 13.0},
            ]
        ),
        tmp_path / "PMT_output_0.parquet",
    )

    result = reconstruct_run(tmp_path, 0)

    assert set(result) == {"charged", "electron_pairs", "photons", "pi0", "events"}
    assert result["events"].iloc[0]["pion_multiplicity"] == 2
    assert result["events"].iloc[0]["pmt_photons"] == 2
    assert result["events"].iloc[0]["n_pmt_hits"] == 2


def test_preliminary_selection_uses_thesis_cutflow(tmp_path: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "Name": "pi+",
                    "Origin": "CarbonPV",
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 0.0,
                    "eDep": 1.0,
                    "trackl": 1.0,
                }
            ]
        ),
        tmp_path / "TPC_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "eDep": 30.0,
                    "y": 10.0,
                    "particle_x": 0.0,
                    "particle_y": 0.0,
                    "particle_z": 0.0,
                }
            ]
        ),
        tmp_path / "Scintillator_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 11, "Name": "gamma", "x": 0.0, "y": 100.0, "z": 0.0, "eDep": 300.0},
                {"Event_ID": 1, "Track_ID": 12, "Name": "gamma", "x": 0.0, "y": 0.0, "z": 100.0, "eDep": 300.0},
                {"Event_ID": 1, "Track_ID": 13, "Name": "gamma", "x": -100.0, "y": 0.0, "z": 0.0, "eDep": 300.0},
            ]
        ),
        tmp_path / "LeadGlass_output_0.parquet",
    )

    events = reconstruct_run(tmp_path, 0)["events"]

    assert bool(events.iloc[0]["pass_scintillator_energy"])
    assert bool(events.iloc[0]["pass_tpc_foil_track"])
    assert bool(events.iloc[0]["pass_pion_count"])
    assert bool(events.iloc[0]["pass_invariant_mass"])
    assert bool(events.iloc[0]["pass_sphericity"])
    assert bool(events.iloc[0]["pass_scintillator_balance"])
    assert bool(events.iloc[0]["passes_preliminary_selection"])


def test_event_summary_includes_thesis_calorimeter_directional_variables(tmp_path: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 1, "eDep": 30.0, "x": 0.0, "y": 0.0, "z": 10.0},
                {"Event_ID": 1, "Track_ID": 2, "eDep": 50.0, "x": 10.0, "y": 0.0, "z": 0.0},
            ]
        ),
        tmp_path / "Scintillator_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 11, "Name": "gamma", "x": 0.0, "y": 0.0, "z": 20.0, "eDep": 100.0},
                {"Event_ID": 1, "Track_ID": 12, "Name": "gamma", "x": 30.0, "y": 0.0, "z": 0.0, "eDep": 200.0},
                {"Event_ID": 1, "Track_ID": 13, "Name": "gamma", "x": 0.0, "y": 25.0, "z": 0.0, "eDep": 7.0},
                {"Event_ID": 1, "Track_ID": 14, "Name": "gamma", "x": 0.0, "y": -25.0, "z": 0.0, "eDep": 11.0},
            ]
        ),
        tmp_path / "LeadGlass_output_0.parquet",
    )

    events = reconstruct_run(tmp_path, 0)["events"]
    row = events.iloc[0]

    assert row["upper_leadglass_edep"] == 7.0
    assert row["lower_leadglass_edep"] == 11.0
    assert row["scintillator_longitudinal_energy"] == 30.0
    assert row["scintillator_transverse_energy"] == 50.0
    assert row["leadglass_longitudinal_energy"] == 100.0
    assert row["leadglass_transverse_energy"] == 218.0
    assert row["calorimeter_longitudinal_energy"] == 130.0
    assert row["calorimeter_transverse_energy"] == 268.0
