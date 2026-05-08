from pathlib import Path

import pandas as pd
import pytest

from nnbar_reconstruction.io import discover_runs, load_run
from nnbar_reconstruction.reconstruction import (
    ReconstructionConfig,
    annotate_timing_windows,
    reconstruct_electron_pair_objects,
    reconstruct_event_vertices,
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


def test_charged_direction_uses_tpc_hit_geometry_before_momentum_columns() -> None:
    tpc = pd.DataFrame(
        [
            {
                "Event_ID": 1,
                "Track_ID": 10,
                "Name": "pi+",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "px": 1.0,
                "py": 0.0,
                "pz": 0.0,
                "t": 1.0,
                "eDep": 0.2,
                "trackl": 1.0,
            },
            {
                "Event_ID": 1,
                "Track_ID": 10,
                "Name": "pi+",
                "x": 0.0,
                "y": 10.0,
                "z": 0.0,
                "px": 1.0,
                "py": 0.0,
                "pz": 0.0,
                "t": 2.0,
                "eDep": 0.3,
                "trackl": 1.0,
            },
        ]
    )

    charged = reconstruct_charged_objects(tpc)
    row = charged.iloc[0]

    assert row["px"] == pytest.approx(0.0)
    assert row["py"] == pytest.approx(1.0)
    assert row["pz"] == pytest.approx(0.0)


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


def test_event_vertex_projects_tpc_tracks_back_to_source_foil() -> None:
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "x": 10.0, "y": 0.0, "z": 10.0, "t": 1.0},
            {"Event_ID": 1, "Track_ID": 1, "x": 20.0, "y": 0.0, "z": 20.0, "t": 2.0},
            {"Event_ID": 1, "Track_ID": 2, "x": 0.0, "y": 20.0, "z": 20.0, "t": 1.0},
            {"Event_ID": 1, "Track_ID": 2, "x": 0.0, "y": 30.0, "z": 30.0, "t": 2.0},
            {"Event_ID": 1, "Track_ID": 3, "x": 9.0, "y": 9.0, "z": 10.0, "t": 1.0},
            {"Event_ID": 1, "Track_ID": 3, "x": 19.0, "y": 9.0, "z": 10.0, "t": 2.0},
        ]
    )

    vertices = reconstruct_event_vertices(tpc)

    assert len(vertices) == 1
    vertex = vertices.iloc[0]
    assert vertex["event_id"] == 1
    assert vertex["n_projected_tracks"] == 2
    assert vertex["vertex_x"] == 0.0
    assert vertex["vertex_y"] == 0.0
    assert vertex["vertex_z"] == 0.0
    assert vertex["vertex_radial_spread"] == 0.0
    assert vertex["n_skipped_tracks"] == 1


def test_timing_windows_filter_scintillator_and_leadglass_hits_from_vertex(tmp_path: Path) -> None:
    config = ReconstructionConfig(
        scintillator_time_resolution_ns=0.05,
        leadglass_time_resolution_ns=0.05,
    )
    vertices = pd.DataFrame(
        [
            {
                "event_id": 1,
                "vertex_x": 0.0,
                "vertex_y": 0.0,
                "vertex_z": 0.0,
                "vertex_time_ns": 0.0,
            }
        ]
    )
    scint = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 1, "x": 0.0, "y": 0.0, "z": 29.9792458, "t": 1.10, "eDep": 30.0},
            {"Event_ID": 1, "Track_ID": 2, "x": 0.0, "y": 0.0, "z": 29.9792458, "t": 1.50, "eDep": 70.0},
        ]
    )
    lead = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 11, "Name": "gamma", "x": 0.0, "y": 0.0, "z": 29.9792458, "t": 1.00, "eDep": 100.0},
            {"Event_ID": 1, "Track_ID": 12, "Name": "gamma", "x": 0.0, "y": 0.0, "z": 29.9792458, "t": 1.25, "eDep": 200.0},
        ]
    )

    annotated_scint = annotate_timing_windows(scint, vertices, "scintillator", config)
    annotated_lead = annotate_timing_windows(lead, vertices, "leadglass", config)

    assert annotated_scint["in_timing_window"].tolist() == [True, False]
    assert annotated_lead["in_timing_window"].tolist() == [True, False]

    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 1, "x": 10.0, "y": 0.0, "z": 10.0, "t": 1.0},
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "Name": "pi+",
                    "x": 20.0,
                    "y": 0.0,
                    "z": 20.0,
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 0.0,
                    "t": 2.0,
                    "eDep": 0.1,
                    "trackl": 1.0,
                },
            ]
        ),
        tmp_path / "TPC_output_0.parquet",
    )
    _write(scint, tmp_path / "Scintillator_output_0.parquet")
    _write(lead, tmp_path / "LeadGlass_output_0.parquet")

    events = reconstruct_run(tmp_path, 0, config=config)["events"]
    row = events.iloc[0]

    assert row["vertex_time_ns"] == 0.0
    assert row["scintillator_timing_edep"] == 30.0
    assert row["scintillator_out_of_time_edep"] == 70.0
    assert row["leadglass_timing_edep"] == 100.0
    assert row["leadglass_out_of_time_edep"] == 200.0
    assert row["calorimeter_timing_edep"] == 130.0


def test_photon_directions_use_reconstructed_vertex_for_pi0_mass(tmp_path: Path) -> None:
    _write(
        pd.DataFrame(
            [
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "x": 20.0,
                    "y": 10.0,
                    "z": 10.0,
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 1.0,
                    "t": 1.0,
                    "eDep": 0.1,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 1,
                    "x": 30.0,
                    "y": 10.0,
                    "z": 20.0,
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 1.0,
                    "t": 2.0,
                    "eDep": 0.1,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 2,
                    "x": 10.0,
                    "y": 20.0,
                    "z": 10.0,
                    "px": 0.0,
                    "py": 1.0,
                    "pz": 1.0,
                    "t": 1.0,
                    "eDep": 0.1,
                    "trackl": 1.0,
                },
                {
                    "Event_ID": 1,
                    "Track_ID": 2,
                    "x": 10.0,
                    "y": 30.0,
                    "z": 20.0,
                    "px": 0.0,
                    "py": 1.0,
                    "pz": 1.0,
                    "t": 2.0,
                    "eDep": 0.1,
                    "trackl": 1.0,
                },
            ]
        ),
        tmp_path / "TPC_output_0.parquet",
    )
    _write(
        pd.DataFrame(
            [
                {"Event_ID": 1, "Track_ID": 101, "Name": "gamma", "x": 110.0, "y": 10.0, "z": 0.0, "eDep": 100.0},
                {"Event_ID": 1, "Track_ID": 102, "Name": "gamma", "x": 10.0, "y": 110.0, "z": 0.0, "eDep": 100.0},
            ]
        ),
        tmp_path / "LeadGlass_output_0.parquet",
    )

    result = reconstruct_run(tmp_path, 0)
    photons = result["photons"].sort_values("source_track_id").reset_index(drop=True)
    pi0 = result["pi0"].iloc[0]

    assert photons.loc[0, "vertex_x"] == pytest.approx(10.0)
    assert photons.loc[0, "vertex_y"] == pytest.approx(10.0)
    assert photons.loc[0, "ux"] == pytest.approx(1.0)
    assert photons.loc[0, "uy"] == pytest.approx(0.0)
    assert photons.loc[1, "ux"] == pytest.approx(0.0)
    assert photons.loc[1, "uy"] == pytest.approx(1.0)
    assert pi0["opening_angle_deg"] == pytest.approx(90.0)
    assert pi0["mass"] == pytest.approx(141.4213562373095)


def test_photon_charged_match_uses_geometry_not_truth_track_id() -> None:
    config = ReconstructionConfig(charged_cluster_match_angle_deg=5.0)
    vertices = pd.DataFrame(
        [
            {
                "event_id": 1,
                "vertex_x": 0.0,
                "vertex_y": 0.0,
                "vertex_z": 0.0,
            }
        ]
    )
    tpc = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 7, "x": 10.0, "y": 0.0, "z": 0.0, "t": 1.0},
            {"Event_ID": 1, "Track_ID": 7, "x": 50.0, "y": 0.0, "z": 0.0, "t": 2.0},
        ]
    )
    lead = pd.DataFrame(
        [
            {"Event_ID": 1, "Track_ID": 101, "Name": "pi+", "x": 100.0, "y": 1.0, "z": 0.0, "eDep": 100.0},
            {"Event_ID": 1, "Track_ID": 7, "Name": "gamma", "x": 0.0, "y": 100.0, "z": 0.0, "eDep": 100.0},
        ]
    )

    photons = reconstruct_photon_objects(lead, tpc=tpc, config=config, vertices=vertices)
    by_source = photons.set_index("source_track_id")

    assert bool(by_source.loc[101, "has_tpc_track"])
    assert by_source.loc[101, "matched_tpc_track_id"] == 7
    assert by_source.loc[101, "charged_match_angle_deg"] < 1.0
    assert not bool(by_source.loc[7, "has_tpc_track"])
    assert pd.isna(by_source.loc[7, "matched_tpc_track_id"])


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

    assert set(result) == {"charged", "electron_pairs", "vertices", "photons", "pi0", "events"}
    assert result["events"].iloc[0]["pion_multiplicity"] == 2
    assert result["events"].iloc[0]["pmt_photons"] == 2
    assert result["events"].iloc[0]["n_pmt_hits"] == 2
    assert "vertex_x" in result["events"].columns


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
