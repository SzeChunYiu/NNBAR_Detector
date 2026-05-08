"""Thesis-aligned offline reconstruction over simulation Parquet outputs."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .io import load_run


@dataclass(frozen=True)
class ReconstructionConfig:
    """Numerical reconstruction cuts.

    Defaults encode the preliminary thesis rules where explicit values are
    available. The charged pion/proton threshold is intentionally simple until
    calibration scans provide a measured range-dependent curve.
    """

    proton_dedx_min: float = 8.0
    short_range_cm: float = 20.0
    short_range_proton_dedx_min: float = 4.5
    min_photon_energy: float = 5.0
    pi0_mass_min: float = 100.0
    pi0_mass_max: float = 180.0
    pi0_total_energy_max: float = 720.0
    pi0_scint_energy_max: float = 250.0
    pi0_leadglass_energy_max: float = 980.0
    pi0_leadglass_fraction_min: float = 0.55
    pi0_opening_angle_min_deg: float = 30.0
    selection_scintillator_energy_min: float = 20.0
    selection_scintillator_energy_max: float = 2000.0
    selection_invariant_mass_min: float = 500.0
    selection_sphericity_min: float = 0.2
    selection_upper_scintillator_max: float = 320.0
    selection_lower_scintillator_max: float = 930.0
    electron_pair_max_entry_separation_cm: float = 5.0


DEFAULT_CONFIG = ReconstructionConfig()


def _empty(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def _safe_sum(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns or df.empty:
        return 0.0
    return float(pd.to_numeric(df[column], errors="coerce").fillna(0.0).sum())


def _pmt_photons_for_event(pmt: pd.DataFrame, event_id: int) -> tuple[int, int]:
    if pmt.empty or "Event_ID" not in pmt.columns:
        return 0, 0
    event_pmt = pmt[pmt["Event_ID"] == event_id]
    if event_pmt.empty:
        return 0, 0
    n_hits = int(len(event_pmt))
    if {"Module_ID", "photons"}.issubset(event_pmt.columns):
        photons = (
            event_pmt.groupby("Module_ID")["photons"]
            .max()
            .fillna(0)
            .astype(int)
            .sum()
        )
        return int(photons), n_hits
    return n_hits, n_hits


def _unit_vector(values: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(values))
    if norm == 0.0 or not np.isfinite(norm):
        return np.zeros(3)
    return values / norm


def _weighted_direction(group: pd.DataFrame, energy_col: str = "eDep") -> np.ndarray:
    points = group[["x", "y", "z"]].to_numpy(dtype=float, copy=False)
    if points.size == 0:
        return np.zeros(3)
    weights = pd.to_numeric(group.get(energy_col, 1.0), errors="coerce").fillna(0.0).to_numpy()
    if np.sum(weights) <= 0:
        centroid = np.mean(points, axis=0)
    else:
        centroid = np.average(points, axis=0, weights=weights)
    return _unit_vector(centroid)


def _span(points: pd.DataFrame) -> float:
    if points.empty or not {"particle_x", "particle_y", "particle_z"}.issubset(points.columns):
        return 0.0
    coords = points[["particle_x", "particle_y", "particle_z"]].to_numpy(dtype=float, copy=False)
    if len(coords) < 2:
        return 0.0
    deltas = coords[:, None, :] - coords[None, :, :]
    return float(np.sqrt(np.max(np.sum(deltas * deltas, axis=2))))


def _directional_energy(df: pd.DataFrame) -> tuple[float, float]:
    """Return signed longitudinal and transverse energy about the beam z-axis."""

    required = {"x", "y", "z", "eDep"}
    if df.empty or not required.issubset(df.columns):
        return 0.0, 0.0

    coords = df[["x", "y", "z"]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    energy = pd.to_numeric(df["eDep"], errors="coerce").to_numpy(dtype=float)
    radius = np.linalg.norm(coords, axis=1)
    valid = np.isfinite(radius) & (radius > 0.0) & np.isfinite(energy)
    if not np.any(valid):
        return 0.0, 0.0

    coords = coords[valid]
    energy = energy[valid]
    radius = radius[valid]
    longitudinal = float(np.sum(energy * coords[:, 2] / radius))
    transverse_radius = np.linalg.norm(coords[:, :2], axis=1)
    transverse = float(np.sum(energy * transverse_radius / radius))
    return longitudinal, transverse


def _has_foil_origin(group: pd.DataFrame) -> bool:
    if "Origin" not in group.columns:
        return True
    origins = group["Origin"].dropna().astype(str)
    if origins.empty:
        return False
    return bool(origins.str.contains("Carbon|Target|Foil", case=False, regex=True).any())


def reconstruct_charged_objects(
    tpc: pd.DataFrame,
    scintillator: pd.DataFrame | None = None,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Reconstruct charged-track objects from TPC hits.

    The PID rule follows the thesis separation idea: protons are high dE/dx and
    short-ranged, while charged pions sit in the lower dE/dx band.
    """

    columns = [
        "event_id",
        "track_id",
        "truth_name",
        "n_tpc_hits",
        "tpc_edep",
        "tpc_path",
        "dedx",
        "scintillator_edep",
        "scintillator_range",
        "px",
        "py",
        "pz",
        "pid_guess",
        "has_foil_origin",
    ]
    if tpc is None or tpc.empty:
        return _empty(columns)

    scintillator = scintillator if scintillator is not None else pd.DataFrame()
    rows: list[dict[str, Any]] = []

    group_cols = ["Event_ID", "Track_ID"]
    for (event_id, track_id), group in tpc.groupby(group_cols, dropna=False):
        path = _safe_sum(group, "trackl")
        edep = _safe_sum(group, "eDep")
        dedx = edep / path if path > 0 else np.nan
        direction = _unit_vector(group[["px", "py", "pz"]].mean().to_numpy(dtype=float))

        scint = pd.DataFrame()
        if not scintillator.empty:
            scint = scintillator[
                (scintillator["Event_ID"] == event_id)
                & (scintillator["Track_ID"] == track_id)
            ]
        scint_edep = _safe_sum(scint, "eDep")
        scint_range = _span(scint)

        is_proton = bool(
            np.isfinite(dedx)
            and (
                dedx >= config.proton_dedx_min
                or (
                    scint_range > 0
                    and scint_range <= config.short_range_cm
                    and dedx >= config.short_range_proton_dedx_min
                )
            )
        )

        rows.append(
            {
                "event_id": int(event_id),
                "track_id": int(track_id),
                "truth_name": str(group["Name"].iloc[0]) if "Name" in group else "",
                "n_tpc_hits": int(len(group)),
                "tpc_edep": edep,
                "tpc_path": path,
                "dedx": dedx,
                "scintillator_edep": scint_edep,
                "scintillator_range": scint_range,
                "px": direction[0],
                "py": direction[1],
                "pz": direction[2],
                "pid_guess": "proton" if is_proton else "charged_pion",
                "has_foil_origin": _has_foil_origin(group),
            }
        )

    return pd.DataFrame(rows, columns=columns)


def reconstruct_photon_objects(
    leadglass: pd.DataFrame,
    scintillator: pd.DataFrame | None = None,
    tpc: pd.DataFrame | None = None,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Reconstruct neutral photon-like objects from calorimeter deposits."""

    columns = [
        "event_id",
        "object_id",
        "source_track_id",
        "truth_name",
        "leadglass_edep",
        "scintillator_edep",
        "total_energy",
        "leadglass_fraction",
        "ux",
        "uy",
        "uz",
        "has_tpc_track",
    ]
    if leadglass is None or leadglass.empty:
        return _empty(columns)

    scintillator = scintillator if scintillator is not None else pd.DataFrame()
    tpc = tpc if tpc is not None else pd.DataFrame()
    tpc_keys = set()
    if not tpc.empty:
        tpc_keys = set(zip(tpc["Event_ID"].astype(int), tpc["Track_ID"].astype(int)))

    rows: list[dict[str, Any]] = []
    for object_id, ((event_id, track_id), group) in enumerate(
        leadglass.groupby(["Event_ID", "Track_ID"], dropna=False)
    ):
        lead_edep = _safe_sum(group, "eDep")
        if lead_edep < config.min_photon_energy:
            continue
        scint = pd.DataFrame()
        if not scintillator.empty:
            scint = scintillator[
                (scintillator["Event_ID"] == event_id)
                & (scintillator["Track_ID"] == track_id)
            ]
        scint_edep = _safe_sum(scint, "eDep")
        total = lead_edep + scint_edep
        direction = _weighted_direction(group)
        has_tpc = (int(event_id), int(track_id)) in tpc_keys

        rows.append(
            {
                "event_id": int(event_id),
                "object_id": int(object_id),
                "source_track_id": int(track_id),
                "truth_name": str(group["Name"].iloc[0]) if "Name" in group else "",
                "leadglass_edep": lead_edep,
                "scintillator_edep": scint_edep,
                "total_energy": total,
                "leadglass_fraction": lead_edep / total if total > 0 else 0.0,
                "ux": direction[0],
                "uy": direction[1],
                "uz": direction[2],
                "has_tpc_track": bool(has_tpc),
            }
        )

    return pd.DataFrame(rows, columns=columns)


def _electron_charge_sign(name: str) -> int:
    normalized = name.strip().lower()
    if normalized in {"e-", "electron"}:
        return -1
    if normalized in {"e+", "positron"}:
        return 1
    return 0


def reconstruct_electron_pair_objects(
    tpc: pd.DataFrame,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Identify close TPC-entry electron-pair candidates.

    The thesis object definition tags an e+e- candidate when two TPC track entry
    points are within 5 cm of each other. Truth names are carried only as a
    validation aid; the geometric candidate rule is the reconstruction rule.
    """

    columns = [
        "event_id",
        "pair_id",
        "track1_id",
        "track2_id",
        "track1_truth_name",
        "track2_truth_name",
        "entry_separation",
        "entry_midpoint_x",
        "entry_midpoint_y",
        "entry_midpoint_z",
        "track1_tpc_edep",
        "track2_tpc_edep",
        "total_tpc_edep",
        "has_opposite_charge_truth",
    ]
    required = {"Event_ID", "Track_ID", "x", "y", "z"}
    if tpc is None or tpc.empty or not required.issubset(tpc.columns):
        return _empty(columns)

    working = tpc.copy()
    working["_input_order"] = np.arange(len(working))
    sort_cols = ["Event_ID", "Track_ID", "t", "_input_order"] if "t" in working else ["Event_ID", "Track_ID", "_input_order"]
    first_hits = (
        working.sort_values(sort_cols)
        .groupby(["Event_ID", "Track_ID"], dropna=False, as_index=False)
        .first()
    )
    edep = (
        working.groupby(["Event_ID", "Track_ID"], dropna=False)["eDep"].sum()
        if "eDep" in working
        else pd.Series(dtype=float)
    )

    rows: list[dict[str, Any]] = []
    pair_id = 0
    for event_id, event_tracks in first_hits.groupby("Event_ID", dropna=False):
        if len(event_tracks) < 2:
            continue
        coords = event_tracks[["x", "y", "z"]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(coords).all(axis=1)
        if np.count_nonzero(valid) < 2:
            continue
        event_tracks = event_tracks.loc[valid].reset_index(drop=True)
        coords = coords[valid]

        left, right = np.triu_indices(len(event_tracks), k=1)
        distances = np.linalg.norm(coords[left] - coords[right], axis=1)
        close = distances <= config.electron_pair_max_entry_separation_cm
        for i, j, distance in zip(left[close], right[close], distances[close]):
            track1 = event_tracks.iloc[int(i)]
            track2 = event_tracks.iloc[int(j)]
            track1_id = int(track1.Track_ID)
            track2_id = int(track2.Track_ID)
            track1_name = str(track1.Name) if "Name" in event_tracks else ""
            track2_name = str(track2.Name) if "Name" in event_tracks else ""
            track1_edep = float(edep.get((event_id, track1.Track_ID), 0.0)) if not edep.empty else 0.0
            track2_edep = float(edep.get((event_id, track2.Track_ID), 0.0)) if not edep.empty else 0.0
            midpoint = 0.5 * (coords[int(i)] + coords[int(j)])
            rows.append(
                {
                    "event_id": int(event_id),
                    "pair_id": pair_id,
                    "track1_id": track1_id,
                    "track2_id": track2_id,
                    "track1_truth_name": track1_name,
                    "track2_truth_name": track2_name,
                    "entry_separation": float(distance),
                    "entry_midpoint_x": float(midpoint[0]),
                    "entry_midpoint_y": float(midpoint[1]),
                    "entry_midpoint_z": float(midpoint[2]),
                    "track1_tpc_edep": track1_edep,
                    "track2_tpc_edep": track2_edep,
                    "total_tpc_edep": track1_edep + track2_edep,
                    "has_opposite_charge_truth": bool(
                        _electron_charge_sign(track1_name) * _electron_charge_sign(track2_name) == -1
                    ),
                }
            )
            pair_id += 1

    return pd.DataFrame(rows, columns=columns)


def find_pi0_candidates(
    photons: pd.DataFrame,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Build neutral-pion candidates from photon-like object pairs."""

    columns = [
        "event_id",
        "photon1_id",
        "photon2_id",
        "mass",
        "opening_angle_deg",
        "total_energy",
        "leadglass_edep",
        "scintillator_edep",
        "leadglass_fraction",
        "passes_selection",
    ]
    if photons is None or photons.empty:
        return _empty(columns)

    rows: list[dict[str, Any]] = []
    for event_id, group in photons.groupby("event_id", dropna=False):
        neutral = group[group["has_tpc_track"] == False]  # noqa: E712
        for _, a in neutral.iterrows():
            for _, b in neutral[neutral["object_id"] > a["object_id"]].iterrows():
                va = np.array([a.ux, a.uy, a.uz], dtype=float)
                vb = np.array([b.ux, b.uy, b.uz], dtype=float)
                cosang = float(np.clip(np.dot(va, vb), -1.0, 1.0))
                angle = float(np.degrees(np.arccos(cosang)))
                mass2 = 2.0 * float(a.total_energy) * float(b.total_energy) * (1.0 - cosang)
                mass = float(np.sqrt(max(mass2, 0.0)))
                lead = float(a.leadglass_edep + b.leadglass_edep)
                scint = float(a.scintillator_edep + b.scintillator_edep)
                total = lead + scint
                fraction = lead / total if total > 0 else 0.0
                passes = (
                    config.pi0_mass_min <= mass <= config.pi0_mass_max
                    and total <= config.pi0_total_energy_max
                    and scint <= config.pi0_scint_energy_max
                    and lead <= config.pi0_leadglass_energy_max
                    and fraction >= config.pi0_leadglass_fraction_min
                    and angle >= config.pi0_opening_angle_min_deg
                )
                rows.append(
                    {
                        "event_id": int(event_id),
                        "photon1_id": int(a.object_id),
                        "photon2_id": int(b.object_id),
                        "mass": mass,
                        "opening_angle_deg": angle,
                        "total_energy": total,
                        "leadglass_edep": lead,
                        "scintillator_edep": scint,
                        "leadglass_fraction": fraction,
                        "passes_selection": bool(passes),
                    }
                )

    return pd.DataFrame(rows, columns=columns)


def _sphericity(vectors: np.ndarray) -> float:
    if vectors.size == 0:
        return float("nan")
    denom = float(np.sum(vectors * vectors))
    if denom <= 0:
        return float("nan")
    tensor = np.zeros((3, 3), dtype=float)
    for vec in vectors:
        tensor += np.outer(vec, vec)
    eigvals = np.linalg.eigvalsh(tensor / denom)
    eigvals.sort()
    return float(1.5 * (eigvals[0] + eigvals[1]))


def _visible_invariant_mass(charged: pd.DataFrame, photons: pd.DataFrame) -> float:
    total_energy = 0.0
    momentum = np.zeros(3, dtype=float)

    masses = {"charged_pion": 139.57039, "proton": 938.272088}
    for _, row in charged.iterrows():
        kinetic = max(float(row.tpc_edep + row.scintillator_edep), 0.0)
        mass = masses.get(str(row.pid_guess), 139.57039)
        energy = kinetic + mass
        p_mag = np.sqrt(max(energy * energy - mass * mass, 0.0))
        direction = _unit_vector(np.array([row.px, row.py, row.pz], dtype=float))
        total_energy += energy
        momentum += p_mag * direction

    for _, row in photons.iterrows():
        energy = max(float(row.total_energy), 0.0)
        direction = _unit_vector(np.array([row.ux, row.uy, row.uz], dtype=float))
        total_energy += energy
        momentum += energy * direction

    mass2 = total_energy * total_energy - float(np.dot(momentum, momentum))
    return float(np.sqrt(max(mass2, 0.0))) if total_energy > 0 else float("nan")


def _selection_flags(row: dict[str, Any], config: ReconstructionConfig) -> dict[str, bool]:
    pass_scint = (
        config.selection_scintillator_energy_min
        <= row["scintillator_edep"]
        <= config.selection_scintillator_energy_max
    )
    pass_tpc = bool(row["has_foil_tpc_track"])
    pass_pion = row["pion_multiplicity"] >= 1
    mass = row["visible_invariant_mass"]
    pass_mass = bool(np.isfinite(mass) and mass >= config.selection_invariant_mass_min)
    sphericity = row["sphericity"]
    pass_sphericity = bool(np.isfinite(sphericity) and sphericity >= config.selection_sphericity_min)
    pass_scint_balance = (
        row["upper_scintillator_edep"] <= config.selection_upper_scintillator_max
        and row["lower_scintillator_edep"] <= config.selection_lower_scintillator_max
    )
    return {
        "pass_scintillator_energy": pass_scint,
        "pass_tpc_foil_track": pass_tpc,
        "pass_pion_count": pass_pion,
        "pass_invariant_mass": pass_mass,
        "pass_sphericity": pass_sphericity,
        "pass_scintillator_balance": pass_scint_balance,
        "passes_preliminary_selection": all(
            [pass_scint, pass_tpc, pass_pion, pass_mass, pass_sphericity, pass_scint_balance]
        ),
    }


def summarize_events(
    data: dict[str, pd.DataFrame],
    charged: pd.DataFrame,
    photons: pd.DataFrame,
    pi0: pd.DataFrame,
    electron_pairs: pd.DataFrame | None = None,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Create per-event variables used by the preliminary selection."""

    event_ids: set[int] = set()
    for table in data.values():
        if table is not None and not table.empty and "Event_ID" in table:
            event_ids.update(table["Event_ID"].dropna().astype(int).tolist())
    event_ids.update(charged["event_id"].dropna().astype(int).tolist() if not charged.empty else [])
    event_ids.update(photons["event_id"].dropna().astype(int).tolist() if not photons.empty else [])
    if electron_pairs is not None and not electron_pairs.empty:
        event_ids.update(electron_pairs["event_id"].dropna().astype(int).tolist())

    rows: list[dict[str, Any]] = []
    for event_id in sorted(event_ids):
        scint = data.get("Scintillator", pd.DataFrame())
        lead = data.get("LeadGlass", pd.DataFrame())
        tpc = data.get("TPC", pd.DataFrame())
        pmt = data.get("PMT", pd.DataFrame())
        scint_event = scint[scint["Event_ID"] == event_id] if not scint.empty and "Event_ID" in scint else pd.DataFrame()
        lead_event = lead[lead["Event_ID"] == event_id] if not lead.empty and "Event_ID" in lead else pd.DataFrame()
        tpc_event = tpc[tpc["Event_ID"] == event_id] if not tpc.empty and "Event_ID" in tpc else pd.DataFrame()
        scint_e = _safe_sum(scint_event, "eDep")
        upper_scint_e = _safe_sum(scint_event[scint_event["y"] > 0], "eDep") if "y" in scint_event else 0.0
        lower_scint_e = _safe_sum(scint_event[scint_event["y"] < 0], "eDep") if "y" in scint_event else 0.0
        lead_e = _safe_sum(lead_event, "eDep")
        upper_lead_e = _safe_sum(lead_event[lead_event["y"] > 0], "eDep") if "y" in lead_event else 0.0
        lower_lead_e = _safe_sum(lead_event[lead_event["y"] < 0], "eDep") if "y" in lead_event else 0.0
        scint_longitudinal_e, scint_transverse_e = _directional_energy(scint_event)
        lead_longitudinal_e, lead_transverse_e = _directional_energy(lead_event)
        tpc_e = _safe_sum(tpc_event, "eDep")
        pmt_photons, pmt_hits = _pmt_photons_for_event(pmt, event_id)

        cands = charged[charged["event_id"] == event_id] if not charged.empty else charged
        gammas = photons[photons["event_id"] == event_id] if not photons.empty else photons
        pi0s = pi0[(pi0["event_id"] == event_id) & (pi0["passes_selection"] == True)] if not pi0.empty else pi0  # noqa: E712
        pairs = (
            electron_pairs[electron_pairs["event_id"] == event_id]
            if electron_pairs is not None and not electron_pairs.empty
            else pd.DataFrame()
        )
        charged_pion_count = int((cands["pid_guess"] == "charged_pion").sum()) if not cands.empty else 0

        vectors: list[np.ndarray] = []
        for _, row in cands.iterrows():
            scale = max(float(row.tpc_edep + row.scintillator_edep), 1.0)
            vectors.append(scale * np.array([row.px, row.py, row.pz], dtype=float))
        for _, row in gammas.iterrows():
            scale = max(float(row.total_energy), 0.0)
            vectors.append(scale * np.array([row.ux, row.uy, row.uz], dtype=float))

        event_row = {
            "event_id": event_id,
            "tpc_edep": tpc_e,
            "scintillator_edep": scint_e,
            "upper_scintillator_edep": upper_scint_e,
            "lower_scintillator_edep": lower_scint_e,
            "leadglass_edep": lead_e,
            "upper_leadglass_edep": upper_lead_e,
            "lower_leadglass_edep": lower_lead_e,
            "calorimeter_edep": scint_e + lead_e,
            "scintillator_longitudinal_energy": scint_longitudinal_e,
            "scintillator_transverse_energy": scint_transverse_e,
            "leadglass_longitudinal_energy": lead_longitudinal_e,
            "leadglass_transverse_energy": lead_transverse_e,
            "calorimeter_longitudinal_energy": scint_longitudinal_e + lead_longitudinal_e,
            "calorimeter_transverse_energy": scint_transverse_e + lead_transverse_e,
            "pmt_photons": pmt_photons,
            "n_pmt_hits": pmt_hits,
            "n_charged_objects": int(len(cands)),
            "n_charged_pion": charged_pion_count,
            "n_proton": int((cands["pid_guess"] == "proton").sum()) if not cands.empty else 0,
            "n_photon_like": int(len(gammas)),
            "n_pi0": int(len(pi0s)),
            "n_electron_pairs": int(len(pairs)),
            "pion_multiplicity": charged_pion_count + int(len(pi0s)),
            "has_tpc_track": int(len(cands)) > 0,
            "has_foil_tpc_track": bool(cands["has_foil_origin"].any()) if not cands.empty else False,
            "visible_invariant_mass": _visible_invariant_mass(cands, gammas),
            "sphericity": _sphericity(np.vstack(vectors)) if vectors else float("nan"),
        }
        event_row.update(_selection_flags(event_row, config))
        rows.append(event_row)

    return pd.DataFrame(rows)


def reconstruct_run(
    output_dir: str | Path,
    run: int = 0,
    config: ReconstructionConfig = DEFAULT_CONFIG,
) -> dict[str, pd.DataFrame]:
    """Load and reconstruct a simulation run."""

    data = load_run(output_dir, run)
    charged = reconstruct_charged_objects(data["TPC"], data["Scintillator"], config)
    electron_pairs = reconstruct_electron_pair_objects(data["TPC"], config)
    photons = reconstruct_photon_objects(data["LeadGlass"], data["Scintillator"], data["TPC"], config)
    pi0 = find_pi0_candidates(photons, config)
    events = summarize_events(data, charged, photons, pi0, electron_pairs, config)
    return {
        "charged": charged,
        "electron_pairs": electron_pairs,
        "photons": photons,
        "pi0": pi0,
        "events": events,
    }
