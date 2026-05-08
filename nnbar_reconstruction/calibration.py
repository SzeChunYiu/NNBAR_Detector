"""Calibration helpers for simulation-truth reconstruction studies."""

from __future__ import annotations

from dataclasses import replace
from itertools import product
from typing import Iterable

import pandas as pd

from .reconstruction import DEFAULT_CONFIG, ReconstructionConfig, reconstruct_charged_objects


def _truth_pid(name: object) -> str | None:
    normalized = str(name).strip().lower()
    if normalized == "proton":
        return "proton"
    if normalized in {"pi+", "pi-", "charged_pion"}:
        return "charged_pion"
    return None


def _safe_ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else 0.0


def scan_charged_pid_thresholds(
    tpc: pd.DataFrame,
    scintillator: pd.DataFrame | None = None,
    *,
    proton_dedx_values: Iterable[float],
    short_range_values: Iterable[float],
    short_range_proton_dedx_values: Iterable[float],
    base_config: ReconstructionConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Score charged pion/proton PID thresholds against simulation truth labels."""

    columns = [
        "proton_dedx_min",
        "short_range_cm",
        "short_range_proton_dedx_min",
        "n_labeled_tracks",
        "true_proton",
        "true_pion",
        "has_both_classes",
        "tp",
        "fp",
        "tn",
        "fn",
        "accuracy",
        "proton_precision",
        "proton_recall",
        "pion_recall",
        "balanced_f1",
    ]
    if tpc is None or tpc.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, float | int]] = []
    for proton_dedx, short_range, short_range_dedx in product(
        proton_dedx_values,
        short_range_values,
        short_range_proton_dedx_values,
    ):
        config = replace(
            base_config,
            proton_dedx_min=float(proton_dedx),
            short_range_cm=float(short_range),
            short_range_proton_dedx_min=float(short_range_dedx),
        )
        charged = reconstruct_charged_objects(tpc, scintillator, config)
        if charged.empty:
            continue

        labeled = charged.copy()
        labeled["truth_pid"] = labeled["truth_name"].map(_truth_pid)
        labeled = labeled[labeled["truth_pid"].notna()]
        if labeled.empty:
            continue

        truth_proton = labeled["truth_pid"] == "proton"
        predicted_proton = labeled["pid_guess"] == "proton"
        tp = int((truth_proton & predicted_proton).sum())
        fp = int((~truth_proton & predicted_proton).sum())
        tn = int((~truth_proton & ~predicted_proton).sum())
        fn = int((truth_proton & ~predicted_proton).sum())
        true_proton = int(truth_proton.sum())
        true_pion = int((~truth_proton).sum())

        proton_precision = _safe_ratio(tp, tp + fp)
        proton_recall = _safe_ratio(tp, tp + fn)
        pion_recall = _safe_ratio(tn, tn + fp)
        proton_f1 = (
            2.0 * proton_precision * proton_recall / (proton_precision + proton_recall)
            if proton_precision + proton_recall > 0.0
            else 0.0
        )
        pion_precision = _safe_ratio(tn, tn + fn)
        pion_f1 = (
            2.0 * pion_precision * pion_recall / (pion_precision + pion_recall)
            if pion_precision + pion_recall > 0.0
            else 0.0
        )
        n_labeled = int(len(labeled))

        rows.append(
            {
                "proton_dedx_min": float(proton_dedx),
                "short_range_cm": float(short_range),
                "short_range_proton_dedx_min": float(short_range_dedx),
                "n_labeled_tracks": n_labeled,
                "true_proton": true_proton,
                "true_pion": true_pion,
                "has_both_classes": bool(true_proton > 0 and true_pion > 0),
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "accuracy": _safe_ratio(tp + tn, n_labeled),
                "proton_precision": proton_precision,
                "proton_recall": proton_recall,
                "pion_recall": pion_recall,
                "balanced_f1": 0.5 * (proton_f1 + pion_f1),
            }
        )

    result = pd.DataFrame(rows, columns=columns)
    if result.empty:
        return result
    return result.sort_values(
        ["has_both_classes", "balanced_f1", "accuracy", "proton_recall", "pion_recall"],
        ascending=[False, False, False, False, False],
        kind="mergesort",
    ).reset_index(drop=True)
