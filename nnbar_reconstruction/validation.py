"""Truth-supported validation metrics for reconstructed tables."""

from __future__ import annotations

from typing import Any

import pandas as pd


def _safe_ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else 0.0


def _charged_pid_truth(name: object) -> str | None:
    normalized = str(name).strip().lower()
    if normalized == "proton":
        return "proton"
    if normalized in {"pi+", "pi-", "charged_pion"}:
        return "charged_pion"
    return None


def _photon_charge_truth(name: object) -> str | None:
    normalized = str(name).strip().lower()
    charged = {
        "e+",
        "e-",
        "mu+",
        "mu-",
        "pi+",
        "pi-",
        "proton",
        "antiproton",
        "deuteron",
        "alpha",
    }
    neutral = {"gamma", "neutron", "pi0", "opticalphoton"}
    if normalized in charged:
        return "charged"
    if normalized in neutral:
        return "neutral"
    return None


def _binary_report(
    truth_positive: pd.Series,
    predicted_positive: pd.Series,
    *,
    positive_name: str,
    negative_name: str,
) -> dict[str, Any]:
    tp = int((truth_positive & predicted_positive).sum())
    fp = int((~truth_positive & predicted_positive).sum())
    tn = int((~truth_positive & ~predicted_positive).sum())
    fn = int((truth_positive & ~predicted_positive).sum())
    positive_count = int(truth_positive.sum())
    negative_count = int((~truth_positive).sum())
    n_labeled = int(len(truth_positive))
    precision = _safe_ratio(tp, tp + fp)
    recall = _safe_ratio(tp, tp + fn)
    negative_recall = _safe_ratio(tn, tn + fp)
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    negative_precision = _safe_ratio(tn, tn + fn)
    negative_f1 = (
        2.0 * negative_precision * negative_recall / (negative_precision + negative_recall)
        if negative_precision + negative_recall > 0
        else 0.0
    )
    return {
        "n_labeled": n_labeled,
        f"true_{positive_name}": positive_count,
        f"true_{negative_name}": negative_count,
        "usable": bool(positive_count > 0 and negative_count > 0),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "accuracy": _safe_ratio(tp + tn, n_labeled),
        f"{positive_name}_precision": precision,
        f"{positive_name}_recall": recall,
        f"{negative_name}_recall": negative_recall,
        "balanced_f1": 0.5 * (f1 + negative_f1),
    }


def _empty_binary_report(*, positive_name: str, negative_name: str) -> dict[str, Any]:
    return {
        "n_labeled": 0,
        f"true_{positive_name}": 0,
        f"true_{negative_name}": 0,
        "usable": False,
        "tp": 0,
        "fp": 0,
        "tn": 0,
        "fn": 0,
        "accuracy": 0.0,
        f"{positive_name}_precision": 0.0,
        f"{positive_name}_recall": 0.0,
        f"{negative_name}_recall": 0.0,
        "balanced_f1": 0.0,
    }


def evaluate_reconstruction_truth(result: dict[str, pd.DataFrame]) -> dict[str, Any]:
    """Evaluate reconstructed tables against simulation truth labels when present."""

    charged = result.get("charged", pd.DataFrame())
    photons = result.get("photons", pd.DataFrame())

    if charged.empty or not {"truth_name", "pid_guess"}.issubset(charged.columns):
        charged_report = _empty_binary_report(positive_name="proton", negative_name="pion")
    else:
        labeled_charged = charged.copy()
        labeled_charged["truth_pid"] = labeled_charged["truth_name"].map(_charged_pid_truth)
        labeled_charged = labeled_charged[labeled_charged["truth_pid"].notna()]
        if labeled_charged.empty:
            charged_report = _empty_binary_report(positive_name="proton", negative_name="pion")
        else:
            charged_report = _binary_report(
                labeled_charged["truth_pid"] == "proton",
                labeled_charged["pid_guess"] == "proton",
                positive_name="proton",
                negative_name="pion",
            )

    if photons.empty or not {"truth_name", "has_tpc_track"}.issubset(photons.columns):
        photon_report = _empty_binary_report(positive_name="charged", negative_name="neutral")
    else:
        labeled_photons = photons.copy()
        labeled_photons["truth_charge"] = labeled_photons["truth_name"].map(_photon_charge_truth)
        labeled_photons = labeled_photons[labeled_photons["truth_charge"].notna()]
        if labeled_photons.empty:
            photon_report = _empty_binary_report(positive_name="charged", negative_name="neutral")
        else:
            photon_report = _binary_report(
                labeled_photons["truth_charge"] == "charged",
                labeled_photons["has_tpc_track"].astype(bool),
                positive_name="charged",
                negative_name="neutral",
            )

    return {
        "charged_pid": charged_report,
        "photon_charged_match": photon_report,
        "overall_usable": bool(charged_report["usable"] and photon_report["usable"]),
    }
