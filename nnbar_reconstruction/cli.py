"""Command-line interface for NNBAR offline reconstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .calibration import scan_charged_pid_thresholds
from .io import load_run
from .reconstruction import reconstruct_run


def _write_tables(result: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in result.items():
        table.to_csv(output_dir / f"{name}.csv", index=False)


def _cutflow(events) -> dict[str, int]:
    cuts = [
        "pass_scintillator_energy",
        "pass_tpc_foil_track",
        "pass_pion_count",
        "pass_invariant_mass",
        "pass_sphericity",
        "pass_scintillator_balance",
    ]
    if events.empty:
        return {cut: 0 for cut in cuts}

    mask = None
    counts: dict[str, int] = {}
    for cut in cuts:
        current = events[cut].astype(bool) if cut in events else False
        mask = current if mask is None else (mask & current)
        counts[cut] = int(mask.sum())
    return counts


def summarize(args: argparse.Namespace) -> int:
    result = reconstruct_run(args.output_dir, run=args.run)
    if args.tables_dir:
        _write_tables(result, Path(args.tables_dir))

    events = result["events"]
    summary = {
        "run": args.run,
        "events": int(len(events)),
        "charged_objects": int(len(result["charged"])),
        "electron_pair_candidates": int(len(result["electron_pairs"])),
        "reconstructed_vertices": int(len(result["vertices"])),
        "photon_like_objects": int(len(result["photons"])),
        "pi0_candidates": int(len(result["pi0"])),
        "selected_pi0": int(result["pi0"]["passes_selection"].sum()) if not result["pi0"].empty else 0,
        "total_calorimeter_edep": float(events["calorimeter_edep"].sum()) if not events.empty else 0.0,
        "timing_calorimeter_edep": (
            float(events["calorimeter_timing_edep"].sum())
            if not events.empty and "calorimeter_timing_edep" in events
            else 0.0
        ),
        "out_of_time_calorimeter_edep": (
            float(events["calorimeter_out_of_time_edep"].sum())
            if not events.empty and "calorimeter_out_of_time_edep" in events
            else 0.0
        ),
        "pmt_photons": int(events["pmt_photons"].sum()) if not events.empty and "pmt_photons" in events else 0,
        "pmt_hits": int(events["n_pmt_hits"].sum()) if not events.empty and "n_pmt_hits" in events else 0,
        "selected_events": int(events["passes_preliminary_selection"].sum()) if not events.empty else 0,
        "cutflow": _cutflow(events),
        "mean_sphericity": float(events["sphericity"].dropna().mean()) if not events.empty else None,
    }

    payload = json.dumps(summary, indent=2, sort_keys=True)
    if args.json:
        Path(args.json).write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0


def _float_grid(value: str) -> list[float]:
    try:
        values = [float(part.strip()) for part in value.split(",") if part.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid float grid: {value}") from exc
    if not values:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return values


def _int_grid(value: str) -> list[int]:
    try:
        values = [int(part.strip()) for part in value.split(",") if part.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid integer grid: {value}") from exc
    if not values:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return values


def _load_pid_scan_tables(output_dir: Path, runs: list[int]) -> tuple:
    tpc_tables = []
    scintillator_tables = []
    for run in runs:
        data = load_run(output_dir, run)
        event_offset = int(run) * 1_000_000_000
        tpc = data["TPC"].copy()
        scintillator = data["Scintillator"].copy()
        if not tpc.empty and "Event_ID" in tpc:
            tpc["Event_ID"] = tpc["Event_ID"].astype(int) + event_offset
        if not scintillator.empty and "Event_ID" in scintillator:
            scintillator["Event_ID"] = scintillator["Event_ID"].astype(int) + event_offset
        tpc_tables.append(tpc)
        scintillator_tables.append(scintillator)
    return (
        pd.concat(tpc_tables, ignore_index=True) if tpc_tables else pd.DataFrame(),
        pd.concat(scintillator_tables, ignore_index=True) if scintillator_tables else pd.DataFrame(),
    )


def scan_pid(args: argparse.Namespace) -> int:
    runs = args.runs if args.runs is not None else [args.run]
    tpc, scintillator = _load_pid_scan_tables(args.output_dir, runs)
    scan = scan_charged_pid_thresholds(
        tpc,
        scintillator,
        proton_dedx_values=args.proton_dedx,
        short_range_values=args.short_range,
        short_range_proton_dedx_values=args.short_range_dedx,
    )
    if args.table:
        args.table.parent.mkdir(parents=True, exist_ok=True)
        scan.to_csv(args.table, index=False)

    top = json.loads(scan.head(args.top).to_json(orient="records")) if not scan.empty else []
    summary = {
        "run": runs[0] if len(runs) == 1 else None,
        "runs": runs,
        "scanned_configs": int(len(scan)),
        "labeled_tracks": int(scan.iloc[0]["n_labeled_tracks"]) if not scan.empty else 0,
        "calibration_usable": bool(top[0]["has_both_classes"]) if top else False,
        "best": top[0] if top else None,
        "top": top,
    }
    payload = json.dumps(summary, indent=2, sort_keys=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Reconstruct NNBAR detector Parquet outputs")
    sub = parser.add_subparsers(dest="command", required=True)

    p_summary = sub.add_parser("summarize", help="Reconstruct one run and print summary metrics")
    p_summary.add_argument("output_dir", type=Path, help="Directory containing *_output_<run>.parquet files")
    p_summary.add_argument("--run", type=int, default=0, help="Run number to reconstruct")
    p_summary.add_argument("--json", type=Path, help="Optional JSON summary output path")
    p_summary.add_argument("--tables-dir", type=Path, help="Optional directory for reconstructed CSV tables")
    p_summary.set_defaults(func=summarize)

    p_scan_pid = sub.add_parser("scan-pid", help="Scan charged pion/proton PID thresholds against truth labels")
    p_scan_pid.add_argument("output_dir", type=Path, help="Directory containing *_output_<run>.parquet files")
    p_scan_pid.add_argument("--run", type=int, default=0, help="Run number to scan")
    p_scan_pid.add_argument("--runs", type=_int_grid, help="Comma-separated run numbers to combine")
    p_scan_pid.add_argument(
        "--proton-dedx",
        type=_float_grid,
        default=_float_grid("4,5,6,7,8,9,10"),
        help="Comma-separated proton dE/dx thresholds",
    )
    p_scan_pid.add_argument(
        "--short-range",
        type=_float_grid,
        default=_float_grid("5,10,15,20,25,30"),
        help="Comma-separated short scintillator range thresholds in cm",
    )
    p_scan_pid.add_argument(
        "--short-range-dedx",
        type=_float_grid,
        default=_float_grid("3,4,5,6,7"),
        help="Comma-separated short-range proton dE/dx thresholds",
    )
    p_scan_pid.add_argument("--top", type=int, default=5, help="Number of top configurations to print")
    p_scan_pid.add_argument("--json", type=Path, help="Optional JSON summary output path")
    p_scan_pid.add_argument("--table", type=Path, help="Optional CSV table for all scanned configurations")
    p_scan_pid.set_defaults(func=scan_pid)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
