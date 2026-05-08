"""Command-line interface for NNBAR offline reconstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Reconstruct NNBAR detector Parquet outputs")
    sub = parser.add_subparsers(dest="command", required=True)

    p_summary = sub.add_parser("summarize", help="Reconstruct one run and print summary metrics")
    p_summary.add_argument("output_dir", type=Path, help="Directory containing *_output_<run>.parquet files")
    p_summary.add_argument("--run", type=int, default=0, help="Run number to reconstruct")
    p_summary.add_argument("--json", type=Path, help="Optional JSON summary output path")
    p_summary.add_argument("--tables-dir", type=Path, help="Optional directory for reconstructed CSV tables")
    p_summary.set_defaults(func=summarize)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
