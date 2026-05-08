"""Parquet I/O utilities for NNBAR simulation outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd
import pyarrow.parquet as pq


OUTPUT_KINDS = (
    "Particle",
    "Interaction",
    "Carbon",
    "Silicon",
    "Beampipe",
    "TPC",
    "Scintillator",
    "LeadGlass",
    "PMT",
    "GPUEnergy",
)


def output_path(output_dir: str | Path, kind: str, run: int = 0) -> Path:
    """Return the expected Parquet path for an output kind and run number."""

    if kind not in OUTPUT_KINDS:
        raise ValueError(f"Unknown output kind {kind!r}; expected one of {OUTPUT_KINDS}")
    return Path(output_dir) / f"{kind}_output_{run}.parquet"


def read_output_table(output_dir: str | Path, kind: str, run: int = 0) -> pd.DataFrame:
    """Read one output table.

    Missing outputs are returned as empty data frames so reconstruction can run
    on partial simulations. macOS AppleDouble files such as ``._TPC...parquet``
    are intentionally ignored by addressing the expected file directly.
    """

    path = output_path(output_dir, kind, run)
    if not path.exists():
        return pd.DataFrame()
    return pq.read_table(path).to_pandas()


def load_run(output_dir: str | Path, run: int = 0) -> dict[str, pd.DataFrame]:
    """Load every known output table for a run."""

    return {kind: read_output_table(output_dir, kind, run) for kind in OUTPUT_KINDS}


def discover_runs(output_dir: str | Path) -> Mapping[str, list[int]]:
    """Discover run numbers by output kind, ignoring non-Parquet sidecar files."""

    runs: dict[str, list[int]] = {kind: [] for kind in OUTPUT_KINDS}
    for path in Path(output_dir).glob("*_output_*.parquet"):
        if path.name.startswith("._"):
            continue
        for kind in OUTPUT_KINDS:
            prefix = f"{kind}_output_"
            if path.name.startswith(prefix):
                suffix = path.name[len(prefix) : -len(".parquet")]
                if suffix.isdigit():
                    runs[kind].append(int(suffix))
                break
    return {kind: sorted(set(values)) for kind, values in runs.items()}
