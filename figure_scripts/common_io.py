"""Read and validate the frozen Paper 1 figure-data contract."""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd


def load_run_metadata(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_gates(meta):
    return meta["gates"]


def _numbered_columns(df, prefix):
    return sorted((c for c in df.columns
                   if re.fullmatch(rf"{prefix}[1-9][0-9]*", c)),
                  key=lambda c: int(c.rsplit("_", 1)[1]))


def _require_columns(df, columns, path):
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if df.empty or not np.isfinite(df[list(columns)].to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: empty or nonfinite numerical data")


def _require_increasing(values, name):
    if len(values) < 2 or not np.isfinite(values).all() or not np.all(np.diff(values) > 0):
        raise ValueError(f"{name}: expected strictly increasing finite coordinates")


def load_mode_profiles(path):
    df = pd.read_csv(path)
    cols = _numbered_columns(df, "u_")
    if not cols:
        raise ValueError(f"{path}: no mode columns")
    _require_columns(df, ["r", *cols], path)
    _require_increasing(df["r"].to_numpy(), path)
    return df, cols


def load_time_series(path):
    df = pd.read_csv(path)
    cols = _numbered_columns(df, "E_mode_")
    if not cols:
        raise ValueError(f"{path}: no mode-energy columns")
    _require_columns(df, ["t", "f_cap", *cols], path)
    _require_increasing(df["t"].to_numpy(), path)
    return df


def resolved_mode_indices(df):
    return [int(c.rsplit("_", 1)[1]) for c in _numbered_columns(df, "E_mode_")]


def load_field_evolution(path):
    with np.load(path, allow_pickle=False) as data:
        r, t, chi = (data[k] for k in ("r", "t_snap", "chi"))
    _require_increasing(r, f"{path}: r")
    _require_increasing(t, f"{path}: t_snap")
    if chi.shape != (len(t), len(r)) or not np.isfinite(chi).all():
        raise ValueError(f"{path}: invalid field array shape or values")
    return r, t, chi
