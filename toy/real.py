"""His recorded traces, read like the synthetic ones (the owner's message 117: "the point is to compare the real traces with
the synthetic ones. So treat them with the same dignity").

They are files of this repository, `data/recorded/<peak>/`, in the form of the synthetic traces: the frequency in MHz
from the top of the trace and the signal in V. `data/recorded/INDEX.json` says where each came from: the file of the
record, its MD5, which repetition it is, and whether its time stamps were put in equal steps. The peak is the one that
config.toml names, 993.4192 nm by default.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import traces

ROOT = Path(__file__).resolve().parents[1]


def index(root: Path = ROOT) -> list:
    """The list of the recorded traces of the repository, one entry a file."""
    return json.loads((Path(root) / "data" / "recorded" / "INDEX.json").read_text())["traces"]


def rows(const: dict, config: dict, root: Path = ROOT) -> list:
    """His 8 recorded traces of the peak of config.toml, the first repetition of each condition, as rows for a fit."""
    peak = config["recorded"]["peak"]
    if config["recorded"].get("repetition", "first") != "first":
        raise ValueError(f"config.toml asks for the repetition {config['recorded']['repetition']!r}: only the first one of each condition is read")
    listed = {}
    for e in sorted(index(root), key=lambda e: -e["repetition"]):          # the smallest number of repetition is kept, whatever the index holds
        if e["peak_nm"] == peak:
            listed[(e["power_mW"], e["temperature_C"])] = e
    if not listed:
        raise ValueError(f"the repository holds no recorded trace of the peak {peak} nm: it holds {', '.join(sorted({e['peak_nm'] for e in index(root)}))}")
    out = []
    for P, T in traces.conditions(const):
        if (P, T) not in listed:
            raise ValueError(f"no recorded trace of the peak {peak} nm at {P} mW and {T} C")
        entry = listed[(P, T)]
        table = np.loadtxt(Path(root) / entry["file"], delimiter=",", skiprows=1)
        out.append({"P": P, "T": T, "nu": table[:, 0], "y": table[:, 1], "file": entry["file"], "repeat": entry["repetition"]})
    return out
