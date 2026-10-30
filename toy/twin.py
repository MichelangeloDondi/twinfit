"""The digital twin of his 8 recorded traces: the model of the toy drawn with what his traces measure of each trace, his
messages 305 and 307, "yes, the twin takes my heights" and "yes, the twin carries a laser noise for each trace".

The model is the toy's, `toy/model.py`, with one laser noise for all 8 traces, and its check on 100 draws and the counts
of the toy are made on its own traces, `fits.synthetic`. The twin draws each trace as `traces.trace` does, from the same
random numbers, with the laser noise that his trace takes alone and with his height at its top, and the noise drawn
again at that height by the law of its condition. Both are read from `toy/twin.toml`, which `python toy/adopt.py heights
--write` and `python toy/adopt.py laser-noise FILE --write` write. The file and this program stand apart from the model
and its two files, so that the counts of the toy, which read those, stay as they are.
"""
from __future__ import annotations

import copy
import json
import tomllib
from pathlib import Path

import numpy as np

from . import model, traces

HERE = Path(__file__).resolve().parent
FILE = "twin.toml"


def load(folder: Path = HERE) -> dict:
    """What his traces measure of each trace, as toy/twin.toml holds it: nothing where the file is not."""
    path = Path(folder) / FILE
    return tomllib.loads(path.read_text()) if path.exists() else {}


def cell(table: dict, key: str, P_mW: float, T_celsius: float) -> float | None:
    """The number of a table of the twin at one condition, or None where the table does not hold it."""
    for P, T, v in zip(table.get("power_mW", []), table.get("temperature_C", []), table.get(key, [])):
        if (float(P), float(T)) == (float(P_mW), float(T_celsius)):
            return float(v)
    return None


def height_of(P_mW: float, T_celsius: float, twin: dict) -> float | None:
    """His height at a condition, in V, or None where the twin holds none."""
    return cell(twin.get("heights", {}), "height_V", P_mW, T_celsius)


def laser_noise_of(P_mW: float, T_celsius: float, value: float, twin: dict) -> float:
    """The laser noise of one trace of the twin: the value given times the ratio of the trace's own to the joint value
    of all 8 that the table was written beside, so that at the twin's value each trace has its own exactly, and a laser
    noise planted at another size keeps the pattern of his traces. The value given where the table holds no trace."""
    table = twin.get("laser_noise_of_each_trace", {})
    own = cell(table, "laser_noise_MHz", P_mW, T_celsius)
    return float(value) if own is None else float(value) * own / float(table["of_the_joint_MHz"])


def trace(detail: int, P_mW: float, T_celsius: float, const: dict, config: dict, draw: int = 0, twin: dict | None = None) -> dict:
    """One trace of the twin: the trace that `traces.trace` draws with the laser noise of that trace, scaled to his
    height, with its noise drawn again at that height from the same random numbers. Where the twin holds nothing, the
    model's own trace."""
    twin = load() if twin is None else twin
    made = config
    if "laser_noise" in config["parameter"]:
        made = copy.deepcopy(config)
        made["parameter"]["laser_noise"]["true"] = laser_noise_of(P_mW, T_celsius, model.truth(config)["laser_noise"], twin)
    row = traces.trace(detail, P_mW, T_celsius, const, made, draw)
    height = height_of(P_mW, T_celsius, twin)
    if height is None:
        return row
    z = (np.asarray(row["y"]) - np.asarray(row["signal"])) / np.asarray(row["sd"])
    signal = np.asarray(row["signal"]) * height / float(np.max(row["signal"]))
    sd = traces.noise_sd(P_mW, T_celsius, signal, row["nu"], const)
    return {**row, "signal": signal, "sd": sd, "y": signal + sd * z}


def synthetic(detail: int, const: dict, config: dict, draw: int = 0, twin: dict | None = None) -> list:
    """The twin's 8 traces of a detail as rows for a fit, as `fits.synthetic` gives the model's own."""
    twin = load() if twin is None else twin
    return [trace(detail, P, T, const, config, draw, twin) for P, T in traces.conditions(const)]


def write(detail: int, folder: Path, const: dict, config: dict, twin: dict | None = None) -> list:
    """The files of the twin's traces of a detail in `folder`, as `traces.write` makes the model's: `traces/` with one CSV
    file a trace, and `truth/truth.json` with the values that made them, the laser noise and the height of each trace
    among them, which an analyst is not given."""
    twin = load() if twin is None else twin
    folder = Path(folder)
    (folder / "traces").mkdir(parents=True, exist_ok=True)
    (folder / "truth").mkdir(parents=True, exist_ok=True)
    written, rows, truth = [], [], model.truth(config)
    for P, T in traces.conditions(const):
        t = trace(detail, P, T, const, config, 0, twin)
        path = folder / "traces" / traces.name_of(P, T)
        path.write_text("frequency_MHz,signal_V\n" + "".join(f"{x:.6f},{y:.9e}\n" for x, y in zip(t["nu"], t["y"])))
        written.append(path)
        rows.append({"file": path.name, "power_mW": P, "temperature_C": T, "height_V": float(t["signal"].max()),
                     "lorentzian_width_MHz": float(model.width(detail, T, truth, const)), "density_per_cm3": float(model.density(T, const))}
                    | ({"laser_noise_MHz": laser_noise_of(P, T, truth["laser_noise"], twin)} if config["parameter"]["laser_noise"]["enters"] <= detail else {}))
    held = {name: value for name, value in truth.items() if config["parameter"][name]["enters"] <= detail}
    (folder / "truth" / "truth.json").write_text(json.dumps({"detail": detail, "parameters": held, "traces": rows}, indent=1) + "\n")
    return written + [folder / "truth" / "truth.json"]
