"""The joint fits of the toy model, made with the package (the owner's message 114: "the fit is always and only a joint fit
across the temperature sweep and/or the power sweep and/or all 8 traces together").

One function fits both sets of traces, the synthetic ones and his recorded ones, with either model of a subfolder: the
model of the detail, or the model that lacks its new term. The parameters that the detail frees are shared by the traces
of the fit. Every other parameter is held at its value of config.toml. Each trace has its own height, centre, offset and
slope of the baseline. The weights are those of the noise law of each condition at the fitted line, settled by the
package in three passes.
"""
from __future__ import annotations

import numpy as np

import twinfit
from twinfit import noise as N

from . import model, traces

OWN = ("height", "centre", "offset", "slope")
TIED = 0.999          # two parameters whose correlation is beyond this are not told apart by the traces of a fit


def tied(cov) -> bool:
    """Whether a covariance says that two parameters are one: a correlation beyond TIED. Where two parameters do the
    same thing to every trace of a fit, the search still ends somewhere, and by the chance of its last steps it may
    give each a finite uncertainty: the correlation says what the uncertainties hide."""
    cov = np.asarray(cov, dtype=float)
    bars = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    if cov.shape[0] < 2 or not np.all(np.isfinite(cov)) or not np.all(bars > 0):
        return False
    corr = cov / np.outer(bars, bars)
    return bool(np.max(np.abs(corr[np.triu_indices(cov.shape[0], 1)])) > TIED)


def start_of(record) -> list:
    """Where the search of one trace starts: the top above the ends, the place of the top, the level of the ends, no slope."""
    y, nu = np.asarray(record.y, dtype=float), np.asarray(record.x, dtype=float)
    ends = np.concatenate([y[:40], y[-40:]])
    smooth = np.convolve(y, np.ones(31) / 31.0, mode="same")
    return [float(smooth.max() - np.median(ends)), float(nu[int(np.argmax(smooth))]), float(np.median(ends)), 0.0]


def curve(detail: int, const: dict, config: dict, free: list):
    """The model of one trace as the package asks it: the shared parameters are those that are free, in their order."""
    held = model.truth(config)

    def of(x, shared, own, design):
        values = dict(held, **{name: float(v) for name, v in zip(free, shared)})
        u = np.asarray(x, dtype=float) - own[1]
        top = float(model.shape(detail, 0.0, design["T"], values, const))
        return own[0] * model.shape(detail, u, design["T"], values, const) / top + own[2] + own[3] * u
    return of


def variance_of(const: dict):
    """The law of the noise of a condition, at the fitted line: what the package sets its weights from."""
    def of(mean, rec, own):
        u = rec.x - own[1]
        signal = np.clip(np.asarray(mean, dtype=float) - own[2] - own[3] * u, 0.0, None)
        return traces.noise_sd(rec.design["P"], rec.design["T"], signal, rec.x, const) ** 2
    return of


def one_value_for_all(rows: list, const: dict, config: dict, free: list) -> None:
    """A value of theory is one value for all 8 traces (his message 135: "make sure that all 8 traces will have always
    the same w_0 and alpha and beta_self"). It is held in every fit, and a fit that frees it as a test holds all 8
    traces: freed on a sweep alone it gave the trace that the two sweeps share two values. The beam radius is no
    parameter of a fit at all: it is in const.toml, and the transit rate that it gives is of kind theory."""
    theory = [name for name in free if config["parameter"][name].get("kind") == "theory"]
    held = sorted((row["P"], row["T"]) for row in rows)
    if theory and held != sorted(traces.conditions(const)):
        raise ValueError(f"{', '.join(theory)} of kind theory would be freed on {len(held)} traces: a value of theory is one for all 8, and is freed only in the joint fit of all 8")


def problem_of(detail: int, rows: list, const: dict, config: dict, free: list | None = None) -> twinfit.Problem:
    """The joint problem of some traces under the model of a detail, with the parameters that the detail frees, or with
    those of `free` alone: the others are held at their values of config.toml."""
    free = list(config["detail"][str(detail)]["free"] if free is None else free)          # an empty list is a model with nothing shared to fit: what theory gives is held, and each trace keeps its own numbers
    one_value_for_all(rows, const, config, free)
    entries = [config["parameter"][name] for name in free]
    records = []
    for row in rows:
        y = np.asarray(row["y"], dtype=float)
        first = float(np.std(np.concatenate([y[:40], y[-40:]]), ddof=1))
        records.append(twinfit.Record(row["nu"], y, max(first, 1e-12), {"P": row["P"], "T": row["T"]}, name=traces.name_of(row["P"], row["T"])))
    return twinfit.Problem(records, curve(detail, const, config, free), tuple(free), OWN, [e["start"] for e in entries], start_of,
                           centre=lambda s, o, d: o[1], name=f"detail {detail}",
                           bounds_shared=([e["lower"] for e in entries], [e["upper"] for e in entries]))


def fit(detail: int, rows: list, const: dict, config: dict, free: list | None = None) -> dict:
    """One joint fit of some traces under the model of a detail. The values and the uncertainties of the free parameters,
    the numbers of each trace, the chi square, and for each trace its fitted line and its residuals over the noise. A
    fit whose parameters these traces do not tell apart says so, and gives no uncertainty for them."""
    problem = problem_of(detail, rows, const, config, free)
    settled, result, moved = N.settle(problem, variance_of(const))
    free = list(settled.shared)
    bars = np.sqrt(np.clip(np.diag(result.cov), 0.0, None))
    told = bool(result.converged and np.all(np.isfinite(bars)) and np.all(bars > 0) and not result.at_bound and not tied(result.cov))
    out = {"detail": detail, "free": free, "conditions": [(row["P"], row["T"]) for row in rows], "converged": bool(result.converged), "told": told,
           "values": {name: float(v) for name, v in zip(free, result.shared)}, "bars": {name: (float(b) if told else None) for name, b in zip(free, bars)},
           "chi2": float(result.chi2), "dof": int(result.dof), "weights_moved_in_bars": None if moved is None else float(moved), "traces": []}
    for k, rec in enumerate(settled.records):
        line = settled.predict(k, result.shared, result.own[k])
        out["traces"].append({"P": rec.design["P"], "T": rec.design["T"], "own": dict(zip(OWN, (float(v) for v in result.own[k]))),
                              "nu": rec.x, "y": rec.y, "line": line, "residual": (rec.y - line) / rec.sd, "sd": rec.sd})
    return out


def three(detail: int, rows: list, const: dict, config: dict, free: list | None = None) -> dict:
    """The three joint fits of a subfolder under one model: the power sweep, the temperature sweep, all 8 traces."""
    by = {(row["P"], row["T"]): row for row in rows}
    return {name: fit(detail, [by[c] for c in which], const, config, free) for name, which in traces.sweeps(const).items()}


def tested(detail: int, rows: list, const: dict, config: dict, sweeps=None) -> dict | None:
    """What the fits of a detail hold from theory, tested on some traces (his choice of 2026-10-14: "Hold theory, free it
    as a test"). None for a model that holds nothing of theory. For each joint fit of `three`, or for those of `sweeps`:

    `freed`    the fit of all 8 traces with the values of theory free as well, and on a sweep alone None: a value of
               theory is one for all 8 traces (his message 135, `one_value_for_all`). For each of them its fitted value
               and uncertainty, and how far it lies from theory, in the two uncertainties together. A fit that cannot
               tell its parameters apart gives no distance.
    `theory`   the uncertainty that theory puts on each value of the main fit. The main fit is made again with one held
               value an uncertainty below theory and an uncertainty above it: half of what a fitted value moves by is
               the share of that held value, in `moves` with its sign, and the shares add in squares. These are fits
               made again and no first order (CLAUDE.md, rule 5b)."""
    import copy
    entry = config["detail"][str(detail)]
    free = list(entry["free"])
    names = [name for name in entry.get("freed", free) if name not in free]
    if not names:
        return None
    by = {(row["P"], row["T"]): row for row in rows}
    of_theory = {name: {"value": float(config["parameter"][name]["true"]), "uncertainty": float(config["parameter"][name]["uncertainty"])} for name in names}
    out = {"detail": detail, "held": of_theory, "free": free, "sweeps": {}}
    for sweep, which in traces.sweeps(const).items():
        if sweeps is not None and sweep not in sweeps:
            continue
        some = [by[c] for c in which]
        other = fit(detail, some, const, config, list(entry["freed"])) if sweep == "all" else None
        apart = {name: (other["values"][name] - t["value"]) / float(np.hypot(other["bars"][name], t["uncertainty"])) if other["told"] else None for name, t in of_theory.items()} if other else None
        moves, told = {}, True
        for name, t in of_theory.items():
            ends = []
            for side in (-1.0, 1.0):
                moved = copy.deepcopy(config)
                moved["parameter"][name]["true"] = t["value"] + side * t["uncertainty"]
                ends.append(fit(detail, some, const, moved, free))
            told = told and all(end["told"] for end in ends)
            moves[name] = {k: (ends[1]["values"][k] - ends[0]["values"][k]) / 2.0 for k in free}
        out["sweeps"][sweep] = {"freed": other and {"told": other["told"], "values": other["values"], "bars": other["bars"], "chi2": other["chi2"], "dof": other["dof"], "from_theory": apart},
                                "moves": moves, "theory": {k: (float(np.sqrt(sum(moves[name][k] ** 2 for name in names))) if told else None) for k in free}}
    return out


def lacking(detail: int, config: dict) -> list:
    """The free parameters of the model that lacks the new term of a detail: those of the detail, less the ones that the
    detail brings. Everything else is as in the model with the term, so that the two differ by the term alone. With the
    natural width free in one and held in the other, the model without the gas was the model with it under another name
    (review 1, the severe finding)."""
    return [name for name in config["detail"][str(detail)]["free"] if config["parameter"][name]["enters"] < detail]


def synthetic(detail: int, const: dict, config: dict, draw: int = 0) -> list:
    """The 8 synthetic traces of a detail as rows for a fit."""
    return [traces.trace(detail, P, T, const, config, draw) for P, T in traces.conditions(const)]


RIGHT = 0.6827          # the share of the draws on which the interval of one right uncertainty holds the truth


def held(draws: int) -> tuple:
    """The band in which the share of intervals that hold the truth must lie, on a number of draws: three standard errors
    about the share of a right uncertainty, and never wider than 0.50 to 0.85."""
    reach = 3.0 * float(np.sqrt(RIGHT * (1.0 - RIGHT) / max(int(draws), 1)))
    return max(0.50, RIGHT - reach), min(0.85, RIGHT + reach)


def check(detail: int, const: dict, config: dict, draws: int = 100, got: list | None = None, free: list | None = None) -> dict:
    """The own fit of a detail on many draws of its synthetic traces, all 8 together: for each free parameter the bias over
    the spread of its estimates, the mean uncertainty over that spread, and the share of the draws whose interval of one
    uncertainty holds the true value. `passes` says whether the model may go on his recorded traces: every search ended,
    a bias under a quarter of the spread, and intervals that hold the truth on a share of the draws that lies within
    three standard errors of 0.6827, the share of a right uncertainty, and never outside 0.50 to 0.85. On 100 draws that
    is 0.543 to 0.822. The band was 0.55 to 0.80, set by hand, until the model of detail 4 gave 0.81 on its first 100
    draws and 0.68 to 0.69 on 300 more: a right model was refused by a band that was narrower above than below. `got`
    holds the fits when they were made one by one and kept, as `draws.py` does: fewer of them than `draws` does not pass."""
    true = model.truth(config)
    free = list(config["detail"][str(detail)]["free"] if free is None else free)
    got = [fit(detail, synthetic(detail, const, config, 1 + k), const, config, free) for k in range(draws)] if got is None else list(got)
    ended = [g for g in got if g["told"]]
    out = {"detail": detail, "draws": draws, "ended": len(ended), "parameters": {}}
    for name in free:
        values = np.array([g["values"][name] for g in ended])
        bars = np.array([g["bars"][name] for g in ended])
        spread = float(values.std(ddof=1)) if len(ended) > 2 else float("nan")
        out["parameters"][name] = {"true": true[name], "mean": float(values.mean()) if len(ended) else float("nan"), "spread": spread,
                                   "bias_over_spread": float((values.mean() - true[name]) / spread) if spread > 0 else float("nan"),
                                   "bar_over_spread": float(bars.mean() / spread) if spread > 0 else float("nan"),
                                   "holds_the_truth": float(np.mean(np.abs(values - true[name]) <= bars)) if len(ended) else float("nan")}
    low, high = held(draws)
    out["band_of_the_share_that_holds_the_truth"] = [low, high]
    out["passes"] = bool(len(got) == draws and len(ended) == draws and all(abs(p["bias_over_spread"]) < 0.25 and low <= p["holds_the_truth"] <= high for p in out["parameters"].values()))
    return out
