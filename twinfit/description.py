"""The description of a problem: what a researcher writes before any fit, for the package to read (the owner's message 140).

The input of the package, with or without an agent, is the description of the data collection and nothing more: how
many traces, in which conditions, what the physical problem is. Optional addenda may follow: an initial model, papers to
read, manuals of the instruments, ansatz about parameters with their uncertainty. The description is one file in TOML,
and its paths are read from the folder that holds it:

    [problem]
    name = "..."
    what = "the physical problem, in words, for the researcher or the agent who chooses the terms"
    abscissa = "time"               the name of the abscissa, and its unit
    abscissa_unit = "s"
    abscissa_offset = 0.0           optional: a number taken from the abscissa of every trace as it is read, so that the
                                    curves sit near zero, as 587.6 for lines near 587.6 nm. The model then reads the
                                    abscissa less it, and so does every place that the output gives
    value = "signal"
    value_unit = "V"

    [covariates.temperature]        each condition under which a trace was taken, with its unit
    unit = "C"

    [[trace]]                       one table a trace: its file, its columns, the value of each covariate
    file = "traces/one.csv"
    columns = ["time_s", "signal_V"]
    temperature = 70

    [noise]                         the standard deviation of each sample: one column of the files, one number, or read
    from = "neighbours"             from the differences of neighbouring samples, as twinfit.noise.flat does

    [setup.gain]                    optional: each instrument and setting of the data collection, in words, with its
    words = "the gain of the        value, its uncertainty and its unit where they are known, and whether it was
      detector"                     measured, assumed or not recorded. The package reads no physics in it: it is for
    value = 1e6                     whoever chooses the terms, and a lever of kind setup names one of them
    unit = "V/A"
    known = "set"                   one of measured, set, from its maker, assumed, not recorded

    [addenda]                       all optional
    model = "model.py"              a file that defines model, SHARED, OWN and start_own, as agent/problem_template.py
    papers = ["paper.pdf"]
    manuals = ["manual.pdf"]
    [addenda.ansatz.width]          a value that the researcher holds before the fit, with its uncertainty: a covariate
    value = 3.5                     that every record gives at that value, or a key of the dictionary HELD of the model
    uncertainty = 0.75              file, which the model reads. The loop refuses an ansatz on anything else, a fitted
                                    parameter among them, since it would put no bar on the others (review 22, finding 6)
    [addenda.lever.gain]            a change of the experiment that the researcher could make, and how it enters the twin:
    kind = "noise"                  one of the kinds of proposals.KINDS, with its factor, or for a covariate its name
    factor = 0.5                    and its new values, and in `at` the conditions of the traces that the new records
                                    start from. The words say what the change is, for the researcher
    words = "a higher gain of the detector, which halves the noise of each sample"
    [addenda.lever.width]           a better characterization of the setup: the uncertainty of a value held in the
    kind = "characterization"       ansatz, `of`, times the factor, under one
    of = "width"
    factor = 0.5
    words = "a measurement of the width before the run, to half its uncertainty"
    [addenda.lever.mirror]          a change of an item of the setup, `of`: measured on the twin where the model holds a
    kind = "setup"                  value of that name, as a characterization with its factor, and otherwise named in
    of = "gain"                     the proposals with the words that no term of the model reads it yet
    words = "the gain measured at each power"

    description = read("description.toml")          refuses what is missing or does not agree, in plain words
    records = description.records()                 one twinfit.Record a trace, its design the covariates
    description.varies(), description.fixed()       the covariates that change across the traces, and those that do not
    description.levers()                            the levers of the experiment, each a dict with its name
    description.setup                               the items of the setup, each a dict with its words

A covariate that does not change across the traces is one along which no reading can tell anything: the report of where
to search says so of it, and it is said here once. The package reads no physics in the description: its words are for
whoever chooses the terms.
"""
from __future__ import annotations

import csv
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .noise import flat
from .problem import Record

NOISE = ("column", "value", "neighbours")          # the three ways a description may give the noise of each sample
LEVERS = ("noise", "samples", "span", "repeats", "covariate", "characterization", "setup")
KNOWN = ("measured", "set", "from its maker", "assumed", "not recorded")          # how an item of the setup is known: a gain is set, an aperture is from its maker          # the ways a lever enters the twin, as proposals.KINDS says each


@dataclass
class Description:
    """A problem as its researcher describes it, read from one file. `home` is the folder its paths are read from."""
    name: str
    what: str
    abscissa: tuple
    value: tuple
    covariates: dict
    traces: list
    noise: dict
    addenda: dict = field(default_factory=dict)
    home: Path = Path(".")
    offset: float = 0.0                  # taken from the abscissa of every trace as it is read (round 5 of the cold tests)
    order: str = ""                      # in which order the curves were taken, which a drift follows (round 10 of the cold tests)
    setup: dict = field(default_factory=dict)          # the instruments and settings of the data collection, his message 184: the agent chose terms knowing the conditions and the noise, and nothing of the setup

    def varies(self) -> list:
        """The covariates whose value changes across the traces: the only ones along which a behaviour can be read."""
        return [name for name in self.covariates if len({t[name] for t in self.traces}) > 1]

    def fixed(self) -> list:
        """The covariates that take one value in every trace: nothing can be read along them."""
        return [name for name in self.covariates if name not in self.varies()]

    def levers(self) -> list:
        """The levers of the experiment that the addenda give, each a dict with its name, its kind and its words."""
        return [dict(entry, name=name) for name, entry in (self.addenda.get("lever") or {}).items()]

    def records(self) -> list:
        """One Record a trace: its abscissa rising, its values, the standard deviation of each sample from the noise of
        the description, and its design, the value of each covariate."""
        out = []
        for k, trace in enumerate(self.traces):
            with open(self.home / trace["file"], newline="") as f:
                rows = list(csv.DictReader(f))
            if not rows:
                raise ValueError(f"trace {k + 1}: {trace['file']} holds no line of values under its first line")
            columns = trace.get("columns") or list(rows[0])[:2]
            lacking = [c for c in columns + ([self.noise["column"]] if self.noise["from"] == "column" else []) if c not in rows[0]]
            if lacking:          # the first line names the columns, separated by commas: a comment, a tab or another name left a KeyError with no file in it (the cold tests of row M169)
                raise ValueError(f"trace {k + 1}: {trace['file']} has the columns {', '.join(c for c in rows[0] if c) or 'none'} on its first line, and the description asks {', '.join(lacking)}: "
                                 "the first line of a file names its columns, separated by commas")
            x = np.array([float(r[columns[0]]) for r in rows]) - self.offset
            order = np.argsort(x, kind="stable")          # sorted before the noise of neighbours is read: rows in another order read it five times too large (round 7 of the cold tests)
            x, y = x[order], np.array([float(r[columns[1]]) for r in rows])[order]
            how = self.noise["from"]
            if how == "column":
                sd = np.array([float(r[self.noise["column"]]) for r in rows])[order]
            elif how == "value":
                sd = np.full(x.shape, float(self.noise["value"]))
            else:
                sd = np.full(x.shape, flat(x, y))
            out.append(Record(x, y, sd, {name: trace[name] for name in self.covariates}, name=trace.get("name", trace["file"])))
        return out


KEYS = {"": {"problem", "covariates", "trace", "noise", "addenda", "setup"}, "setup.*": {"words", "value", "uncertainty", "unit", "known"},
        "problem": {"name", "what", "abscissa", "abscissa_unit", "value", "value_unit", "abscissa_offset", "order"},
        "covariates.*": {"unit"}, "noise": {"from", "value", "column"},
        "addenda": {"model", "ansatz", "lever", "papers", "manuals", "notes"},
        "addenda.ansatz.*": {"value", "uncertainty"}, "addenda.lever.*": {"kind", "factor", "words", "covariate", "values", "of", "at"}}          # what a description may hold: a lever under [addenda.lvr] was dropped in silence (round 10 of the cold tests)


def unknown(data: dict) -> list:
    """The keys and tables of a description that the package does not read, each with its place: a misspelt one is
    otherwise dropped in silence, as a lever under [addenda.lvr], and a key that the package does not know, as an order
    of the curves under a name of one's own, is read by nothing."""
    out = [k for k in data if k not in KEYS[""]]
    for table in ("problem", "noise", "addenda"):
        out += [f"{table}.{k}" for k in (data.get(table) or {}) if k not in KEYS[table]]
    for name, entry in (data.get("covariates") or {}).items():
        out += [f"covariates.{name}.{k}" for k in (entry or {}) if k not in KEYS["covariates.*"]]
    for name, entry in (data.get("setup") or {}).items():
        out += [f"setup.{name}.{k}" for k in (entry or {}) if k not in KEYS["setup.*"]]
    for kind in ("ansatz", "lever"):
        for name, entry in ((data.get("addenda") or {}).get(kind) or {}).items():
            out += [f"addenda.{kind}.{name}.{k}" for k in (entry or {}) if k not in KEYS[f"addenda.{kind}.*"]]
    return out


def read(path) -> Description:
    """The description in a file, checked: every trace has its file and a value of every covariate, every covariate has
    its unit, the noise is given one of the three ways, and every ansatz has its uncertainty. A refusal says what is
    missing, so that the researcher can write it."""
    path = Path(path)
    data = tomllib.loads(path.read_text())
    strange = unknown(data)
    if strange:
        raise ValueError(f"the description {path.name} holds {', '.join(strange)}, which the package reads nowhere: a misspelt name, or a key of its own. "
                         "The order in which the curves were taken goes in [problem] as order = \"...\", and a note of your own in [addenda] as notes = \"...\"")
    problem = data.get("problem") or {}
    missing = [key for key in ("name", "what", "abscissa", "value") if not str(problem.get(key, "")).strip()]
    if missing:
        raise ValueError(f"the description {path.name} gives no {', '.join(missing)} in [problem]: the package and whoever chooses the terms read them")
    covariates = data.get("covariates") or {}
    for name, entry in covariates.items():
        if not str((entry or {}).get("unit", "")).strip():
            raise ValueError(f"the covariate {name!r} has no unit: write unit = \"...\", or \"none\"")
    traces = data.get("trace") or []
    if not traces:
        raise ValueError(f"the description {path.name} has no [[trace]]: it describes no data")
    for k, trace in enumerate(traces):
        if not (path.parent / str(trace.get("file", ""))).is_file():
            raise ValueError(f"trace {k + 1}: its file {trace.get('file')!r} is not found next to the description")
        lacking = [name for name in covariates if name not in trace]
        if lacking:
            raise ValueError(f"trace {k + 1}, {trace['file']}: no value of {', '.join(lacking)}. Every trace gives every covariate")
    noise = data.get("noise") or {"from": "neighbours"}
    if noise.get("from") not in NOISE:
        raise ValueError(f"the noise is given from {noise.get('from')!r}: it is one of {', '.join(NOISE)}")
    if noise["from"] == "column" and not noise.get("column"):
        raise ValueError("the noise is given from a column, and no column is named: write column = \"...\"")
    if noise["from"] == "value" and not float(noise.get("value", 0.0)) > 0.0:
        raise ValueError("the noise is given as one value, and none above zero is written: write value = ...")
    setup = data.get("setup") or {}
    for name, item in setup.items():
        if not isinstance(item, dict) or not str(item.get("words", "")).strip():
            raise ValueError(f"the item {name!r} of the setup has no words: write words = \"...\", what the instrument or the setting is")
        if any(k in item and not isinstance(item[k], (int, float)) for k in ("value", "uncertainty")) or ("uncertainty" in item and "value" not in item):
            raise ValueError(f"the item {name!r} of the setup: value and uncertainty are numbers, and an uncertainty goes with a value")
        if item.get("known", "measured") not in KNOWN:
            raise ValueError(f"the item {name!r} of the setup is known as {item.get('known')!r}: write known = one of {', '.join(KNOWN)}")
    addenda = data.get("addenda") or {}
    for name, ansatz in (addenda.get("ansatz") or {}).items():
        if "value" not in ansatz or not float(ansatz.get("uncertainty", 0.0)) > 0.0:
            raise ValueError(f"the ansatz on {name!r} has no value or no uncertainty above zero: a value held before the fit carries the bar that it puts on the others")
    for name, lever in (addenda.get("lever") or {}).items():
        kind = lever.get("kind")
        if kind not in LEVERS:
            raise ValueError(f"the lever {name!r} is of kind {kind!r}: it is one of {', '.join(LEVERS)}")
        if not str(lever.get("words", "")).strip():
            raise ValueError(f"the lever {name!r} has no words: write words = \"...\", what the change of the experiment is")
        if kind == "covariate":
            if lever.get("covariate") not in covariates or not lever.get("values") or not all(isinstance(v, (int, float)) for v in lever["values"]):
                raise ValueError(f"the lever {name!r} adds records along a covariate: write covariate = one of {', '.join(covariates) or 'none, since none is given'}, and values = [...], numbers")
            at = lever.get("at") or {}
            if any(key not in covariates for key in at) or not any(all(t.get(key) == value for key, value in at.items()) for t in traces):
                raise ValueError(f"the lever {name!r} starts from the trace at {at}: no trace is at those conditions")
        elif kind == "characterization":
            if lever.get("of") not in (addenda.get("ansatz") or {}) or not (isinstance(lever.get("factor"), (int, float)) and 0 < lever["factor"] < 1):
                held = sorted(addenda.get("ansatz") or {})          # which ansatz there are, or that there is none: a tester wrote of = a fitted parameter and the refusal did not say why (round 6 of the cold tests)
                raise ValueError(f"the lever {name!r} narrows a value held by the description: write of = the name of an ansatz, " + (f"one of {', '.join(held)}" if held else "and the description holds none: an ansatz under [addenda.ansatz.NAME] first")
                                 + (f", where it names {lever.get('of')!r}" if lever.get("of") is not None and lever.get("of") not in held else "") + ", and factor = the share of its uncertainty that a better characterization leaves, above zero and under one")
        elif kind == "setup":
            if lever.get("of") not in setup:
                raise ValueError(f"the lever {name!r} changes an item of the setup: write of = the name of an item under [setup], " + (f"one of {', '.join(sorted(setup))}" if setup else "and the description gives none"))
            if "factor" in lever and not (isinstance(lever["factor"], (int, float)) and 0.0 < lever["factor"] < 1.0):
                raise ValueError(f"the lever {name!r}: its factor is the share of the uncertainty of the item that the change leaves, above zero and under one")
        elif kind == "repeats":
            if not (isinstance(lever.get("factor"), int) and lever["factor"] >= 2):
                raise ValueError(f"the lever {name!r} repeats each record: write factor = a whole number, two or more")
        elif not (isinstance(lever.get("factor"), (int, float)) and lever["factor"] > 0):
            raise ValueError(f"the lever {name!r} multiplies by a factor: write factor = a number above zero")
    if addenda.get("model") and not (path.parent / addenda["model"]).is_file():
        raise ValueError(f"the initial model {addenda['model']!r} is not found next to the description")
    return Description(name=problem["name"], what=problem["what"], abscissa=(problem["abscissa"], problem.get("abscissa_unit", "")), value=(problem["value"], problem.get("value_unit", "")),
                       covariates=covariates, traces=traces, noise=noise, addenda=addenda, home=path.parent, offset=float(problem.get("abscissa_offset", 0.0)), order=str(problem.get("order", "")),
                       setup={name: dict(item, known=item.get("known", "measured")) for name, item in setup.items()})
