"""Two worked projects of other fields, written into a folder as a researcher writes one: the curves, the description of
their collection, the first model, and no term yet. The curves are drawn here from a truth that the first model lacks a
part of. What that part is, the folder does not say: the report of the package is there to find it.

    python -m twinfit example chromatography DIR
    python -m twinfit example nmr DIR

chromatography   the peak of one compound at four temperatures of the column, each on its own height, place and
                 baseline, of one width. An impurity elutes just after it, small, which the first model lacks
nmr              one line at four temperatures of the sample, each on its own height, place and baseline, of one width.
                 Its width grows with the temperature, which the first model does not hold

Each folder runs whole, `python -m twinfit loop DIR`, in about a minute on one core.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

FIELDS = {
    "chromatography": {"what": "the peak of one compound eluting from a column held at four temperatures. Each chromatogram has its own height, retention time and baseline. The first model is one Gaussian peak of one width for all four",
                       "abscissa": ("time", "s"), "value": ("absorbance", "mAU"), "covariate": ("temperature", "C"), "levels": (30.0, 40.0, 50.0, 60.0),
                       "span": 30.0, "step": 0.1, "noise": 0.02, "start": 0.0, "at": 12.0},          # a time from the injection, so that no retention time is negative (the cold tests of row M169)
    "nmr": {"what": "one resonance of a sample measured at four temperatures. Each spectrum has its own height, frequency and baseline. The first model is one Lorentzian line of one width for all four",
            "abscissa": ("frequency", "Hz"), "value": ("intensity", "a.u."), "covariate": ("temperature", "K"), "levels": (290.0, 300.0, 310.0, 320.0),
            "span": 40.0, "step": 0.1, "noise": 0.01, "start": -20.0, "at": 0.0},          # a frequency from that of the reference
}
MODELS = {
    "chromatography": '''"""The first model of the chromatograms: one Gaussian peak of one width, on its own height, place and baseline."""
import numpy as np

SHARED = ("width",)                     # the standard deviation of the peak, in s, the same for every chromatogram
OWN = ("height", "place", "baseline")   # of each chromatogram
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.05], [10.0])
CENTRE = "place"
NAME = "one Gaussian peak of one width"


def model(x, shared, own, design):
    return own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    return [float(record.y.max() - record.y.min()), float(record.x[np.argmax(record.y)]), float(np.median(record.y))]


def centre(shared, own, design):
    return own[1]
''',
    "nmr": '''"""The first model of the spectra: one Lorentzian line of one width, on its own height, place and baseline."""
import numpy as np

SHARED = ("width",)                     # the full width of the line at half its height, in Hz, the same for every spectrum
OWN = ("height", "place", "baseline")   # of each spectrum
START_SHARED = [2.0]
BOUNDS_SHARED = ([0.05], [20.0])
CENTRE = "place"
NAME = "one Lorentzian line of one width"


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    return [float(record.y.max() - record.y.min()), float(record.x[np.argmax(record.y)]), float(np.median(record.y))]


def centre(shared, own, design):
    return own[1]
''',
}


TRUTH = {"chromatography": {"width": 1.0, "impurity": 0.06, "after": 1.4},          # the impurity: its height over that of the peak, and how many widths later it elutes
         "nmr": {"width": 2.0, "at": 290.0, "growth": 0.004}}                       # the width at 290 K, and its growth, per kelvin over that width


def truth(field: str, x: np.ndarray, level: float, k: int) -> np.ndarray:
    """The curve of one record, with the part that the first model lacks, from the values of TRUTH."""
    height, place, base = (1.0 + 0.1 * k, FIELDS[field]["at"] + 0.3 * (k - 1.5), 0.05)
    if field == "chromatography":
        c = TRUTH["chromatography"]
        main = np.exp(-0.5 * ((x - place) / c["width"]) ** 2)
        late = c["impurity"] * np.exp(-0.5 * ((x - place - c["after"] * c["width"]) / c["width"]) ** 2)          # an impurity that elutes later
        return height * (main + late) + base
    n = TRUTH["nmr"]
    width = n["width"] * (1.0 + n["growth"] * (level - n["at"]))          # the line grows by 12 per cent from 290 to 320 K
    return height / (1.0 + (2.0 * (x - place) / width) ** 2) + base


def in_the_units_of_the_fit(field: str) -> dict:
    """The truth of `field` in the parameters that the first model and the term that it lacks fit: for nmr, the width at
    the mean temperature and the a of along("width", "temperature"), whose z is the temperature centred and in units of
    its standard deviation over the records. For chromatography, the size of second(1.4 widths)."""
    if field == "chromatography":
        c = TRUTH["chromatography"]
        return {"width": c["width"], f"second component at {c['after'] * c['width']:g}": c["impurity"]}
    n, levels = TRUTH["nmr"], np.array(FIELDS["nmr"]["levels"])
    middle = n["width"] * (1.0 + n["growth"] * (levels.mean() - n["at"]))
    return {"width": float(middle), "width along temperature": float(n["width"] * n["growth"] * levels.std() / middle)}


def write(field: str, folder, seed: int = 11) -> list:
    """Writes the project of `field` into `folder`: traces/, description.toml, model.py and an empty terms/. A file that
    is already there is refused, so that no work of a researcher is written over."""
    if field not in FIELDS:
        raise ValueError(f"no example of {field!r}: the examples are {', '.join(FIELDS)}")
    f, home = FIELDS[field], Path(folder)
    for name in ("description.toml", "model.py", "traces", "terms"):
        if (home / name).exists():
            raise ValueError(f"{home / name} is there already: an example is written into a folder that does not hold one")
    (home / "traces").mkdir(parents=True)
    (home / "terms").mkdir()
    rng = np.random.default_rng(seed)
    x = np.round(f["start"] + np.arange(0.0, f["span"] + f["step"] / 2.0, f["step"]), 6)
    key, unit = f["covariate"]
    text = ["[problem]", f'name = "{field}"', f'what = "{f["what"]}"', f'abscissa = "{f["abscissa"][0]}"', f'abscissa_unit = "{f["abscissa"][1]}"',
            f'value = "{f["value"][0]}"', f'value_unit = "{f["value"][1]}"', "", f"[covariates.{key}]", f'unit = "{unit}"', ""]
    written = []
    for k, level in enumerate(f["levels"]):
        y = truth(field, x, level, k) + rng.normal(0.0, f["noise"], x.size)
        path = home / "traces" / f"{key}_{level:g}.csv"
        path.write_text(f"{f['abscissa'][0]},{f['value'][0]}\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
        written.append(path)
        text += ["[[trace]]", f'file = "traces/{path.name}"', f'columns = ["{f["abscissa"][0]}", "{f["value"][0]}"]', f"{key} = {level:g}", ""]
    text += ["[noise]", 'from = "value"', f"value = {f['noise']}", "", "[addenda]", 'model = "model.py"', "",
             "[addenda.lever.quieter]", 'kind = "noise"', "factor = 0.5", 'words = "a quieter detector, which halves the noise of each sample"', "",
             "[addenda.lever.repeats]", 'kind = "repeats"', "factor = 3", 'words = "each run made three times"', "",
             "[addenda.lever.slower]", 'kind = "samples"', "factor = 2", 'words = "twice the samples over the same span"', ""]
    (home / "description.toml").write_text("\n".join(text))
    (home / "model.py").write_text(MODELS[field])
    return written + [home / "description.toml", home / "model.py", home / "terms"]
