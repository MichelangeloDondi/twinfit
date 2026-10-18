"""A project for the agentic cold tester of round 2, whose truth it is not told: absorption lines of a gas cell at five
pressures, each a Gaussian of one width on its own depth, place and baseline, with a satellite line 12 MHz above each, of 8 per cent of its area,
which the first model lacks. Writes DIR with description.toml, model.py, traces/ and
terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.laser]', 'words = "a diode laser swept across the line"', 'known = "measured"', '', '[setup.gas]', 'words = "the gas of the cell, whose isotopes and their lines were not checked"', 'known = "not recorded"', '']          # the instruments and settings of the data collection, a table of the description

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261017)
x = np.round(np.arange(-40.0, 40.0001, 0.1), 4)          # detuning, MHz
pressures = (10.0, 20.0, 30.0, 40.0, 50.0)
text = ['[problem]', 'name = "gas cell"',
        'what = "the absorption of one line of a gas in a cell, measured at five pressures by a laser swept across it. Each record has its own depth, centre and baseline. The first model is one Gaussian line of one width for all five"',
        'abscissa = "detuning"', 'abscissa_unit = "MHz"', 'value = "absorbance"', 'value_unit = "a.u."', '', '[covariates.pressure]', 'unit = "Pa"', '']
for k, p in enumerate(pressures):
    depth, place, base, width = 0.5 + 0.05 * k, rng.uniform(-1.0, 1.0), 0.02, 4.0
    core = np.exp(-0.5 * ((x - place) / width) ** 2)
    wide = np.exp(-0.5 * ((x - place - 12.0) / width) ** 2)
    y = depth * (core + 0.08 * wide) + base + noise(rng, 0.002, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/p{int(p)}.csv"
    (home / name).write_text("detuning,absorbance\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["detuning", "absorbance"]', f'pressure = {p:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.002):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Gaussian line of one width, on its own depth, centre and baseline."""
import numpy as np

SHARED = ("width",)                     # the standard deviation of the line, in MHz, the same at every pressure
OWN = ("depth", "centre", "baseline")   # of each record
START_SHARED = [3.0]
BOUNDS_SHARED = ([0.1], [30.0])
CENTRE = "centre"
NAME = "one Gaussian line of one width"


def model(x, shared, own, design):
    return own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    return [float(record.y.max() - record.y.min()), float(record.x[np.argmax(record.y)]), float(np.median(record.y))]


def centre(shared, own, design):
    return own[1]
''')
TRUTH = "truth: a second component 12 MHz above, of relative area 0.08, beside a Gaussian of width 4 MHz"
