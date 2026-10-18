"""An example of the package, from a project of the cold tests: one resonance driven by a laser at
five powers, each a Lorentzian line whose full width grows with the power as 1.0 MHz times the root of 1 + P / 50 mW,
the broadening by saturation, and whose height grows as P / (P + 50 mW), with three glitches in the record at 40 mW.
The first model has one width for all five. Writes DIR with description.toml, model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.laser]', 'words = "the laser that drives the resonance, at five powers read by a power meter"', 'known = "measured"', '', '[setup.detector]', 'words = "a photomultiplier counting the fluorescence, whose record at 40 mW was taken while a lamp of the room flickered"', 'known = "measured"', '']          # the instruments and settings of the data collection, a table of the description

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261016)
x = np.round(np.arange(-15.0, 15.0001, 0.02), 3)          # detuning, MHz
powers = (10.0, 20.0, 40.0, 80.0, 160.0)
text = ['[problem]', 'name = "driven resonance"',
        'what = "the fluorescence of a vapour across one resonance, a laser swept across it, at five powers of the laser. Each record has its own height, place and background. The first model is one Lorentzian line of one width for all five"',
        'abscissa = "detuning"', 'abscissa_unit = "MHz"', 'value = "fluorescence"', 'value_unit = "V"', '', '[covariates.power]', 'unit = "mW"', '']
for p in powers:
    width = 1.0 * np.sqrt(1.0 + p / 50.0)
    y = 0.05 + 1.5 * p / (p + 50.0) / (1.0 + (2.0 * (x - 0.3) / width) ** 2) + noise(rng, 0.004, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    if p == 40.0:
        for at in (-9.3, 2.1, 7.7):
            y[int(np.argmin(np.abs(x - at)))] += 0.08
    name = f"traces/p{int(p)}.csv"
    (home / name).write_text("detuning,fluorescence\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["detuning", "fluorescence"]', f'power = {p:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.004):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Lorentzian line of one width, on its own height, place and background."""
import numpy as np

SHARED = ("width",)                         # the full width of the line at half its height, in MHz, the same at every power
OWN = ("height", "place", "background")     # of each record
START_SHARED = [2.0]
BOUNDS_SHARED = ([0.05], [20.0])
CENTRE = "place"
NAME = "one Lorentzian line of one width"


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), float(record.x[np.argmax(record.y)]), base]


def centre(shared, own, design):
    return own[1]
''')
TRUTH = "truth: full width 1.0 MHz times the root of 1 + P / 50 mW, the broadening by saturation, and height 1.5 P / (P + 50 mW), place 0.3 MHz in every record, background 0.05, three glitches of +0.08 in the record at 40 mW, at -9.3, 2.1 and 7.7 MHz, noise 0.004"
