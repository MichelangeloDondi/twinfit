"""An example of the package, from a project of the cold tests: the absorption of one line at four
temperatures, each a Gaussian dip of one width on its own depth, place and background, complete, recorded behind a
detector whose filter correlates the noise of neighbouring samples, 0.8 from one sample to the next, while the
description states the noise of each sample alone. The first model holds the line whole: what it lacks is in the noise,
correlated from sample to sample. Writes DIR with description.toml, model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.detector]', 'words = "a photodiode behind a low-pass filter of the amplifier, whose time constant spans several samples"', 'value = 3.0', 'unit = "samples"', 'known = "set"', '']          # the instruments and settings of the data collection, a table of the description

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261023)
x = np.round(np.arange(-10.0, 10.0001, 0.01), 3)          # wavelength less the centre of the scan, pm
temps = (280.0, 300.0, 320.0, 340.0)
sd = 0.002          # the spread of one sample, the correlation of neighbouring samples being RHO of the gallery, 0.8
text = ['[problem]', 'name = "absorption dips"',
        'what = "the transmission of a gas cell across one absorption line, a tunable laser swept across it, at four temperatures of the cell. Each record has its own depth, place and background. The first model is one Gaussian dip of one width for all four"',
        'abscissa = "wavelength offset"', 'abscissa_unit = "pm"', 'value = "transmission"', 'value_unit = "a fraction"', '', '[covariates.temperature]', 'unit = "K"', '']
for k, t in enumerate(temps):
    e = noise(rng, sd, x.size, "correlated")          # noqa: F821  a filter of the detector: each sample keeps 0.8 of the last, its spread the same, at the level and of the kind that the gallery gives
    y = 0.95 - (0.3 + 0.05 * k) * np.exp(-0.5 * ((x - rng.uniform(-0.05, 0.05)) / 1.2) ** 2) + e
    name = f"traces/t{int(t)}.csv"
    (home / name).write_text("offset,transmission\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["offset", "transmission"]', f'temperature = {t:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(sd):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Gaussian dip of one width, on its own depth, place and background."""
import numpy as np

SHARED = ("width",)                         # the standard deviation of the dip, in pm, the same at every temperature
OWN = ("depth", "place", "background")      # of each record
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.05], [10.0])
CENTRE = "place"
NAME = "one Gaussian dip of one width"


def model(x, shared, own, design):
    return own[2] - own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2)


def start_own(record):
    return [float(np.median(record.y) - record.y.min()), float(record.x[np.argmin(record.y)]), float(np.median(record.y))]


def centre(shared, own, design):
    return own[1]
''')
TRUTH = "truth: the model is whole, width 1.2 pm, and the noise is correlated from sample to sample, 0.8 from one to the next, of spread 0.002 as stated"
