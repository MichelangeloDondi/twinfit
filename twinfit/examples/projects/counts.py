"""An example of the package, from a project of the cold tests: photon counts of one emission line at
four powers of the pump, each a Lorentzian of one half width on its own height, place and background, complete, drawn
with the noise of a counter, the root of the counts, while the description states one noise for every sample, that of
the background. The first model holds the line whole: what it lacks is in the noise, which grows with the signal.
Writes DIR with description.toml, model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.detector]', 'words = "a single photon counter, whose counts in each bin of the sweep are recorded"', 'known = "set"', '']          # the instruments and settings of the data collection, a table of the description

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261022)
x = np.round(np.arange(-15.0, 15.0001, 0.02), 3)          # detuning, GHz
powers = (1.0, 2.0, 3.0, 4.0)
text = ['[problem]', 'name = "photon counts"',
        'what = "the photons counted from one emission line of an atomic vapour, the frequency of a probe laser swept across it, at four powers of a pump. Each record has its own height, place and background. The first model is one Lorentzian line of one half width for all four"',
        'abscissa = "detuning"', 'abscissa_unit = "GHz"', 'value = "counts"', 'value_unit = "photons in each bin"', '', '[covariates.power]', 'unit = "mW"', '']
for k, p in enumerate(powers):
    mean = 400.0 * p / (1.0 + ((x - rng.uniform(-0.05, 0.05)) / 1.5) ** 2) + 40.0
    y = counted(rng, mean, 6.3)          # noqa: F821  the counts of a counter, at the level and of the kind that the gallery gives, his message 194
    name = f"traces/p{int(p)}.csv"
    (home / name).write_text("detuning,counts\n" + "".join(f"{a:g},{b:g}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["detuning", "counts"]', f'power = {p:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(6.3):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Lorentzian line of one half width, on its own height, place and background."""
SHARED = ("half_width",)                    # the half width at half maximum of the line, in GHz, the same at every power
OWN = ("height", "place", "background")     # of each record
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.05], [20.0])
CENTRE = "place"
NAME = "one Lorentzian line of one half width"


def model(x, shared, own, design):
    return own[0] / (1.0 + ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    import numpy as np
    return [float(record.y.max() - record.y.min()), float(record.x[np.argmax(record.y)]), float(np.median(record.y[:50]))]


def centre(shared, own, design):
    return own[1]
''')
TRUTH = "truth: the model is whole, half width 1.5 GHz, and the noise is that of a counter, the root of the counts, 6.3 on the background and 40 at the top of the strongest line"
