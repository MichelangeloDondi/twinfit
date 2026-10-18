"""An example of the package, from a project of the cold tests: the photoionization yield across one
autoionizing resonance at four pressures, each a Fano profile of one width and one asymmetry q = 3 on its own height,
place and background, while the first model is a Lorentzian line: what it lacks is the interference that makes the line
asymmetric, a shape that no term of the package holds. Writes DIR with description.toml, model.py, traces/ and terms/.
Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.source]', 'words = "a monochromator swept across the resonance, of a resolution far finer than the line"', 'known = "from its maker"', '']          # the instruments and settings of the data collection, a table of the description

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261024)
x = np.round(np.arange(-20.0, 20.0001, 0.02), 3)          # photon energy less the centre of the scan, meV
pressures = (1.0, 2.0, 3.0, 4.0)
q, gamma = 3.0, 2.0
text = ['[problem]', 'name = "resonance yield"',
        'what = "the ion yield across one resonance of an atomic gas, the photon energy swept across it, at four pressures of the gas. Each record has its own height, place and background. The first model is one Lorentzian line of one half width for all four"',
        'abscissa = "energy offset"', 'abscissa_unit = "meV"', 'value = "ion counts"', 'value_unit = "counts, scaled"', '', '[covariates.pressure]', 'unit = "mbar"', '']
for k, p in enumerate(pressures):
    e = (x - rng.uniform(-0.1, 0.1)) / gamma
    shape = ((q + e) ** 2 / (1.0 + e ** 2)) / (q * q + 1.0)          # the Fano profile, its top near one
    y = (0.5 + 0.2 * k) * shape + 0.1 + noise(rng, 0.004, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/p{int(p)}.csv"
    (home / name).write_text("energy,counts\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["energy", "counts"]', f'pressure = {p:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.004):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Lorentzian line of one half width, on its own height, place and background."""
import numpy as np

SHARED = ("half_width",)                    # the half width at half maximum of the line, in meV, the same at every pressure
OWN = ("height", "place", "background")     # of each record
START_SHARED = [1.5]
BOUNDS_SHARED = ([0.05], [20.0])
CENTRE = "place"
NAME = "one Lorentzian line of one half width"


def model(x, shared, own, design):
    return own[0] / (1.0 + ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    top = int(np.argmax(np.abs(record.y - base)))
    return [float(record.y[top] - base), float(record.x[top]), base]


def centre(shared, own, design):
    return own[1]
''')
TRUTH = "truth: a Fano profile, q = 3 and half width 2 meV at every pressure, its top near one, and the first model is a Lorentzian, which lacks the interference"
