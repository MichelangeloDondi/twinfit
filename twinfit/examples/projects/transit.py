"""An example of the package, the search for exoplanets: the light of a star across one transit of its planet, recorded
through four filters, each a dip of one half duration on its own depth, centre and baseline. The star is darker at its
limb than at its centre, by a coefficient u that falls with the wavelength, 0.6 in the blue and 0.3 in the near
infrared, so that the bottom of each dip is round, the rounder the bluer the filter. The first model is a dip of a
uniform star, flat at its bottom: what it lacks is the darkening of the limb. Writes DIR with description.toml,
model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.telescope]', 'words = "a telescope of one aperture, with a filter wheel of four filters at 450, 550, 650 and 800 nm"', 'known = "from its maker"', '',
         '[setup.star]', 'words = "the host star, of one spectral type, whose darkening of the limb was not measured"', 'known = "not recorded"', '']

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261020)
x = np.round(np.arange(-3.0, 3.0001, 0.005), 3)          # time from the middle of the transit, hours
half, edge = 1.2, 0.05          # the half duration, and the time over which the planet enters the disc of the star, hours
filters = (450.0, 550.0, 650.0, 800.0)          # nm
text = ['[problem]', 'name = "planet transit"',
        'what = "the relative flux of a star across one transit of its planet, recorded through four filters. Each record has its own depth, centre and baseline. The first model is a dip of one half duration, whose ingress takes 0.05 hours, flat at its bottom, as a star of uniform brightness gives it"',
        'abscissa = "time"', 'abscissa_unit = "h"', 'value = "relative flux"', 'value_unit = "1"', '', '[covariates.wavelength]', 'unit = "nm"', '']
for lam in filters:
    u = 0.6 - 0.3 * (lam - 450.0) / 350.0
    depth, centre, base = 0.010, rng.uniform(-0.02, 0.02), 1.0
    t = x - centre
    inside = 1.0 / (1.0 + np.exp(-(half - np.abs(t)) / edge))
    mu = np.sqrt(np.clip(1.0 - (t / half) ** 2, 0.0, 1.0))
    y = base - depth * inside * (1.0 - u * (1.0 - mu)) / (1.0 - u / 3.0) + noise(rng, 0.0004, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/f{int(lam)}.csv"
    (home / name).write_text("time,flux\n" + "".join(f"{a:g},{b:.7f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["time", "flux"]', f'wavelength = {lam:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.0004):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: the dip of a star of uniform brightness, of one half duration, on its own depth, centre and baseline."""
import numpy as np

SHARED = ("half_duration",)                 # half the duration of the transit, in hours, the same in every filter
OWN = ("depth", "centre", "baseline")       # of each record
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.2], [2.5])
CENTRE = "centre"
NAME = "the dip of a uniform star"
EDGE = 0.05                                 # the time over which the planet enters the disc of the star, in hours


def model(x, shared, own, design):
    inside = 1.0 / (1.0 + np.exp(-(shared[0] - np.abs(x - own[1])) / EDGE))
    return own[2] - own[0] * inside


def start_own(record):
    base = float(np.median(record.y))
    return [float(base - record.y.min()), float(record.x[np.argmin(record.y)]), base]


def centre(shared, own, design):
    return own[1]
''')

TRUTH = "truth: a star darker at its limb, by a coefficient u = 0.6 - 0.3 (wavelength - 450 nm) / 350 nm, 0.6 in the blue and 0.3 at 800 nm, under a planet of depth 0.010 and half duration 1.2 h"
