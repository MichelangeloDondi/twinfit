"""An example of the package, optics: the light at the output of a Michelson interferometer against the position of
its moving mirror, which a piezo scan drives, at five wavelengths of the source, each a cosine fringe on its own mean,
visibility and phase. The piezo does not respond in proportion to its voltage: its displacement is bent by a square of
the position, 6 nm off a linear scan at both ends of a scan of 4000 nm, so that the phase of each curve drifts by a
square of the position, 0.08 to 0.17 radians at the ends. The first model holds a displacement in proportion to the
scan, with one calibration of the piezo for all wavelengths: what it lacks is the bend of the scan. Writes DIR with
description.toml, model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.source]', 'words = "a source tuned to five wavelengths, read by a wavemeter to a part in a thousand"', 'known = "measured"', '',
         '[setup.piezo]', 'words = "a piezo stage that moves the mirror, driven by a voltage ramp, whose calibration was taken at one wavelength"', 'known = "from its maker"', '']

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261027)
x = np.round(np.arange(-2000.0, 2000.0001, 10.0), 1)          # the position of the mirror from the middle of its scan, nm
half = 2000.0          # half the span of the scan, nm
scale, bend = 1.012, 0.003          # the calibration of the piezo, and the bend of its response
wavelengths = (450.0, 532.0, 633.0, 780.0, 980.0)          # nm
text = ['[problem]', 'name = "interferometer fringes"',
        'what = "the intensity at the output of a Michelson interferometer against the position of its moving mirror, which a piezo scan drives, read at five wavelengths of the source. Each record has its own mean, visibility and phase. The first model is a cosine fringe of period half the wavelength, the mirror displacement in proportion to the scan with one calibration of the piezo for all five"',
        'abscissa = "mirror position"', 'abscissa_unit = "nm"', 'value = "intensity"', 'value_unit = "a.u."', '', '[covariates.wavelength]', 'unit = "nm"', '']
for k, lam in enumerate(wavelengths):
    mean, visibility, phase = 1.0 + 0.05 * k, rng.uniform(0.7, 0.9), rng.uniform(0.0, 2.0 * np.pi)
    d = scale * (x + bend * x ** 2 / half)          # the displacement that the mirror makes, nm
    y = mean * (1.0 + visibility * np.cos(4.0 * np.pi * d / lam + phase)) + noise(rng, 0.02, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/w{int(lam)}.csv"
    (home / name).write_text("position,intensity\n" + "".join(f"{a:g},{b:.5f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["position", "intensity"]', f'wavelength = {lam:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.02):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: a cosine fringe of period half the wavelength, on its own mean, visibility and phase, the mirror displacement in proportion to the scan."""
import numpy as np

SHARED = ("piezo scale",)                   # the displacement of the mirror for each nm of the scan, the same at every wavelength
OWN = ("mean", "visibility", "phase")       # of each record
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.9], [1.1])
NAME = "one cosine fringe of period half the wavelength"


def model(x, shared, own, design):
    return own[0] * (1.0 + own[1] * np.cos(4.0 * np.pi * shared[0] * x / design["wavelength"] + own[2]))


def start_own(record):
    k = 4.0 * np.pi / record.design["wavelength"]
    basis = np.column_stack([np.ones_like(record.x), np.cos(k * record.x), np.sin(k * record.x)])
    a, b, c = np.linalg.lstsq(basis, record.y, rcond=None)[0]
    return [float(a), float(np.hypot(b, c) / a), float(np.arctan2(-c, b))]
''')

TRUTH = "truth: the piezo moves the mirror 1.012 nm for each nm of the scan, and its displacement is bent by a square of the position, 1.012 (x + 0.003 x^2 / 2000 nm), 6 nm off a linear scan at both ends of the scan"
