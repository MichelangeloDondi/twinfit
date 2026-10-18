"""An example of the package, electrophoresis: the band of one fragment in a gel, read by its fluorescence along the lane
at four voltages, each a Gaussian band of one width on its own height, place and background, with a tail on the side of
slower migration, a copy of the band three times wider on that side alone, of size 0.15, which the
adsorption of the fragment on the gel leaves and which the first model lacks. Writes DIR with description.toml,
model.py, traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.gel]', 'words = "an agarose gel of one concentration, cast for the four runs together"', 'known = "set"', '',
         '[setup.imager]', 'words = "a fluorescence imager read along the lane, of a pixel far finer than the band"', 'known = "from its maker"', '']

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261019)
x = np.round(np.arange(-12.0, 12.0001, 0.02), 3)          # distance along the lane from the expected place of the band, mm
voltages = (60.0, 80.0, 100.0, 120.0)
width = 1.0          # the standard deviation of the band, mm
text = ['[problem]', 'name = "gel bands"',
        'what = "the fluorescence of one band of a fragment along the lane of an electrophoresis gel, run at four voltages. Each record has its own height, place and background. The first model is one Gaussian band of one width, the same at every voltage"',
        'abscissa = "distance along the lane"', 'abscissa_unit = "mm"', 'value = "fluorescence"', 'value_unit = "a.u."', '', '[covariates.voltage]', 'unit = "V"', '']
for k, v in enumerate(voltages):
    height, place, base = 1.0 + 0.2 * k, rng.uniform(-0.5, 0.5), 0.05
    d = x - place
    band = np.exp(-0.5 * (d / width) ** 2)
    tail = np.exp(-0.5 * (d / (3.0 * width)) ** 2) / 3.0 * 0.5 * (1.0 + np.tanh(-d / 0.02))          # the band three times wider, over three, on the side of slower migration alone, rising over one sample as the tail of the package does
    y = height * (band + 0.15 * tail) + base + noise(rng, 0.004, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/v{int(v)}.csv"
    (home / name).write_text("distance,fluorescence\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["distance", "fluorescence"]', f'voltage = {v:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.004):g}', '', '[addenda]', 'model = "model.py"', '']          # noqa: F821
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Gaussian band of one width, on its own height, place and background."""
import numpy as np

SHARED = ("width",)                         # the standard deviation of the band, in mm, the same at every voltage
OWN = ("height", "place", "background")     # of each record
START_SHARED = [1.5]
BOUNDS_SHARED = ([0.05], [10.0])
CENTRE = "place"
NAME = "one Gaussian band of one width"


def model(x, shared, own, design):
    return own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), float(record.x[np.argmax(record.y)]), base]


def centre(shared, own, design):
    return own[1]
''')

TRUTH = "truth: a tail on the side of slower migration, a copy of the band three times wider below its place, of size 0.15 as the tail of the package counts it, beside a Gaussian band of standard deviation 1 mm"
