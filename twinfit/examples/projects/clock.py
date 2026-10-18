"""An example of the package, an atomic clock: the resonance of the clock transition read by a probe laser at five powers,
each a Lorentzian of one width on its own height and background, while the centre of the line moves with the power of
the probe, the light shift, 0.025 Hz for each microwatt, from the unperturbed frequency of the clock. The first model
holds one centre for every power, which is what the clock measures. Writes DIR with description.toml, model.py,
traces/ and terms/. Its truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np

SETUP = ['[setup.probe]', 'words = "the probe laser, at five powers read by a photodiode calibrated to 2 per cent"', 'known = "measured"', '',
         '[setup.reference]', 'words = "the frequency reference that the sweep is read against, of a stability far below the width of the line"', 'known = "from its maker"', '']

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261018)
x = np.round(np.arange(-6.0, 6.0001, 0.01), 3)          # detuning from the nominal frequency of the clock, Hz
powers = (2.0, 4.0, 6.0, 8.0, 10.0)          # microwatts
text = ['[problem]', 'name = "clock resonance"',
        'what = "the excitation of the clock transition of trapped atoms against the detuning of the probe laser, read at five powers of the probe. Each record has its own height and background. The first model is one Lorentzian line of one width and one centre, the frequency of the clock, the same at every power"',
        'abscissa = "detuning"', 'abscissa_unit = "Hz"', 'value = "excitation"', 'value_unit = "a fraction"', '', '[covariates.power]', 'unit = "uW"', '']
for p in powers:
    centre = 0.12 + 0.025 * p
    y = 0.02 + 0.8 / (1.0 + (2.0 * (x - centre) / 1.5) ** 2) + noise(rng, 0.003, x.size)          # noqa: F821  the noise at the level and of the kind that the gallery gives, his message 194
    name = f"traces/p{int(p)}.csv"
    (home / name).write_text("detuning,excitation\n" + "".join(f"{a:g},{b:.6f}\n" for a, b in zip(x, y)))
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["detuning", "excitation"]', f'power = {p:g}', '']
text += SETUP + ['[noise]', 'from = "value"', f'value = {stated(0.003):g}', '', '[addenda]', 'model = "model.py"', '',          # noqa: F821
                 '[addenda.ansatz."power scale"]', 'value = 1.0', 'uncertainty = 0.02', '']          # the calibration of the probe's power, 2 per cent, held: a light shift read along the power carries it (round 21 of the tests of the tutorial)
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: one Lorentzian line of one width and one centre, the frequency of the clock, on its own height and background."""
import numpy as np

SHARED = ("width", "centre")                # the full width at half maximum, in Hz, and the centre, the frequency of the clock less the nominal one, in Hz
OWN = ("height", "background")              # of each record
START_SHARED = [1.0, 0.0]
BOUNDS_SHARED = ([0.05, -3.0], [10.0, 3.0])
HELD = {"power scale": 1.0}                # the powers read times this, calibrated within 2 per cent: description.toml holds it, and a model whose centre moves with the power reads it
NAME = "one Lorentzian line of one width and one centre"


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - shared[1]) / shared[0]) ** 2) + own[1]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), base]

''')

TRUTH = "truth: the centre moves with the power of the probe, 0.12 Hz from the nominal frequency at no power and 0.025 Hz for each microwatt, the light shift, and the width is 1.5 Hz at every power"
