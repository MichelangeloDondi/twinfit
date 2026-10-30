"""An example of the package, finance: the daily log returns of a stock index in five periods of different volatility,
5000 returns each, counted in bins of 0.2 per cent within four times the volatility of the period, so that each record is
the histogram of the returns of one period. The returns follow a Student law of 4 degrees of freedom, whose standard
deviation is 0.9 of the volatility that the options of the period quoted, while the first model is a normal law, its
standard deviation a fraction of the quoted volatility of each period: what it lacks is the fat tails. The noise is that
of a counter, the root of the counts of each bin. Writes DIR with description.toml, model.py, traces/ and terms/. Its
truth is TRUTH, written nowhere in DIR."""
from pathlib import Path

import numpy as np
from scipy.special import gammaln

SETUP = ['[setup.index]', 'words = "daily closing prices of a stock index, as log returns in per cent, five periods of 5000 daily returns each, pooled over the stocks of the index"', 'known = "measured"', '',
         '[setup.options]', 'words = "the volatility of each period as the options on the index quoted it at its start, in per cent a day"', 'known = "measured"', '',
         '[setup.histogram]', 'words = "the returns within four times the quoted volatility counted in bins of 0.2 per cent, the noise of each bin stated as the root of its counts as smoothed over neighbouring bins, one at least"', 'known = "set"', '']

home = Path(home)  # noqa: F821  the folder of the project, given by twinfit.examples.gallery
(home / "traces").mkdir(parents=True)
(home / "terms").mkdir()
rng = np.random.default_rng(20261024)
BIN = 0.2          # the width of a bin of the histogram, per cent
returns, ratio, nu = 5000.0, 0.9, 4.0          # the returns of a period, the standard deviation of the returns over the quoted volatility, and the degrees of freedom of the law
volatilities = (0.6, 0.9, 1.3, 1.9, 2.8)          # quoted, per cent a day
text = ['[problem]', 'name = "daily returns"',
        'what = "the number of daily returns in each bin of 0.2 per cent, for five periods of 5000 returns of different quoted volatility. Each record has its own number of returns and its own mean return. The first model is a normal law whose standard deviation is one fraction of the quoted volatility of the period, the same for all five"',
        'abscissa = "log return"', 'abscissa_unit = "%"', 'value = "returns"', 'value_unit = "returns in the bin"', '', '[covariates.volatility]', 'unit = "% a day"', '']
for vol in volatilities:
    edge = np.floor(4.0 * vol / BIN) * BIN
    x = np.round(np.arange(-edge, edge + 1e-9, BIN), 3)          # the log return of a day, at the middle of each bin, per cent
    mean = rng.uniform(-0.03, 0.05)          # the mean return of the period, per cent a day
    c = ratio * vol / np.sqrt(nu / (nu - 2.0))          # the scale of the Student law that has the standard deviation ratio * vol
    t = (x - mean) / c
    pdf = np.exp(gammaln(0.5 * (nu + 1.0)) - gammaln(0.5 * nu) - 0.5 * np.log(nu * np.pi) - 0.5 * (nu + 1.0) * np.log1p(t ** 2 / nu)) / c
    expected = returns * BIN * pdf          # the returns expected in each bin
    sd = np.sqrt(np.maximum(expected, 1.0))
    y = counted(rng, expected, sd)          # noqa: F821  the counts of a counter, at the level and of the kind that the gallery gives, his message 194
    name = f"traces/v{int(round(10 * vol))}.csv"
    (home / name).write_text("return,returns,sd\n" + "".join(f"{a:g},{b:g},{stated(s):.6g}\n" for a, b, s in zip(x, y, sd)))          # noqa: F821
    text += ['[[trace]]', f'file = "{name}"', 'columns = ["return", "returns"]', f'volatility = {vol:g}', '']
text += SETUP + ['[noise]', 'from = "column"', 'column = "sd"', '', '[addenda]', 'model = "model.py"', '']
(home / "description.toml").write_text("\n".join(text))
(home / "model.py").write_text('''"""The first model: a normal law of the returns of a period, its standard deviation a fraction of the volatility that the options quoted, on its own number of returns and its own mean."""
import numpy as np

SHARED = ("volatility ratio",)              # the standard deviation of the returns over the volatility that the options quoted, the same in every period
OWN = ("returns", "mean")                      # of each record
START_SHARED = [0.9]
BOUNDS_SHARED = ([0.2], [3.0])
CENTRE = "mean"
NAME = "one normal law of the returns"
BIN = 0.2                                   # the width of a bin of the histogram, in per cent


def model(x, shared, own, design):
    sd = shared[0] * design["volatility"]
    return own[0] * BIN * np.exp(-0.5 * ((x - own[1]) / sd) ** 2) / (sd * np.sqrt(2.0 * np.pi))


def start_own(record):
    return [float(record.y.sum()), float(np.average(record.x, weights=np.maximum(record.y, 0.0)))]


def centre(shared, own, design):
    return own[1]
''')

TRUTH = "truth: the returns follow a Student law whose tail index, its degrees of freedom, is 4: fat tails, with a standard deviation 0.9 of the quoted volatility of the period"
