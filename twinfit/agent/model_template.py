"""The first model of this project, which `python -m twinfit loop` reads. Fill the places marked TO FILL.

Start with the small model that is understood in every part. The loop fits it to every curve at once, reads what it
leaves against its twin, and says where a term is missing: each term that you then write in terms/ is fitted beside it.
"""
import numpy as np

SHARED = ("width",)                     # TO FILL: the parameters that every curve shares, which the analysis is after
OWN = ("height", "place", "baseline")   # TO FILL: the parameters of each curve alone
START_SHARED = [1.0]                    # TO FILL: a start for each shared parameter, in the order of SHARED
BOUNDS_SHARED = ([1e-9], [1e6])         # optional: the lower bounds and the upper bounds of the shared parameters, a width above zero since a curve divides by it
CENTRE = "place"                        # an own parameter, which `centre` below goes with: the report reads the core and the wings of a curve about the place that `centre` returns. Where the place is shared, as one that a setting moves, name any own parameter here, as the height, and let `centre` return the place from `shared`
NAME = "one peak of one width"          # TO FILL: the model in a few words, as the report names it
EQUATION = ""                           # the model in LaTeX between dollar signs, as r"$h\,e^{-(x - x_0)^2/(2 w^2)} + b$", which the visual log typesets at each turn, `python -m twinfit log`. Left empty, the log writes the Python line that `model` returns
LETTERS = {}                            # each letter of EQUATION with its meaning, as {"h": "the height", "x_0": "the place", "w": "the width", "b": "the baseline"}, which the log writes under it


def model(x, shared, own, design):
    """TO FILL: the mean of one curve. `x` is its abscissa, `shared` and `own` hold the parameters in the order of
    SHARED and OWN, and `design` the settings of the curve, as the description gives them."""
    return own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    """TO FILL: a start for the own parameters of one curve, read from its values, in the order of OWN. This one reads a
    peak or a dip on a baseline that may slope: the baseline as a straight line through the means of the first and the
    last twentieth of the curve, the place as the first sample from the left whose distance from that line reaches four
    fifths of the largest, that distance as the height, negative for a dip, and the line there as the baseline. So every
    curve of a doublet starts on its lower line, and a dip on a slope is not taken for an end of the slope, where the
    median as a baseline took the end of a slope as steep as the dip. A curve started away
    from its peak or its dip is moved by no shared parameter, and the first fit is then refused."""
    order = np.argsort(record.x)
    x, y = np.asarray(record.x, dtype=float)[order], np.asarray(record.y, dtype=float)[order]
    n = max(len(y) // 20, 1)
    left, right = float(np.mean(y[:n])), float(np.mean(y[-n:]))
    line = left + (right - left) * (x - x[0]) / ((x[-1] - x[0]) or 1.0)
    up = np.abs(y - line)
    top = int(np.flatnonzero(up >= 0.8 * up.max())[0])
    return [float(y[top] - line[top]), float(x[top]), float(line[top])]


def start_shared(records):
    """optional, and where it is defined it takes the place of START_SHARED: a start for each shared parameter, read
    from the curves. This one reads the width of the peak or the dip of each curve where its distance from the straight
    line through its ends is half the largest, and gives their mean as a standard deviation. Change it or remove it when
    SHARED changes."""
    widths = []
    for record in records:
        order = np.argsort(record.x)
        x, y = np.asarray(record.x, dtype=float)[order], np.asarray(record.y, dtype=float)[order]
        n = max(len(y) // 20, 1)
        left, right = float(np.mean(y[:n])), float(np.mean(y[-n:]))
        line = left + (right - left) * (x - x[0]) / ((x[-1] - x[0]) or 1.0)
        up = np.abs(y - line)
        above = x[up >= 0.5 * up.max()] if up.max() > 0 else x
        widths.append((above.max() - above.min()) / 2.3548 if above.size > 1 and above.max() > above.min() else float(np.ptp(x)) / 10.0)
    return [float(np.mean(widths))]


def centre(shared, own, design):
    return own[1]
