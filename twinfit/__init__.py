"""twinfit: where a model lacks a term. A model is fitted to curves recorded across a design, the structure of what
the fit leaves is read curve by curve with families of statistics, of the curve, of its values and of their
transformations, and how that structure changes along each setting of the design is set against a twin of the fit.
The package says where a term is missing, of which kind and along which setting, fits the terms that a researcher or
an agent writes, and gives the best model with the uncertainty of each parameter and three proposals for the model and
three for the experiment at each iteration. It knows nothing of any field.

The path of a researcher:

    description.py the description of the data collection, its traces, their conditions and the optional addenda
    loop.py        the loop: each model fitted, read against its twin, reported, and each term file fitted and kept or
                   refused. `python -m twinfit loop DIR`
    where.py       the report of where to search: the place, the kind, the setting and the family that sees it
    proposals.py   three proposals for the model and three for the experiment, each measured on the twin
    best.py        the best model after the iterations, and five parts of the uncertainty of each parameter
    unsuspected.py what can be said of the known parameters before the missing term is known
    agent/         the protocol for an agent that drives the package, and what a project needs to start
    examples/      worked projects of other fields, and the gallery, EXAMPLES.md: projects with their truth, what to do and
                   what not to do on each, and what the loop gave. `python -m twinfit example NAME DIR`
    __main__.py    the command line, `python -m twinfit`, and TUTORIAL.md, from nothing to the best model, written by tutorial.py

What lies under it:

    problem.py     records across a design, a working model with shared and own parameters
    fit.py         the canonical fit, weighted least squares of every record at once, within bounds if the problem gives any
    shapes.py      a curve from its characteristic function: a convolution or a mixture that has no formula
    linear.py      the linear picture of a problem at a fit, the own parameters taken out
    ledger.py      first order: what a missing term does to the estimates and what it leaves in the residuals
    identify.py    what the design can fix, and which block each estimate rests on
    families.py    families of statistics of a curve: moments, windows along both axes, the two sides, the abscissa against
                   the value, a transformation of the value. What each fixes by itself, what it reads in the residuals,
                   its estimate set against the fit's, and a combination of estimates free of the bias of named terms
    losses.py      chi functions: the fit under a cost of the residual that is not its square, with its first order factor, and
                   the windows along both axes of the residuals under each cost
    matched.py     statistics that are formed on the noisy values of a record, and matched to synthetic records of the model
    candidates.py  a dictionary of terms that a model may lack, each as a discrepancy of the mean
    propose.py     the candidates ranked on the residuals, with the bias that each would mean and what cannot be told apart
    twin.py        synthetic records from a model, to give every reading its law before it is read on the records
    mechanisms.py  a guessed term as the working model with one parameter more, or held at the size that theory gives
    refine.py      the candidates of a stage fitted, the law of their gains, the default choice, a candidate on synthetic records
    criteria.py    the information criteria, the level of a gain in a family of candidates, and the rule for a term on theory
    report.py      the report of a stage, in five parts
    noise.py       the weights of a fit from a law of the noise, and what the bars of that law do to the chi square
    pool.py        a function over a list of jobs in several forked processes, for the twin and for the law of the gains
    held.py        a term whose number is read from the dictionary HELD of the model at every call, so that a value held from
                   theory carries its bar
    readings.py    readings of the residuals by position, pairs of slices about the centre of each curve
    plane.py       the plane of position and residual value, cut on both axes, as one estimator
    standard.py    the standard checks of the residuals, beside the readings of the package
    notebook.py    the record of an analysis: its steps, its predictions stamped before the runs, its refusals
    version.py     the number of the package, and for each version the digest of its code and what it changed
    notebook.py    the record of an analysis: its steps, its stamped predictions, and what it refuses
    agent/         the protocol for an agent that drives the package, and what a project needs to start

The procedure and its limits are in agent/AGENT.md. This is a first version. Its table of candidates is first order, and a
candidate that is not small must be fitted to be judged: refine.py does that.
"""
import sys

if sys.version_info < (3, 11):          # the description is read with tomllib, of Python 3.11
    raise ImportError(f"twinfit needs Python 3.11 or later, and this is Python {sys.version_info.major}.{sys.version_info.minor}")
try:
    import numpy  # noqa: F401
    import scipy  # noqa: F401
except ImportError as missing:          # said in plain words, before any module of the package fails on it
    raise ImportError(f"twinfit needs numpy and scipy, and this Python lacks {missing.name}: pip install numpy scipy") from None

from .fit import Fit, fit  # noqa: F401,E402
from .ledger import Ledger  # noqa: F401,E402
from .linear import Linear  # noqa: F401,E402
from .plane import PlaneFit, columns, plane_fit  # noqa: F401,E402
from .problem import Problem, Record, even  # noqa: F401,E402

from .version import __version__  # noqa: E402,F401
