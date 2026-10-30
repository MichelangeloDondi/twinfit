# twinfit

Where a model lacks a term: the structure of the residuals of a fit, read against its twin.

twinfit is a Python package for any researcher who fits one model to a set of measured curves. It fits the model to all the curves at once, reads what the fit leaves through families of statistics, and sets each against a twin, synthetic curves drawn from the fitted model with the law of noise that the records are given and read in the same way. Its report says where the model misses, of which kind and along which setting of the experiment. It names no term of any field: the researcher writes the term, and the package keeps it only where the records ask for it.

```bash
pip install twinfit
```

The page of the package, with five minutes of a first run, is [`twinfit/README.md`](twinfit/README.md), which is also its page on PyPI.

## What this repository holds

| folder | what it holds |
|---|---|
| [`twinfit/`](twinfit/) | the package, as PyPI ships it, with its tutorial, its gallery of other fields and the card of an agent |
| [`toy/`](toy/README.md) | the bench of the package: a model of five terms of the author's own line, the two photon line of rubidium 5S to 6S, from which synthetic traces with a known truth are drawn, the programs that count what the package does on them, and their tables |
| [`runs/`](runs/) | the kept row of each draw of three of the counts behind the shares of the README, the place and the setting, the right models and the gallery, so that a table can be counted again and a draw made again |
| [`tools/met.py`](tools/met.py) | the reader of every criterion of `toy/thresholds.toml` on the tables |
| [`docs/`](docs/) | [`VALIDATION.md`](docs/VALIDATION.md), what the package was measured to do, [`ROADMAP.md`](docs/ROADMAP.md), what comes next, and [`DECISIONS.md`](docs/DECISIONS.md), the decisions that shape it and the questions left open |
| [`data/noise/`](data/noise/) | the law of the noise of each recorded trace of the author, read from them: no sample of a trace |

The recorded traces of the author are not here. The package does not need them, and the bench draws its traces from the toy model.

## Taking the work over

[`twinfit/CONTINUING.md`](twinfit/CONTINUING.md) is the page for a researcher who continues this work: what is where, what to run in the first hour, and where each line of the work starts. In an installed package it stands beside the README, in the folder that `python -m twinfit guide readme` names.

## Licence

MIT, [`LICENSE`](LICENSE).
