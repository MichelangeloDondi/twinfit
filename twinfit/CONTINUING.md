# Continuing this work

This page is written for a researcher who takes twinfit over from its author. It states what the package is, where each of its parts is, what was measured, what to run in the first hour, and where each line of work starts. It builds on the README of the package. Section 1 of `TUTORIAL.md` defines every term of the package.

## The package and its example

twinfit is a general purpose tool. A researcher fits one model to several curves, for example spectra at several temperatures. twinfit reads the residuals of that fit through families of statistics and compares each statistic with its value on a twin, that is, on synthetic curves drawn from the fitted model with the law of noise assigned to the records, fitted and read in the same way. It then states where the model misses, of which kind the miss is, and along which setting it varies. The physics enters after that reading: the researcher writes the term that accounts for it.

The worked example, and the bench on which the package was measured, is a line of the author's own experiment. A model of that line in five terms draws synthetic curves whose truth is known. The figures of the README were counted on those curves, and on a gallery of eleven synthetic projects of other fields. Of the five values of the model, two were fitted to the author's recorded curves, a constant width, 0.037 ± 0.015, and the spread of the frequency noise of the light source, 0.369 ± 0.028, both in the unit of the model's widths (`toy/config.toml`), two are held at their theoretical values with their uncertainties, and one is taken from the literature. The recorded curves are private, and the tool needs none of them. The model, its terms, their values and the source of each are given in `toy/README.md` of the repository. The condition of each value, the fit that gave it and its uncertainty are given in the `source` fields of `toy/config.toml`, `toy/const.toml` and `toy/twin.toml`.

## Where each part is

The wheel and the source distribution on PyPI carry the package and its pages: this page, the README, `TUTORIAL.md`, `EXAMPLES.md` (eleven projects of other fields), `RUBIDIUM.md` (the author's own case, as plots), the card of an agent in `agent/`, and `version.py`, which lists every version with its digest and its changes. They do not carry the bench.

The repository, https://github.com/MichelangeloDondi/twinfit, carries the rest.

| folder | contents | use by a successor |
|---|---|---|
| `twinfit/` | the package | use and change |
| `toy/` | the model of five terms, the programs of the counts, their tables, the thresholds of the criteria, `again.py` and `count.py` | checks, one draw made again, counts made again |
| `runs/` | the kept row of each draw of the counts of pointing, power and the gallery, and the twin of each of their cases | reading, recount of a table |
| `tools/met.py` | the reader of every criterion of `toy/thresholds.toml` on the tables | reading of the criteria |
| `docs/` | `VALIDATION.md`, `ROADMAP.md`, `DECISIONS.md` | reading before any change |
| `data/noise/` | the law of the noise of each recorded curve, one for each peak and condition: the variance of a sample, a² + b h + c h² + d s², in the height h of the curve at that point and its slope s. No sample of any curve | the source of the noise of the model's draws |

The recorded curves and the author's working records are not public. The comments of the code refer to entries of those records ("his message N", "review N", "row MNNN", "candidate N", "decision N", "lesson N") and to paths such as `PLAN.md`, `reviews/` and `runs/reports/`. None of them is needed to run anything. `docs/DECISIONS.md` gives the reasons that matter.

## What was measured

`docs/VALIDATION.md` lists every criterion, written from the tables by `docs/validation.py`. Two are not met at this release.

| criterion | share | threshold |
|---|---|---|
| C1, the place of a width that grows steeply along a setting, put first where the width is hidden at the size of the model | 0.440 ± 0.050 of 100 draws | at least 0.60 |
| C5c, a false alarm on a right model at the default of 30 sets | 0.188 ± 0.025 of 250 draws | at most 0.12 |

The other criteria are met, among them the place of a hidden term pooled over four terms, 0.780 ± 0.021 of 400 draws, and the eleven projects of other fields, 0.759 ± 0.029 of 220 draws. Not yet validated: the agentic version, use of the package by anyone but the author, and terms beyond the fifth of the model.

## The first hour

The package runs on Python 3.12 or later. The bench needs Python 3.14, because the digest of a count hashes the syntax tree of its code, which differs between versions of Python. On another version each check but that of the gallery reports that its table was made from other files.

1. Install the package with matplotlib, which the tutorial and the figures need: `pip install "twinfit[log]"`.
2. Run the five minutes of the README, the clock at five probe powers.
3. Run every command of the tutorial again, on a quiet machine, since a loaded one can make the check refuse. It takes one to five minutes:

```bash
python -m twinfit.tutorial --check
```

It ends with "the tutorial: every output of TUTORIAL.md is what its command prints now".

4. Clone the repository and run its checks, on Python 3.14:

```bash
git clone https://github.com/MichelangeloDondi/twinfit
cd twinfit
python toy/pointing.py --check
python toy/power.py --check
python toy/gallery_count.py --check
python toy/draws.py --check
python tools/met.py M456,M505,M457,M458,M441,M444,M459 --counted
python toy/again.py "2 hidden under 1" 7
```

Each of the first three commands compares the digest of its code with that of its table, sums again the rows that the table holds, and reports that its table was counted, or made, "from the files as they are". A pass means that the table matches the code. The kept rows of `runs/` are read by `toy/again.py` and by `toy/count.py`. The fourth command counts the kept fits of each model on its 100 synthetic draws, makes three of them again, and ends with "every model passes, and three fits of each were made again". `tools/met.py` prints the readings of `docs/VALIDATION.md`. Each label of its command is a row of the author's plan and names the criteria whose `row` field in `toy/thresholds.toml` holds it. `toy/again.py` makes draw 7 of that case again from its seed, in seconds, and ends with "the kept row, field for field".

## Three lines of work

- **Use.** The README, `TUTORIAL.md`, and, for an agent, the card of `agent/AGENT.md`.
- **Test.** `toy/README.md` states how each count is made, by `python toy/count.py PROGRAM`, what it costs, and how a draw is made again.
- **Development.** `docs/ROADMAP.md` gives the versions that follow and a first step that takes an hour. `docs/DECISIONS.md` gives the decisions taken and the questions left open.

## Three rules of the work

1. A change of a file of `twinfit/` is a new version. `twinfit.version.check()` raises until `version.py` holds the new number, its digest and what it changed. A post release, as this one, changes pages alone and keeps the tables of its release.
2. A count is made again after a change of the code that it runs. Its `--check` says when.
3. A criterion is read with `tools/met.py`, against thresholds that the author drafted before the counts. The header of `toy/thresholds.toml` and `docs/DECISIONS.md` say which were set again after tables were read.
