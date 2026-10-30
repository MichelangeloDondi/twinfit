# Decisions of the author and questions left open

This page is `docs/DECISIONS.md` of https://github.com/MichelangeloDondi/twinfit. See [`CONTINUING.md`](../twinfit/CONTINUING.md) for what is where, [`VALIDATION.md`](VALIDATION.md) and [`ROADMAP.md`](ROADMAP.md).

**Words used.** Section 1 of `twinfit/TUTORIAL.md` defines those of the package. A *draw* is one synthetic experiment of the toy, a model of the author's line with five terms (natural width, collisions, permeated gas, laser noise, transit), and a *count* runs the package on many draws. A *criterion* is a share that `tools/met.py` reads against a line of `toy/thresholds.toml`.

## Decisions taken

| decision, why and where | cost of changing it |
|---|---|
| **1. A general tool, rubidium its example:** no term of a field is named and the recorded traces are not needed. Absent: `toy/adopt.py`, `toy/noise.py`, `toy/build.py`, `toy/pages.py`, `toy/meeting.py`. Present but needing the absent `data/recorded/`: `toy/real.py`, `toy/noiseless.py make OUT --recorded`, the `his` units of `toy/chosen.py`. | None for counts on synthetic draws, but `runs/` holds no row of `toy/chosen.json`. |
| **2. Blind:** since 0.1.0rc16 the report names no term, and the package's list of generic terms, its dictionary, is read only on `--dictionary` (`twinfit/loop.py`). | A table of another mode meets no criterion. The dictionary has no count. |
| **3. Joint fits only,** the shared parameters being the quantities sought (`twinfit/fit.py`), **and the researcher writes each term** as `terms/NN.py` with `term` and `why`, six words or more (`twinfit/loop.py`). | A curve by curve mode is new code. Nothing tests `why` against the term. |
| **4. A theory value is held with its bar and freed only as a test** (`theory_moves` in `twinfit/best.py`): adding the transit term at theory moved the laser noise of the recorded traces from 1.148 MHz, fitted by the author outside this release under the model without that term, to 0.369 ± 0.028 MHz (`toy/config.toml`, whose comment keeps 0.455 ± 0.039, the value that the toy held before). | C7 scores the fit that frees the held collisions and transit. |
| **5. First order only where it agrees with the fit made again** within half a bar on 0.95 of the sets (`SHARE`, `SAME` in `twinfit/unsuspected.py`), and no bar for a parameter correlated beyond 0.999 with another, which the records fix only as a combination (`TOLD_APART` in `twinfit/best.py`). | A lower limit drops bars of strong correlations, a higher one prints bars the records do not fix. |
| **6. Two places tie** where their excesses differ by under 1.96 standard deviations of their difference: the records cannot tell them apart (`TIE` in `twinfit/where.py`). | C1 rests on it: in `toy/pointing.json` the place first alone is 173 of 400, 0.432 ± 0.025, with the tie 312, 0.780 ± 0.021. |
| **7. A term is kept** when its fit ends inside its bounds, the Bayesian criterion falls by more than 2 and the twin admits the larger model, so a term pays for its parameters (`TIED` in `twinfit/best.py`). The author has not ruled on the margin 2. | Every count goes stale (its `--check` says its table was made from other files), and C6 and C7, which run this rule, are counted again. |
| **8. numpy 2.5 is the floor** (`pyproject.toml`), for a speed gain measured on a Mac only. | By the author's reading the package calls only what Python 3.11, numpy 1.17 and scipy 1.4 hold. |
| **9. Every change is a new version,** with a digest, a hash, of the code and pages of the package folder but `version.py` (`check()` in `twinfit/version.py`), and an uploaded version is never replaced. A post release such as 0.1.0.post1 changes pages alone, but `check()` raises until `version.py` is stamped. | A count's digest is of the syntax tree of its code, without comments or docstrings, and of the text of the functions its program lists in `MAKES`, the gallery's (C11) of the bytes of code and data: pages move neither. The tree differs between Python versions, so the other counts check on Python 3.14 only. `tools/met.py` reads 0.1.0.post1 as 0.1.0, but a table of 0.1.1 meets nothing until `release` in `toy/thresholds.toml` moves. |

## Thresholds, and which were set after tables were read

Drafted on 2026-10-23 before any count at 0.1.0rc15, the lines of `toy/thresholds.toml` were written again on 2026-10-24, when the author delegated them. Set after tables were read:

- C6, C7 and C10, after a review read the loop count's progress log, and C7's band after the tables of 0.1.0rc16.
- C2's minimum draws, from 30 to 36.
- C5c's 0.12 repeats C5's, set from 30 sets (100 sets give 0.09), with the reading 0.188 ± 0.025 of 250 beside it, dated 2026-10-27, before any count of 0.1.1.

## Questions left open

Until answered, keep what 0.1.0 does.

**1. Thirty sets or 100?** The default of 30 (`sets` in `twinfit/loop.py`) has no recorded reason, and the author has not ruled on it. C5c reads 0.188 ± 0.025 of 250 at 30 sets against at most 0.12 (`toy/pointing.json`), 0.076 ± 0.017 of 250 at 100, from rows no table holds. These print 47 and 19:

```bash
cat runs/pointing_own_*_right_*.jsonl | grep -c '"beyond": \[\['
cat runs/pointing_own100_*_right_*.jsonl | grep -c '"beyond": \[\['
```

Each draw against its own twin, the place first or tied reads 0.765 ± 0.030 of 200 at 30 and 0.760 ± 0.021 of 400 at 100 (`runs/pointing_own_*`, `runs/pointing_own100_*`). A new default in `twinfit/loop.py` leaves only C6, C7 and C11 stale. The 30 also stands in `--sets` of `twinfit/__main__.py`, `OWN_SETS` of `toy/pointing.py` and `SETS` of `toy/chosen.py` and `toy/gallery_count.py`, and `OWN_SETS` makes all of `toy/pointing.json` stale (248 hours of waves in its `queue` block, 135 of them at 100 sets).

**2. What should 0.1.1 do about C1 on the collisions?** They read 0.440 ± 0.050 of 100 against at least 0.60 (`toy/pointing.json`), the pooled place being met. At twice and three times the size of the twin, the true value of a term in `toy/config.toml`, the core comes first or tied on 0.610 ± 0.049 and 0.910 ± 0.029 of 100 draws (`runs/pointing_2_hidden_under_1_at_{2,3}_of_the_size_of_the_twin_*.jsonl`, in no table and remade by nothing here), but the transit's place is wrong on 100 of 100 draws at twice and the laser noise's at three times. Options: keep C1 as 0.1.0 does, move the line or `tools/met.py` to the pooled reading, redefine C1 before counting, or add a statistic of the small term, a broadening of the whole line. Counted by the author toward 0.1.1, by programs that are not in this release, that statistic sees the collisions on 100 of 100 draws and places nothing, so C1 would have to score it as a place.

**3. In which order do the versions come?** `ROADMAP.md` gives one. The other puts the twin's noise first: the twin now draws each sample from a normal law of the description's deviation (`twinfit/TUTORIAL.md`, section 1), where that order would draw the noise of the records. An edit of `twinfit/twin.py` leaves every count stale.

**4. When is the agentic version validated?** The README quotes the best of three early runs, 3 of 4 terms found, where one scored run, counted by the author outside this release, found 2 of 4. Three runs now give a gate, not a rate, with no pass criterion written. Or validate at 0.2.0, a criterion first.

**5. What becomes of the manual path?** `twinfit/agent/MANUAL.md` documents the library entry outside the loop and never names `dictionary`, `log`, `reuse` or `unsuspected`. Drop it, rewrite it, or keep it as 0.1.0 does.

## The labels in the comments

The comments and `twinfit/version.py` cite the author's private records, not in this repository: *his message N*, *review N* (a finding by a reading agent), *round N of the cold tests*, *row MNNN* (a deliverable of the plan), *candidate N*, *rule N*, *decision N*, *lesson N*, and paths such as `PLAN.md`, `reviews/`, `runs/reports/` and `runs/pages/`. *His decision 13 of PLAN.md* in `toy/thresholds.toml` is question 1 here.
