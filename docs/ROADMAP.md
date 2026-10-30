# What comes after 0.1.0.post1

This page builds on [`CONTINUING.md`](../twinfit/CONTINUING.md), [`VALIDATION.md`](VALIDATION.md) and [`DECISIONS.md`](DECISIONS.md). Every command on this page is run from the root of the repository, on Python 3.14.6. Three words are used throughout. A **draw** is one synthetic data set of the model of `toy/`. The **sets** of a draw, 30 by default, are the synthetic copies drawn from the fitted model, that is, its twin. The **line** of the report is the largest excess, over all the answers together, that a set of a right model, one that lacks nothing, reaches on 95 of 100 sets (`twinfit/TUTORIAL.md`). The line of a single statistic, in section 3.1, is its own 95th percentile over the sets.

## 1. Where 0.1.0.post1 stands

The code is that of 0.1.0. The pages are new. `python -m twinfit --version` prints the number of this release, and `twinfit.version.check()` returns it, because its digest, which covers the pages as well as the code, is stamped in `twinfit/version.py`.

Two criteria are not met (`docs/VALIDATION.md`). C1 on the collisions, the share of draws on which the report puts their place first, is 0.440 ± 0.050 of 100 draws against a threshold of at least 0.60. C5c, the share of false alarms, that is, of answers past their line on a right model, at 30 sets, is 0.188 ± 0.025 of 250 draws against a threshold of at most 0.12.

## 2. The versions proposed

| version | holds | closes when |
|---|---|---|
| 0.1.1 | small term, repairs, default of sets (section 3), C13, a criterion not yet written on the usual checks of `toy/benchmark.py` | C1 on the collisions at 0.60 under a scoring fixed first, C5c met or named not met, counts made again (section 5) |
| 0.2.0 | bounds from the term nobody names, an agent in the researcher's place validated by blind runs | C12, a criterion not yet written on how often a bound holds the truth, two public cases |
| 0.3.0 | correlated noise, or counted events | a criterion for each noise, C7 held on each |
| 0.4.0 | a loop in seconds, ten curves or more, a model as a formula | the clock project (`twinfit/EXAMPLES.md`) with ten curves in a time to fix (about 10 s on five on the author's machine) |
| 1.0 | files, commands and pages fixed, a paper | a paper submitted, three outside cases |

The author has not ruled on this order.

## 3. What 0.1.1 holds

### 3.1 The small term

The author has counted, toward 0.1.1 and with `toy/broad.py`, a program that ships with 0.1.1, a statistic for a broadening that is small against the noise. It is computed from `twinfit.linear.Linear` of a fit, whose `J` and `r` are the whitened derivatives of the mean and the whitened residuals with the own parameters of each record taken out. For a shared parameter $p$, column $j$ of `J`, and record $k$, with $a_k$ the rows of record $k$ in that column and $r_k$ its residuals:

$$w_k = a_k \cdot a_k, \qquad e_k = \frac{a_k \cdot r_k}{w_k},$$

the excess of $p$ that record $k$ asks on its own. Records with $w_k = 0$ are left out, and fewer than three leave $p$ unread. Along a numeric setting $s$ of the design, with $z_k = s_k - \bar{s}$ over the records kept:

$$S_{p,s} = \frac{\left|\sum_k z_k w_k e_k\right|}{\sqrt{\sum_k z_k^2 w_k}}, \qquad C_p = \sum_k w_k \left(e_k - \bar{e}\right)^2, \quad \bar{e} = \frac{\sum_k w_k e_k}{\sum_k w_k}.$$

$S_{p,s}$ is a standard normal and $C_p$ a chi square of one degree of freedom less than the records, at first order where $p$ is one value across them. The line of each is its 95th percentile over the twin's sets, 20 at least, as the report's own lines are. $S$ along the temperature and $C$ of the width each stood past its line on 100 of 100 draws of the collisions. Either it or the report stood past its line on 0.320 ± 0.030 of 250 right models, against C5c's 0.12. The author keeps the statistic if it reaches the 0.60 of C1 on the collisions, which it did. Its share, however, counts the draws past the line and not the place put first, so its scoring under C1 is still open.

### 3.2 Repairs from the blind red teams

- `along` (`twinfit/mechanisms.py`) takes no `lower` and no `upper`, so its factor 1 + a z can pass zero and a width falling with power was kept: add both, and have `best.md` say whether a parameter rises or falls.
- `best.md` and `stop_sentence` (`twinfit/loop.py`) name all four ends of the loop, but the stop rule of `twinfit/agent/AGENT.md` names two: add the noise first and no term worth keeping.
- `scan` warns of uneven gaps in its fine branch alone: warn in the coarse branch too.
- `python -m twinfit guide continuing` is refused, and `twinfit/CONTINUING.md` is reached through the folder that `guide readme` names: add it to the guides.
- `python -m twinfit.tutorial --check` says that the outputs differ and not which: name the first line that differs.

### 3.3 The default number of sets, 30 or 100

On the same draws `runs/pointing_own_*` (30 sets, each draw against its own twin) and `runs/pointing_own100_*` (100 sets, made by no program here, read by `docs/validation.py`) give a false alarm of 0.188 ± 0.025 at 30 sets and of 0.076 ± 0.017 at 100 sets, of 250 right models each. The second lies under the 0.09 that `toy/thresholds.toml` derives for 100 sets. The share of the collisions put first or tied on draws 1 to 50 is 0.520 ± 0.071 of 50 at both. Therefore 100 sets meet C5c and leave C1 unchanged, at 3.2 times the time: one right-model draw takes 66 s at 30 sets and 212 s at 100 (`toy/pointing.py`). The author has not ruled on the default.

The default stands in `twinfit/loop.py`, `twinfit/__main__.py` and, as the counts' own 30, `OWN_SETS` of `toy/pointing.py` and `SETS` of `toy/chosen.py` and `toy/gallery_count.py`. Only `OWN_SETS` moves C5c, and moving it makes 1060 units owed (`python toy/pointing.py units`).

## 4. The first step

The table printed by `python -m twinfit terms` lacks the entry `exponential_width`. The first step is a test that catches this, `tests/test_terms_lists.py`, in a new folder. It runs with plain Python, and with pytest where pytest is installed:

```python
import sys

sys.path.insert(0, ".")          # the repository's package, run from its root, and not an installed one

from twinfit.__main__ import NAMED, terms  # noqa: E402


def test_every_named_term_is_in_the_table():
    table = terms().split("\n\n")[1]
    for call, _ in NAMED:
        assert call.split("(")[0] + "(" in table, call


if __name__ == "__main__":
    test_every_named_term_is_in_the_table()
    print("1 passed")
```

Before the edit, `python tests/test_terms_lists.py` fails with `AssertionError: exponential_width()`. The README's line on `python -m twinfit terms`, which says that the table lacks `exponential_width`, goes with this edit. In the docstring of `twinfit/mechanisms.py`, between the `lorentzian_width` and `pedestal(stretch)` entries, with no blank line, insert:

```text
    exponential_width(along=None)
                             the same with the kernel exp(-|u| / b) / (2 b), b the decay length in the unit of the
                             abscissa. along scales b as for the two widths above
```

In `twinfit/version.py` set `__version__ = "0.1.1rc1"` and add after the last entry an entry `"0.1.1rc1": {"digest": "PENDING", "changes": ("one line",)}`.

```bash
python -c "import matplotlib.pyplot"   # builds the font cache once
python -m twinfit.tutorial --write
python -m twinfit.tutorial --check
python tests/test_terms_lists.py
python toy/pointing.py --check
python toy/gallery_count.py --check
```

`--write` also draws again the figures of `twinfit/examples/figures`, so a change of one line shows binary changes there, and takes one to four minutes. The `changes` of an entry is a tuple of strings, as in the entries before it. `--check` ends `the tutorial: every output of TUTORIAL.md is what its command prints now`. The test prints `1 passed`, and `pointing --check` exits 0. `gallery_count --check` exits 1, `made from other files than those on disk now`: its 220 draws, `python toy/count.py gallery_count`, are owed before release, about 12 minutes on nine cores of a laptop.

After `--write`, since the tutorial prints the version and is hashed, run the stamp `python -m twinfit.version --stamp LAST`, LAST being the number that `python -m twinfit --version` printed before the change. Python prints a `RuntimeWarning` about `twinfit.version` before the answer, which is harmless. It writes the digest into the current entry, so that `twinfit.version.check()` returns `0.1.1rc1`. It was not run for this page.

## 5. The rules of the work

1. **A change of the package is a new number**, `0.1.1rcN`. The digest hashes every `.py`, `.md`, `.json`, `.toml`, `.txt` and `.csv` file of `twinfit/` but `version.py`. `tools/met.py` reads a table made at another number than `release` in `toy/thresholds.toml` as meeting nothing, a `.postN` table counting as its release: count every table, then move `release`, then run `python docs/validation.py`.
2. **A count is made again after a change of the code it runs,** by `python toy/count.py PROGRAM` (`tools/autoqueue.py add` does nothing here). `gallery_count` hashes bytes, so even a comment in `twinfit/loop.py` stales it, the others ignoring docstrings and comments. `runs/` keeps no rows of `chosen` or `scoreboard`, and three units of `chosen` read the author's recorded traces, which are not public: C6 and C7 cannot be counted again on his traces. The loop of C6 and the bars of C7 can be run on the project that `python -m twinfit example rubidium DIR` writes, synthetic traces of the same model with their truth beside them, which is one draw and gives no share.

## 6. What is not on the roadmap

- The author's own rubidium line, on traces that are private. His aim with it is the measurement of four quantities: β_self, the self broadening coefficient, α, the differential polarizability of the two levels, the beam radius, and the efficiency of the retroreflection.
- A claim that no other tool does this. The only comparison here is C10, 0.917 ± 0.025 of 120.
