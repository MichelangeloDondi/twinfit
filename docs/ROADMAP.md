# What comes after 0.1.0.post1

This page builds on [`CONTINUING.md`](../twinfit/CONTINUING.md), [`VALIDATION.md`](VALIDATION.md) and [`DECISIONS.md`](DECISIONS.md). Run from the repository root, on Python 3.14.6. A **draw** is one synthetic data set of the toy of `toy/`, its twin the **sets** (30 by default) drawn from the fitted model, a **line** the largest excess over all answers together that a set of a right model, one lacking nothing, reaches on 95 of 100 sets (`twinfit/TUTORIAL.md`).

## 1. Where 0.1.0.post1 stands

The code is that of 0.1.0, with new pages. `python -m twinfit --version` prints the number of this release, and `twinfit.version.check()` returns it: its digest, which hashes the pages too, is stamped in `twinfit/version.py`.

Not met (`docs/VALIDATION.md`): C1 on the collisions put first, 0.440 ± 0.050 of 100 draws (at least 0.60), and C5c, a false alarm (an answer past its line on a right model) at 30 sets, 0.188 ± 0.025 of 250 (at most 0.12).

## 2. The versions proposed, undecided

| version | holds | closes when |
|---|---|---|
| 0.1.1 | small term, repairs, default of sets (section 3), C13, a criterion not yet written on the usual checks of `toy/benchmark.py` | C1 on the collisions at 0.60 under a scoring fixed first, C5c met or named not met, counts made again (section 5) |
| 0.2.0 | bounds from the term nobody names, an agent in the researcher's place validated by blind runs | C12, a criterion not yet written on how often a bound holds the truth, two public cases |
| 0.3.0 | correlated noise, or counted events | a criterion for each noise, C7 held on each |
| 0.4.0 | a loop in seconds, ten curves or more, a model as a formula | the clock project (`twinfit/EXAMPLES.md`) with ten curves in a time to fix (about 10 s on five on the author's machine) |
| 1.0 | files, commands and pages fixed, a paper | a paper submitted, three outside cases |

## 3. What 0.1.1 holds

### 3.1 The small term

Counted by the author toward 0.1.1, by a program that is not in this release: a statistic for a broadening small against the noise. It is computed from `twinfit.linear.Linear` of a fit, whose `J` and `r` are the whitened derivatives of the mean and the whitened residuals with the own parameters of each record taken out. For a shared parameter $p$, column $j$ of `J`, and record $k$, with $a_k$ the rows of record $k$ in that column and $r_k$ its residuals:

$$w_k = a_k \cdot a_k, \qquad e_k = \frac{a_k \cdot r_k}{w_k},$$

the excess of $p$ that record $k$ asks on its own. Records with $w_k = 0$ are left out, and fewer than three leave $p$ unread. Along a numeric setting $s$ of the design, with $z_k = s_k - \bar{s}$ over the records kept:

$$S_{p,s} = \frac{\left|\sum_k z_k w_k e_k\right|}{\sqrt{\sum_k z_k^2 w_k}}, \qquad C_p = \sum_k w_k \left(e_k - \bar{e}\right)^2, \quad \bar{e} = \frac{\sum_k w_k e_k}{\sum_k w_k}.$$

$S_{p,s}$ is a standard normal and $C_p$ a chi square of one degree of freedom less than the records, at first order where $p$ is one value across them. The line of each is its 95th percentile over the twin's sets, 20 at least, as the report's own lines are. $S$ along the temperature and $C$ of the width each stood past its line on 100 of 100 draws of the collisions. Either it or the report stood past its line on 0.320 ± 0.030 of 250 right models, against C5c's 0.12. The author keeps it if it reaches C1's 0.60 on the collisions, which it did, but its share counts draws past the line, not the place put first: its scoring is open.

### 3.2 Repairs from the blind red teams

- `along` (`twinfit/mechanisms.py`) takes no `lower` and no `upper`, so its factor 1 + a z can pass zero and a width falling with power was kept: add both, and have `best.md` say whether a parameter rises or falls.
- `best.md` and `stop_sentence` (`twinfit/loop.py`) name all four ends of the loop, but the stop rule of `twinfit/agent/AGENT.md` names two: add the noise first and no term worth keeping.
- `scan` warns of uneven gaps in its fine branch alone: warn in the coarse branch too.
- `python -m twinfit guide continuing` is refused, and `twinfit/CONTINUING.md` is reached through the folder that `guide readme` names: add it to the guides.
- `python -m twinfit.tutorial --check` says that the outputs differ and not which: name the first line that differs.

### 3.3 The default number of sets, 30 or 100, undecided

On the same draws `runs/pointing_own_*` (30 sets, each draw against its own twin) and `runs/pointing_own100_*` (100 sets, made by no program here, read by `docs/validation.py`) give a false alarm of 0.188 ± 0.025 against 0.076 ± 0.017 of 250 right models, the second under the 0.09 that `toy/thresholds.toml` derives for 100 sets. The collisions, put first or tied on draws 1 to 50, do not move: 0.520 ± 0.071 of 50 at both. So 100 sets meet C5c and leave C1 as it is, at 3.2 times the time (`toy/pointing.py`, one right-model draw: 66 s against 212 s).

The default stands in `twinfit/loop.py`, `twinfit/__main__.py` and, as the counts' own 30, `OWN_SETS` of `toy/pointing.py` and `SETS` of `toy/chosen.py` and `toy/gallery_count.py`. Only `OWN_SETS` moves C5c, and moving it makes 1060 units owed (`python toy/pointing.py units`).

## 4. The first step

`python -m twinfit terms` lacks `exponential_width` in its table. Write `tests/test_terms_lists.py`, in a new folder. It runs with plain Python, and with pytest where it is installed:

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

Before the edit, `python tests/test_terms_lists.py` fails with `AssertionError: exponential_width()`. In the docstring of `twinfit/mechanisms.py`, between the `lorentzian_width` and `pedestal(stretch)` entries, with no blank line, insert:

```text
    exponential_width(along=None)
                             the same with the kernel exp(-|u| / a) / (2 a), a the decay length in the unit of the
                             abscissa. along scales a as for the two widths above
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

## 5. The rule of the work, in detail

1. **A change of the package is a new number**, `0.1.1rcN`. The digest hashes every `.py`, `.md`, `.json`, `.toml`, `.txt` and `.csv` file of `twinfit/` but `version.py`. `tools/met.py` reads a table made at another number than `release` in `toy/thresholds.toml` as meeting nothing, a `.postN` table counting as its release: count every table, then move `release`, then run `python docs/validation.py`.
2. **A count is made again after a change of the code it runs,** by `python toy/count.py PROGRAM` (`tools/autoqueue.py add` does nothing here). `gallery_count` hashes bytes, so even a comment in `twinfit/loop.py` stales it, the others ignoring docstrings and comments. `runs/` keeps no rows of `chosen` or `scoreboard`, and three units of `chosen` read the author's recorded traces, which are not public: C6 and C7 cannot be counted again here.

## 6. What is not on the roadmap

- The author's own rubidium line (β_self, α, the beam radius, the retroreflection), on traces that are private.
- A claim that no other tool does this: only C10 compares here, 0.917 ± 0.025 of 120.
