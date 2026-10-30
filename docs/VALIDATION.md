# What the package was measured to do

Every share below is read by `tools/met.py` from the tables of the counts in `toy/`, against the thresholds of `toy/thresholds.toml`, for the release 0.1.0. This page is written by `python docs/validation.py` from those tables, and `python docs/validation.py --check` says whether it still is the page they give. No number here is typed.

A count draws synthetic traces from the toy model of `toy/` with a known truth, hides a term or none, runs the package on each draw and reads what its report says. A share is the fraction of draws on which the report did what the criterion asks, with its standard error, the square root of p(1 - p)/n over its n draws. A threshold was drafted before the counts, and the header of `toy/thresholds.toml` and `docs/DECISIONS.md` say which were set again after tables were read. `toy/README.md` says how each count is made and how to make one draw of it again.

## Not met at this release

- C1, the place, collisions: 0.440 ± 0.050 of 100 draws, against at least 0.60
- C5c, a false alarm on a right model against its own twin at the default sets, pooled: 0.188 ± 0.025 of 250 draws, against at most 0.12

## Every criterion

| criterion | reading | share | draws | threshold | |
|---|---|---|---|---|---|
| C1 | the place, pooled over the hidden terms | 0.780 ± 0.021 | 400 | at least 0.60 | met |
| C1 | the place, collisions | 0.440 ± 0.050 | 100 | at least 0.60 | **not met** |
| C1 | the place, gas_width | 0.730 ± 0.044 | 100 | at least 0.60 | met |
| C1 | the place, laser_noise | 0.9900 ± 0.0099 | 100 | at least 0.60 | met |
| C1 | the place, transit | 0.960 ± 0.020 | 100 | at least 0.60 | met |
| C4 | the setting, pooled over the terms that change along one | 0.600 ± 0.035 | 200 | at least 0.60 | met |
| C5 | a false alarm on a right model, pooled | 0.083 ± 0.016 | 300 | at most 0.12 | met |
| C5c | a false alarm on a right model against its own twin at the default sets, pooled | 0.188 ± 0.025 | 250 | at most 0.12 | **not met** |
| C2 | the kind on the bench | 0.500 ± 0.079 | 40 | at least 0.50 | met |
| C2b | the kind on the held-out questions | 0.736 ± 0.040 | 120 | at least 0.30 | met |
| C5b | something missing on a right model, pooled | 0.040 ± 0.014 | 200 | at most 0.08 | met |
| C6 | reached_the_fifth | 1.000 ± 0 | 20 | at least 0.80 | met |
| C6 | refused_the_absent | 1.000 ± 0 | 20 | at least 0.85 | met |
| C7 | one bar holds the truth, collisions, the noise from neighbouring samples | 0.750 ± 0.097 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, gas_width, the noise from neighbouring samples | 0.65 ± 0.11 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, laser_noise, the noise from neighbouring samples | 0.70 ± 0.10 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, transit, the noise from neighbouring samples | 0.800 ± 0.089 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, collisions, the toy's law of the noise | 0.800 ± 0.089 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, gas_width, the toy's law of the noise | 0.65 ± 0.11 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, laser_noise, the toy's law of the noise | 0.65 ± 0.11 | 20 | 0.40 to 0.97 | met |
| C7 | one bar holds the truth, transit, the toy's law of the noise | 0.800 ± 0.089 | 20 | 0.40 to 0.97 | met |
| C10 | the package, pooled | 0.917 ± 0.025 | 120 | at least 0.76 | met |
| C11 | pooled | 0.759 ± 0.029 | 220 | at least 0.50 | met |

## What each criterion measures, and where it is counted

| criterion | what it measures | the count |
|---|---|---|
| C1 | the place of each term of tiers 2 to 5, hidden at the size of the twin, put first by the report | `toy/pointing.py`, its table `toy/pointing.json` |
| C4 | the setting along which the hidden term changes, named first where it changes along one | `toy/pointing.py`, its table `toy/pointing.json` |
| C5 | a right model read against the twin of its case: some answer beyond the twin, a false alarm | `toy/pointing.py`, its table `toy/pointing.json` |
| C5c | a right model read against the twin of its own fit at the loop's default of 30 sets, as a researcher runs it: some answer beyond the twin, the false alarm the default gives | `toy/pointing.py`, its table `toy/pointing.json` |
| C2 | what the blind report puts first on the bench of planted shapes at the toy's noise, the place, the shape, the noise and the setting pooled over its questions, against keys that name no term | `toy/scoreboard.py`, its table `toy/scoreboard.json` |
| C2b | the same on the held-out questions, at most 0.20 under the bench | `toy/scoreboard.py`, its table `toy/scoreboard.json` |
| C5b | a right model: the package says that something is missing, a false alarm of toy/power.py | `toy/power.py`, its table `toy/power.json` |
| C6 | the loop with the researcher's choices reaches the fifth model, and refuses the term that the toy lacks | `toy/chosen.py`, its table `toy/chosen.json` |
| C7 | one bar of each parameter of the final model holds its truth, on each arm of the loop count apart, the noise read from neighbouring samples and the toy's law of the noise: the collisions and the transit, held at theory, by the fit that frees them as a test, his message 340, row M480 | `toy/chosen.py`, its table `toy/chosen.json` |
| C10 | the package sees a hidden term more often than the best usual check of an analyst by two standard errors of their difference, each usual check set at a false alarm of 5 per cent on the right model of its case and the package at its own, pooled over the terms at the size of the twin, on the draws of toy/power.py | `toy/benchmark.py`, its table `toy/benchmark.json` |
| C11 | another field: place, shape of the residuals and setting pooled over the projects of the gallery, against keys that name no term | `toy/gallery_count.py`, its table `toy/gallery.json` |

## Shares that the README quotes beside the criteria

No criterion of `toy/thresholds.toml` reads these, so `tools/met.py` does not print them. Each is read here from its source.

| reading | share | draws | source |
|---|---|---|---|
| the place of a hidden term at half the size of the twin, pooled over the four | 0.128 ± 0.017 | 400 | `toy/pointing.json` |
| a false alarm on a right model against its own twin at 100 sets | 0.076 ± 0.017 | 250 | `runs/pointing_own100_*` |
| the place on the gallery of eleven projects of other fields | 0.588 ± 0.055 | 80 | `toy/gallery.json` |

What was not measured at all, an outside user, the agentic version and terms after the fifth, is said in `twinfit/CONTINUING.md`.

## Read it again

```bash
python tools/met.py M456,M505,M457,M458,M441,M444,M459 --counted
```

```bash
python toy/pointing.py --check
```

```bash
python toy/again.py "2 hidden under 1" 7
```

The first prints every reading of this page. The second compares the digest of the code of the count of the place and the setting with that of its table and sums again the rows that the table holds. The third makes one draw of that count again from its seed and sets it beside its kept row.
