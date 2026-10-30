# The toy model, the bench of twinfit

This folder holds the bench on which [`twinfit`](../twinfit/README.md) was measured: a model of the author's two photon line of rubidium from 5S to 6S and the programs that count what the package does on traces drawn from it. [`../twinfit/CONTINUING.md`](../twinfit/CONTINUING.md) is the successor's first page.

## 1. The toy and its words

The toy has five *details* (`const.toml` [details]): detail $k$ holds the first $k$ terms of section 2. A *trace* is the signal in V against the frequency $\nu$ in MHz at one power $P$ in mW and one temperature $T_\mathrm{C}$ in °C, a *draw* the 8 synthetic traces of one number. A *right model* is fitted to traces of its own detail. A *hidden term* is one the fitted model lacks: the model of detail $k-1$, keeping the free parameters of detail $k$ but the new term (none for $k=3$: `3 hidden under 2 with nothing free`), fitted to traces of detail $k$. [`TUTORIAL.md`](../twinfit/TUTORIAL.md) section 1 defines the package's words, except its *place*, where the report points: the core, the flanks, the wings or the asymmetry. A *set* is a synthetic copy of the traces from a twin, a *false alarm* some answer beyond the twin on a right model, a *share* the fraction of draws on which something happens. *The twin of a case* is made once from the case's first draw, *the twin of a draw's own fit* for each draw. A *unit* is a draw or a part of a twin, a *condition* a named group of units, a *row* one JSON line of [`../runs/`](../runs/) in a file named with the *digest*, a hash of the code and values it depends on. The author's *queue*, not public, ran units in *waves*.

**Environment.** Every command ran on an arm64 Mac, not on Linux, with Python 3.14.6, numpy 2.5.2 and scipy 1.18.0. The wheel on PyPI holds `twinfit/` alone, so the bench is obtained by cloning https://github.com/MichelangeloDondi/twinfit and running from its root. The digest of a count hashes the syntax tree of its code, which differs between versions of Python, so under Python 3.12.13 every check but that of `gallery_count.py` reports that its table was made from other files.

## 2. The line and its five terms

Widths are in MHz: $W$ and $\Gamma_\mathrm{nat}$ are Lorentzian full widths at half maximum, $\sigma_L$ is the standard deviation of a Gaussian, $\Gamma_\mathrm{tr}$ a rate. A trace is built on the Fourier transform $L(\nu)$ of an even function $C(t)$, $t$ in µs, of three factors.

$$C(t)=e^{-\pi W(T)\lvert t\rvert}\,e^{-2\pi^2\sigma_L^2t^2}\,\frac{1}{1+(2\pi b(T)t)^2},\qquad L(\nu)=\int C(t)\,e^{2\pi i\nu t}\,dt$$

$$W(T)=\Gamma_\mathrm{nat}+\beta_\mathrm{self}\,\frac{n(T)}{n_0}+\beta,\qquad b(T)=\Gamma_\mathrm{tr}\sqrt{T/T_0}$$

Here $T=T_\mathrm{C}+273.15$ is in kelvin. The density is $n(T)=p(T)/k_BT$ in cm$^{-3}$, $p(T)=101325\text{ Pa}\times10^{A-B/T}$, $A=4.312$, $B=4040$ K (`const.toml` [density]), with $n_0=10^{15}$ cm$^{-3}$ and $T_0=403.15$ K. The area of a trace goes as $P^2n(T)$, set so that the trace at 225 mW and 130 °C has the height 4.24 V (`const.toml` [reference_height]).

|detail|term, symbol, `config.toml` name|value, unit|source|
|---|---|---|---|
|1|natural width, $\Gamma_\mathrm{nat}$, `natural_width`|3.493 ± 0.013 MHz|$1/(2\pi\tau)$, $\tau=45.57\pm0.17$ ns, Gomez 2005 (`const.toml`)|
|2|collisions between rubidium atoms, $\beta_\mathrm{self}$, `collisions`|3.50 ± 0.75 MHz per $10^{15}$ atoms per cm³|129 ± 13 kHz per mTorr on 7S, Zameroski and others 2014, carried to 6S (`config.toml`)|
|3|permeated gas, $\beta$, `gas_width`|0.037 ± 0.015 MHz|joint fit of the 8 recorded traces (`config.toml`)|
|4|laser noise, $\sigma_L$, `laser_noise`|0.369 ± 0.028 MHz|the same fit|
|5|transit, $\Gamma_\mathrm{tr}$, `transit`|1.040 ± 0.024 MHz|$u/(2\pi w_0)$ at $T_0$, $u$ the most probable speed of $^{85}$Rb (`const.toml` [atom]), $w_0=43.0\pm1.0$ µm adopted by the author in place of a fit of the recorded traces (`const.toml` [beam]), the bar from the waist's alone|

The main fit frees the natural width in details 1 and 2, the gas width from 3 and the laser noise from 4, and holds collisions and the transit at theory, freeing them only in a fit beside it, as a test, on all 8 traces (`config.toml` [detail.N]). The uncertainty of $\beta_\mathrm{self}$ is the error of the carry from 7S to 6S alone: the carry is also off by a third on average, in one direction, which no uncertainty holds (`config.toml`). The gas width, constant to about 5 per cent in [`DERIVATION.md`](DERIVATION.md) section 3, is exactly constant here and absorbs every constant width the model lacks (`DERIVATION.md` section 6(b)). `DERIVATION.md` derives details 1 to 3 alone, and its section 7 lists what that model leaves out, among it the toy's transit and laser noise.

## 3. The eight traces

The power sweep is 25, 75, 125, 175 and 225 mW at 130 °C, the temperature sweep 70, 90, 110 and 130 °C at 225 mW, sharing the trace at 225 mW and 130 °C: 8 traces (`const.toml` [conditions]) of one hyperfine peak, 993.4192 nm in `config.toml` [recorded], each of 1177 points from $-25$ to $25$ MHz (`const.toml` [axis]). Three heights of the trace at 225 mW and 130 °C appear in the files: 4.24 V, the model's reference (`const.toml` [reference_height]), 4.48 V, the first repetition fitted alone, and 4.277 V, fitted jointly under the fifth model (`twin.toml` [heights]), which `twin.py` draws with.

## 4. The twin and how a trace is drawn

`traces.py` draws the model's own traces, which `pointing.py`, `power.py`, `draws.py` and the synthetic runs of `chosen.py` take through `fits.synthetic`: the signal of section 2 plus Gaussian noise of variance $a^2+bh+ch^2+ds^2$, with $h$ the signal in V and $s$ its slope in V per MHz. `scoreboard.py` adds noise of its own to its noiseless signal (seed $[2026,10,19]$ and the draw number), its field group and `gallery_count.py` draw the package's examples, and `benchmark.py` draws nothing. The 32 laws of $a$ (V), $b$ (V) and $c$, one for each peak and each of the 8 traces, and $d=0.00135\pm0.00026$ MHz², one for all, are in `../data/noise/noise_laws.json`. The seed of `traces.py` is $[2026,10,12]$ (`const.toml` [seed]), then the detail, the power, the temperature and, for a draw other than 0, its number.

`twin.py` draws synthetic traces that follow the recorded traces in height and noise (`twin.toml`). This is not the twin of the package, and no count uses it. `twin.toml` also holds a laser noise for each trace, 0.169 to 0.627 MHz, seven of them fitted alone and the one at 225 mW and 70 °C the joint value, and `twin.py` draws each trace with its own, while the model and its fits hold the one value of the joint fit, 0.369 ± 0.028 MHz. The author has since chosen the one value for all eight traces, which 0.1.1 carries. [`../twinfit/examples/rubidium_twin`](../twinfit/examples/rubidium_twin) holds its draw 0 at detail 5, up to the rounding of the files.

## 5. The counts

`pointing.py` has 13 cases: for each of details 2 to 5 the term hidden at the size of the twin and at half of it, and a right model for each of details 1 to 5. `EXPECTED` says where: the core for collisions, gas and laser noise, the flanks for the transit, and a change along $T$ for collisions and transit.

|program|what it counts, criteria of `thresholds.toml`|draws and sets|table, kept rows|cost|
|---|---|---|---|---|
|`pointing.py`|where the report points for each hidden term, and the false alarm of each right model: C1, C4, C5, C5c|30 hidden draws, then 100 where a share lies within two standard errors of 0.60: the table holds 100 on each hidden case. Right 60 against the twin of the case, 50 against its own at 30 sets. Twin of a case 100 sets|`pointing.json`. `../runs/pointing_*`|668 waves, 893636.3 s in `pointing.json`|
|`power.py`|12 of those cases, no right model of detail 5: how often the package says something is missing: C5b|hidden 30, right 50, cut window 30 (width read on the whole window and parts of it). 200 sets|`power.json`. `../runs/power_*`, without the cut window rows|35 waves, 13199.5 s in `power.json`|
|`benchmark.py`|the package, the usual checks, the chi square and the package with them, on the draws of `power.py`, all but the package at a false alarm of 5 per cent: C10|those of `power.py`|`benchmark.json`|not recorded|
|`scoreboard.py`|the blind report (no term named) on 36 questions in five groups (bench 8, held out 6, toy 13, field 8, waist 1): C2, C2b|5 draws a question in bench and toy, 20 in the rest. 100 sets|`scoreboard.json`. Readings not shipped|1602 units, not recorded|
|`chosen.py`|the loop offered the toy's four terms in order, then a pedestal that the toy lacks and it should refuse: C6, C7|46 units in 8 conditions: 20 draws with the noise read from neighbouring samples, 20 with the toy's law, and 6 single runs on one draw or on the recorded traces. 30 sets|`chosen.json`. Units not shipped|11 waves, 43 of the 46 units, 21465.2 s in `chosen.json`|
|`gallery_count.py`|another field, the 11 projects of [`../twinfit/examples`](../twinfit/examples): C11|20 draws a project, 30 sets|`gallery.json`. `../runs/gallery_count_*`|not recorded|

The cost is the sum of the seconds of the queue's waves, not the elapsed time. The cost of `pointing.json` also covers rows at sizes and numbers of sets that no table here counts. The author drafted the thresholds, and some were set again after tables were read (`../docs/DECISIONS.md`). Of the 250 right draws read against their own twin, 47 stand beyond it at 30 sets (`pointing.json`) and 19 at 100 sets (`../runs/pointing_own100_*`, counted by the author toward 0.1.1 with programs that are not in this release).

## 6. Checking and reproducing

```bash
python toy/pointing.py --check
python toy/power.py --check
python toy/chosen.py --check
python toy/scoreboard.py --check
python toy/gallery_count.py --check
python toy/draws.py --check
python docs/validation.py --check
```

Each prints one line: five contain `from the files as they are`, that of `draws.py` says `every model passes`, the last names `docs/VALIDATION.md`. Each of the five compares the digest of its code with that of its table and sums again the rows that the table holds, so a pass says that the table matches the code and reads nothing of `runs/`, whose kept rows `again.py` and `count.py` read.

```bash
python tools/met.py M456,M505,M457,M458,M441,M444,M459 --counted
```

It prints the 23 readings of [`../docs/VALIDATION.md`](../docs/VALIDATION.md). Each label is a row of the author's plan and names the criteria whose `row` field in `thresholds.toml` holds it, so one label reads one group of criteria. Two end NOT MET:

```text
  C1, the place, collisions: 0.440 +- 0.050 on 100, at least 0.60: NOT MET
  C5c, a false alarm on a right model against its own twin at the default sets, pooled: 0.188 +- 0.025 on 250, at most 0.12: NOT MET
```

`again.py` makes one draw of `pointing.py` again from its seed, and lists the 13 cases given no argument.

```bash
python toy/again.py "2 hidden under 1" 7
```

It ends `draw 7 of '2 hidden under 1' made again in 4 s: the kept row, field for field`. The seconds here and below are those of the author's machine.

`count.py` makes a count without the queue, `--draws K` keeping units 1 to $K$ of each condition of draws, and `--only TEXT` the conditions whose name holds TEXT. In a copy, move the rows to remake aside:

```bash
cp -R . ../twinfit-try
cd ../twinfit-try
mv runs/gallery_count_pedestal_*.jsonl ../pedestal_kept.jsonl
python toy/count.py gallery_count --workers 4 --only pedestal --draws 3
grep -c -F -x -f ../pedestal_kept.jsonl runs/gallery_count_pedestal_*.jsonl
```

The last line begins `3 units made in 13 s. The count is not whole, so toy/gallery.json is left as it is`, and `grep` prints `3`: the rows made are among the 20 kept. A whole recount by `count.py` equals the committed table draw by draw, not byte for byte, since it writes the rows in the order they end. In the copy again, since these write tables:

```bash
python toy/count.py power --workers 8 --only cut_window
python toy/pointing.py table
python toy/benchmark.py table .
```

The first makes the 30 cut window draws that `../runs/` lacks (46 s, 79 s with the table) and writes `power.json`, which lacks the section `cut_window` if written before. Each table equals the committed one but for `queue` (`tools/autoqueue.py` stands in for the queue, which counted no waves), and the column `before` of `pointing.json`.

## 7. Adding a sixth term

1. `toy/terms/NAME.py`: `width(T_celsius, values, const)` (a width in MHz added to $W$), `factor(t, T_celsius, values, const)` (a factor of $C(t)$) or both, reading `values["NAME"]`.
2. `const.toml`: `NAME = 6` under [details].
3. `config.toml`: `details = [..., 6]` under [folder], `[detail.6]` (`name`, `free`, `freed`) and `[parameter.NAME]` (`enters = 6`, `true`, `uncertainty`, `start`, `lower`, `upper`, `unit`, `what`, and `kind = "theory"` for a held value).
4. `pointing.py`: `EXPECTED[NAME]` and `WRITTEN[NAME]`, as for the four terms.
5. `chosen.py`: `NAME` in `ORDER` and `WHY`, `HELD` for a value of theory, and the 5 of `fits.synthetic(5, …)` changed to 6: it, `scoreboard.py` and `noiseless.py` draw at detail 5.

A missing `[detail.6]` raises `KeyError` at the first fit, a missing `EXPECTED` entry only when `pointing.py` writes its table, after its draws. Every count but that of `gallery_count.py` must then be made again, since the digest moves: their checks say the tables were made from other files, and the check of `draws.py` says model 6 has none, which `draws.py` makes when given 6 as its argument. Make the counts with `count.py`: `pointing.py table` alone writes a table of no case, which its `--check` passes.

## 8. What is not shipped

- `real.py` reads the author's recorded traces, which are private: the conditions `his`, `his law` and `his report` of `chosen.py` need them, so its table cannot be written.
- Not here: `adopt.py`, `noise.py`, `build.py` and `pages.py` (they fit the recorded traces and build the author's pages), `noiseless.json`, `rehearsal.json` and `places.json`, the package's own report on noiseless traces, from which `EXPECTED` of `pointing.py` was read. `noiseless.py` refuses to run outside the author's queue.
- The table of `scoreboard.py` needs no recorded trace, but 0 of its 72 reports are kept: `count.py` makes them when given `scoreboard`.
