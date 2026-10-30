"""The page of row M474, his messages 327, 331 and 343: the five tiers at the release, one figure, each criterion of
`toy/thresholds.toml` as `tools/met.py` reads it on the tables of the counting copy, its share with two standard errors
against the band of its threshold, and each criterion not read yet named with the reason.

    python toy/tiers.py COPY OUT          the page into OUT, a .pdf or a .png, from the tables of ~/twinfit-copies/COPY
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import met  # noqa: E402
import niceness  # noqa: E402  every process of this work at nice 10, his message 184

ROWS = ("M456", "M457", "M458", "M459", "M441", "M444", "M505")          # the rows of M474, the benchmark beside them, and M505, the own twin at the default (review 112, minor 6)
SHORT = {"C1": "place first", "C4": "setting first", "C5": "false alarm, the twin", "C2": "shape and noise, bench", "C2b": "the same, held out",
         "C5b": "false alarm, the package", "C5c": "false alarm, the twin of its own fit at the default", "C6": "the loop", "C7": "one bar holds its truth", "C11": "another field", "C10": "the package against the best usual check plus two errors"}


LITERATURE = (          # the names of the tables in the notation of the literature, and no word that the criterion above says (review 102, minor 4)
              ("collisions", "β_self"), ("gas_width", "β of the permeated gas"), ("laser_noise", "σ_L"), ("transit", "Γ_tr"),
              ("reached_the_fifth", "reached the fifth detail"), ("refused_the_absent", "refused the absent term"),
              ("one bar holds the truth, ", ""), ("the noise from neighbouring samples", "noise of neighbouring samples"), ("the toy's law of the noise", "the toy's law of the noise"),
              (" over the hidden terms", ""), (" over the terms that change along one", ""),
              ("a false alarm on a right model against its own twin at the default sets, ", ""))          # said by the label of C5c already, and whole it pushed the title off the page (review 112, minor 6)


def label(key: str, reading: str, of=None) -> str:
    """The label of a reading: the criterion on one line and the reading on the next, whole, in the notation of the
    literature and with its draws. Cut at 58 letters, four labels lost their ends (reviews 99 to 101), and with no draws
    a share of 20 and one of 400 read alike (review 102, minor 4)."""
    for name, said in LITERATURE:
        reading = reading.replace(name, said)
    return f"{key} {SHORT.get(key, key)}\n{reading}" + (f", {of} draws" if isinstance(of, int) else "")


def lines(copy: Path, rows: tuple = ROWS) -> list:
    """One line a reading: the row, the criterion, the reading, its share and standard error, its bounds, met or not, or
    the reason why it is not read."""
    out = []
    for row in rows:
        got = met.met(row, copy)
        for key, c in got["criteria"].items():
            if c["why_not_read"]:
                out.append({"row": row, "key": key, "label": f"{key} {SHORT.get(key, key)}", "why": c["why_not_read"]})
                continue
            for x in c["readings"]:
                if not isinstance(x["share"], (int, float)):          # a parameter that the table does not score: said, not drawn (review 99, moderate 1)
                    out.append({"row": row, "key": key, "label": label(key, x["reading"]), "why": "not scored by the table"})
                    continue
                out.append({"row": row, "key": key, "label": label(key, x["reading"], x["of"]), "share": x["share"], "se": x["standard_error"],
                            "lo": x["at_least"], "hi": x["at_most"], "met": x["met"], "of": x["of"], "short": isinstance(x["of"], int) and x["of"] < got["min_of"].get(key, 1)})
    return out


def grouped(rows: list) -> list:
    """The labels of the page, the criterion said once at the head of its rows: said on each, the page held 266 words of
    labels where his message 149 asks less text (review 103, minor 5)."""
    return [r["label"] if k == 0 or rows[k - 1]["key"] != r["key"] else r["label"].split("\n", 1)[-1] for k, r in enumerate(rows)]


def mark(r: dict) -> str:
    """The colour of a reading on the page: green where it meets its threshold, red where not, grey where not read."""
    if "why" in r:
        return "0.45"
    return "#1a7f37" if r["met"] else "#b42318"


def page(copy: Path, out: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = lines(copy)
    fig, ax = plt.subplots(figsize=(15, 0.8 * len(rows) + 2.2))          # two lines a label
    for k, r in enumerate(rows):
        y = len(rows) - 1 - k
        if "why" in r:
            ax.text(0.01, y, "not scored" if r["why"] == "not scored by the table" else "not read yet", va="center", fontsize=15, color=mark(r))
        else:
            if r["lo"] is not None or r["hi"] is not None:          # a band where the threshold has a bound, none over the whole axis
                lo, hi = (r["lo"] if r["lo"] is not None else 0.0), (r["hi"] if r["hi"] is not None else 1.0)
                ax.barh(y, hi - lo, left=lo, height=0.6, color="#cfe8cf", zorder=0)
            ax.errorbar(r["share"], y, xerr=2.0 * (r["se"] or 0.0), fmt="o" if not r["short"] else "s", color=mark(r), ms=9, capsize=4, zorder=2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(list(reversed(grouped(rows))), fontsize=15)
    ax.set_xlim(-0.08, 1.03)          # a dot at 1.0 was half cut by the frame (review 102, minor 4), and a bar under 0 too (review 103, minor 5)
    ax.tick_params(axis="x", labelsize=15)
    ax.set_xlabel("share of draws, two standard errors, its threshold in green\na square: fewer draws than its threshold asks", fontsize=15)          # the square had no legend (review 102, minor 4)
    ax.set_title(f"The five tiers at {met.tomllib.loads((ROOT / 'toy' / 'thresholds.toml').read_text())['release']}: green met, red not met", fontsize=17)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return Path(out)


def main(argv: list) -> int:
    niceness.to_ten()
    if len(argv) == 3:
        print(f"written: {page(met.COPIES / argv[1] if not Path(argv[1]).is_dir() else Path(argv[1]), Path(argv[2]))}")
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
