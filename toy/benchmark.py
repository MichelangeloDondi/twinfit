"""The benchmark beside the methods of the field, row M444 of the plan, his messages 262 to 270 and 352 (6a): on the
same draws as the count of `toy/power.py`, how often a missing term is seen by the package, by the usual checks of an
analyst, the runs test, the binned check and the check of neighbouring samples, which is what the statistic of Durbin
and Watson reads, by the chi square of the fit, and by their hybrid, the package with the usual checks at the level
that it leaves them. Each check is set at the level at which it reports on 5 draws of 100 of the right model of its
case, so that every method is read at one false alarm. And which named term the largest gain of a fit, a scan of the
Bayesian criterion over terms of one parameter each, puts first. Nothing is counted again: the rows are those of the
power count.

    python toy/benchmark.py table COPY        writes COPY/toy/benchmark.json from COPY/toy/power.json
    python toy/benchmark.py page COPY OUT     draws the page of the benchmark, one figure, into OUT (a .pdf or a .png)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import niceness  # noqa: E402  every process of this work at nice 10, his message 184

LEVEL = 0.05          # the false alarm at which every method is read, on the right model of its case
USUAL = ("runs", "bins", "neighbours")          # the usual checks that each row of the power count holds, as p-values
SAID = {"package": "the package", "runs": "the runs test", "bins": "the binned check", "neighbours": "neighbouring samples, as Durbin and Watson",
        "chi2": "the chi square of the fit", "hybrid": "the package and the usual checks"}
HALF = ", at 0.5 of the size of the twin"


def right_of(name: str) -> str:
    """The right model of a case whose term is hidden: "3 hidden under 2 with nothing free, at 0.5 of the size of the
    twin" is read against "2 with nothing free, right"."""
    return name.replace(HALF, "").split(" hidden under ", 1)[1] + ", right"


def share(hits: int, of: int) -> dict:
    p = hits / of if of else float("nan")
    return {"share": p, "standard_error": float(np.sqrt(p * (1.0 - p) / of)) if of else float("nan"), "of": of}


def good(rows: list) -> list:
    return [r for r in rows if "failed" not in r and r.get("usual")]


def levels(right: list) -> dict:
    """The level of each method at which it reports on `LEVEL` of the draws of a right model: the quantile of its own
    p-values there, and for the chi square the quantile of the chi square per degree of freedom from above."""
    out = {name: float(np.quantile([r["usual"][name] for r in right], LEVEL)) for name in USUAL}
    out["chi2"] = float(np.quantile([r["chi2"] / r["dof"] for r in right], 1.0 - LEVEL))
    alone = np.mean([bool(r["missing"]) for r in right])          # the package's own false alarm, which the hybrid adds to
    least = [min(r["usual"].values()) for r in right if not r["missing"]]
    room = max(0.0, LEVEL - alone) * len(right) / max(1, len(least))          # the share of the rest that the usual checks may add
    out["hybrid"] = float(np.quantile(least, room)) if least and room > 0.0 else 0.0
    return out


def seen(r: dict, lv: dict) -> dict:
    """Which method reports a missing term on one draw."""
    out = {"package": bool(r["missing"]), "chi2": r["chi2"] / r["dof"] > lv["chi2"]}
    out.update({name: r["usual"][name] < lv[name] for name in USUAL})
    out["hybrid"] = out["package"] or min(r["usual"].values()) < lv["hybrid"]
    return out


def table(power: dict, version: str) -> dict:
    """The benchmark from the table of the power count: for each case, each method's share of draws on which it reports,
    the right models read at their own false alarm, and the summary that `tools/met.py` reads for C10."""
    cases = power.get("unsuspected", {})
    out = {"package": version, "level": LEVEL, "methods": SAID, "cases": {}}
    for name, c in cases.items():
        if name.endswith(", right"):
            continue
        right = good(cases.get(right_of(name), {}).get("rows", []))
        rows = good(c.get("rows", []))
        if len(right) < 20 or not rows:          # a level read on fewer than 20 draws of the right model says nothing at 5 per cent
            continue
        lv = levels(right)
        on = [seen(r, lv) for r in rows]
        named = {}
        for r in rows:
            term = (r.get("term") or {}).get("name") or "no term is fitted"
            named[term] = named.get(term, 0) + 1
        out["cases"][name] = {"right": right_of(name), "levels": lv, **{m: share(sum(s[m] for s in on), len(on)) for m in SAID},
                              "false_alarm": {m: share(sum(seen(r, lv)[m] for r in right), len(right)) for m in SAID},
                              "named_first": {k: v / len(rows) for k, v in sorted(named.items(), key=lambda kv: -kv[1])}}
    whole = [c for n, c in out["cases"].items() if HALF not in n]
    pooled = lambda m: share(sum(round(c[m]["share"] * c[m]["of"]) for c in whole), sum(c[m]["of"] for c in whole))          # noqa: E731
    best = max(USUAL + ("chi2",), key=lambda m: pooled(m)["share"]) if whole else None
    out["summary"] = {"the package, pooled": pooled("package"), "the hybrid, pooled": pooled("hybrid"),
                      **({f"the best usual check, pooled, {SAID[best]}": pooled(best)} if best else {})} if whole else {}
    return out


def write(copy: Path) -> Path:
    copy = Path(copy)
    power = json.loads((copy / "toy" / "power.json").read_text())
    sys.path.insert(0, str(copy))
    from twinfit.version import __version__
    path = copy / "toy" / "benchmark.json"
    path.write_text(json.dumps(table(power, power.get("package") or __version__), indent=1) + "\n")
    return path


def pages_fill(text: str, width: int = 14) -> str:
    """A label cut into lines of about `width` letters at its spaces."""
    lines, now = [], ""
    for word in str(text).split():
        if now and len(now) + 1 + len(word) > width:
            lines, now = lines + [now], word
        else:
            now = f"{now} {word}".strip()
    return "\n".join(lines + [now])


def page(copy: Path, out: Path) -> Path:
    """One figure: for each term hidden at the size of the twin, the share of draws on which each method reports it, at
    one false alarm of 5 per cent on the right model, and the words of the page under 60."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    held = json.loads((Path(copy) / "toy" / "benchmark.json").read_text())
    cases = [(n, c) for n, c in held["cases"].items() if HALF not in n]
    methods = list(SAID)
    fig, ax = plt.subplots(figsize=(13, 7.5))
    width = 0.8 / len(methods)
    for k, m in enumerate(methods):
        ax.bar([i + k * width for i in range(len(cases))], [c[m]["share"] for _, c in cases], width, yerr=[c[m]["standard_error"] for _, c in cases], label=SAID[m], capsize=3)
    ax.set_xticks([i + 0.4 - width / 2 for i in range(len(cases))])
    names = {}
    try:          # the name of each detail of the toy, where the copy holds its values: the page of another program names the cases as they are
        import tomllib
        names = {k: v.get("name", k) for k, v in tomllib.loads((Path(copy) / "toy" / "config.toml").read_text()).get("detail", {}).items()}
    except (OSError, ValueError):
        pass
    ax.set_xticklabels([pages_fill(names.get(n.split(" hidden")[0], n.split(" hidden")[0])) for n, _ in cases], fontsize=15)
    ax.set_ylabel("share of draws that see it", fontsize=15)
    ax.set_ylim(0, 1.05)
    ax.tick_params(labelsize=15)
    ax.legend(fontsize=15, frameon=True, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    alarm = np.mean([c["false_alarm"]["package"]["share"] for _, c in cases]) if cases else float("nan")
    ax.set_title(f"A hidden term seen, {held['package']}", fontsize=17)
    fig.text(0.01, 0.015, f"The usual checks are set at 5 per cent of false alarm, the package stands at its own, {100 * alarm:.1f} per cent.", ha="left", fontsize=15)          # the package is not set at 5 per cent (review 99, moderate 2)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    fig.savefig(out)
    plt.close(fig)
    return Path(out)


def main(argv: list) -> int:
    niceness.to_ten()
    if argv[1:2] == ["table"] and len(argv) == 3:
        print(f"written: {write(Path(argv[2]))}")
        return 0
    if argv[1:2] == ["page"] and len(argv) == 4:
        print(f"written: {page(Path(argv[2]), Path(argv[3]))}")
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
