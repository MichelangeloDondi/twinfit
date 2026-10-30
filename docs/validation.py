"""docs/VALIDATION.md written from the tables of the repository: every criterion of toy/thresholds.toml as tools/met.py
reads it, so that no number of the page is typed. Run again after a count is made again.

    python docs/validation.py           writes docs/VALIDATION.md from the repository's own tables
    python docs/validation.py --check   ends without error when docs/VALIDATION.md is the page that the tables give
"""
from __future__ import annotations

import math
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import met  # noqa: E402

PAGE = ROOT / "docs" / "VALIDATION.md"
ROWS = ("M456", "M505", "M457", "M458", "M441", "M444", "M459")          # the place and the setting, the false alarm at the default, the bench of shapes, the false alarm of the package, the loop and its bars, the package against canonical checks, the gallery
PROGRAMS = {"pointing": "toy/pointing.py", "power": "toy/power.py", "scoreboard": "toy/scoreboard.py", "chosen": "toy/chosen.py", "benchmark": "toy/benchmark.py", "gallery": "toy/gallery_count.py"}


def pm(share: float, se: float) -> str:
    if se <= 0:
        return f"{share:.3f} ± 0"
    places = max(0, 1 - math.floor(math.log10(se)))
    return f"{share:.{places}f} ± {se:.{places}f}"


def bound(r: dict) -> str:
    if r.get("at_least") is not None and r.get("at_most") is not None:
        return f"{r['at_least']:.2f} to {r['at_most']:.2f}"
    if r.get("at_least") is not None:
        return f"at least {r['at_least']:.2f}"
    if r.get("at_most") is not None:
        return f"at most {r['at_most']:.2f}"
    return "see its definition"


def quoted(root: Path) -> list:
    """The shares that the README quotes beside the criteria, which no criterion reads: the place pooled over the four
    hidden terms at half the size of the twin, from toy/pointing.json, the false alarm of the right models against the
    twin of each draw's own fit at 100 sets, from the kept rows runs/pointing_own100_*, and the place on the gallery of
    other fields, from toy/gallery.json."""
    import json
    cases = json.loads((root / "toy" / "pointing.json").read_text())["cases"]
    half = [v["place"] for k, v in cases.items() if " hidden under " in k and k.endswith(", at 0.5 of the size of the twin")]
    n, hits = sum(x["of"] for x in half), sum(x["share"] * x["of"] for x in half)
    rows = [json.loads(line) for f in sorted((root / "runs").glob("pointing_own100_*_right_*.jsonl")) for line in f.read_text().splitlines() if line.strip()]
    beyond = sum(1 for r in rows if r.get("beyond"))
    where = json.loads((root / "toy" / "gallery.json").read_text())["by_question"]["where"]
    se = lambda s, m: math.sqrt(s * (1 - s) / m)          # noqa: E731
    return [("the place of a hidden term at half the size of the twin, pooled over the four", pm(hits / n, se(hits / n, n)), n, "`toy/pointing.json`"),
            ("a false alarm on a right model against its own twin at 100 sets", pm(beyond / len(rows), se(beyond / len(rows), len(rows))), len(rows), "`runs/pointing_own100_*`"),
            ("the place on the gallery of eleven projects of other fields", pm(where["share"], se(where["share"], where["of"])), where["of"], "`toy/gallery.json`")]


def page(root: Path = ROOT) -> str:
    thresholds = tomllib.loads((root / "toy" / "thresholds.toml").read_text())
    release = thresholds["release"]
    rows, defs, unmet = [], [], []
    for row in ROWS:
        got = met.met(row, root)
        for name, c in got["criteria"].items():
            t = thresholds.get(name, {})
            program = PROGRAMS.get(t.get("table", ""), "")
            defs.append(f"| {name} | {c['what']} | `{program}`, its table `toy/{t.get('table', '')}.json` |")
            if c.get("why_not_read"):
                rows.append(f"| {name} | {c['why_not_read']} | | | not read |")
                continue
            for r in c["readings"]:
                share = pm(r["share"], r["standard_error"]) if isinstance(r.get("share"), (int, float)) and isinstance(r.get("standard_error"), (int, float)) else "not scored"
                mark = "met" if r["met"] else "**not met**"
                rows.append(f"| {name} | {r['reading']} | {share} | {r['of']} | {bound(r)} | {mark} |")
                if not r["met"]:
                    unmet.append(f"{name}, {r['reading']}: {share} of {r['of']} draws, against {bound(r)}")
    return f"""# What the package was measured to do

Every share below is read by `tools/met.py` from the tables of the counts in `toy/`, against the thresholds of `toy/thresholds.toml`, for the release {release}. This page is written by `python docs/validation.py` from those tables, and `python docs/validation.py --check` says whether it still is the page they give. No number here is typed.

A count draws synthetic traces from the toy model of `toy/` with a known truth, hides a term or none, runs the package on each draw and reads what its report says. A share is the fraction of draws on which the report did what the criterion asks, with its standard error, the square root of p(1 - p)/n over its n draws. A threshold was drafted before the counts, and the header of `toy/thresholds.toml` and `docs/DECISIONS.md` say which were set again after tables were read. `toy/README.md` says how each count is made and how to make one draw of it again.

## Not met at this release

{chr(10).join("- " + u for u in unmet) or "- none"}

## Every criterion

| criterion | reading | share | draws | threshold | |
|---|---|---|---|---|---|
{chr(10).join(rows)}

## What each criterion measures, and where it is counted

| criterion | what it measures | the count |
|---|---|---|
{chr(10).join(dict.fromkeys(defs))}

## Shares that the README quotes beside the criteria

No criterion of `toy/thresholds.toml` reads these, so `tools/met.py` does not print them. Each is read here from its source.

| reading | share | draws | source |
|---|---|---|---|
{chr(10).join(f"| {a} | {b} | {c} | {d} |" for a, b, c, d in quoted(root))}

What was not measured at all, an outside user, the agentic version and terms after the fifth, is said in `twinfit/CONTINUING.md`.

## Read it again

```bash
python tools/met.py {",".join(ROWS)} --counted
```

```bash
python toy/pointing.py --check
```

```bash
python toy/again.py "2 hidden under 1" 7
```

The first prints every reading of this page. The second compares the digest of the code of the count of the place and the setting with that of its table and sums again the rows that the table holds. The third makes one draw of that count again from its seed and sets it beside its kept row.
"""


if __name__ == "__main__":
    text = page()
    if sys.argv[1:] == ["--check"]:
        same = PAGE.exists() and PAGE.read_text() == text
        print("docs/VALIDATION.md is the page that the tables give" if same else "docs/VALIDATION.md is not the page that the tables give: python docs/validation.py")
        sys.exit(0 if same else 1)
    PAGE.write_text(text)
    print(f"written {PAGE.relative_to(ROOT)}")
