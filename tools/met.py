"""Whether a row of section B of the plan is met: its tables read against the thresholds of toy/thresholds.toml, which
were written before the counts (the red team of his message 331, severe 1: the proofs read the label of a table, and
passed whatever its rates said).

Each criterion of the thresholds names its row and its table. The table is read in the counting copy that
runs/queue_at names, or in the copy given, from its toy/ folder. A table that is missing, made at another version than
the release, or counted in the other mode, meets nothing, and so does a criterion whose rates the table does not hold
yet. The mode of a criterion is "default", the package with no dictionary, or "flagged", the dictionary on the flag of
the user (his message 337), and a table that names no mode was counted with the dictionary on. Each share is said with
its standard error and its draws.

    python3 tools/met.py M456               every criterion of row M456, and exit 0 only where all are met
    python3 tools/met.py M456 --counted     the same readings, and exit 0 where every criterion is read at the release, met
                                            or not: the proof of a row of a count, whose verdict against the drafts is his
                                            to judge (the slap of 2026-10-23 04:27:23: a row proven by its verdict alone waits
                                            on the package, not on its count)
    python3 tools/met.py M456 --copy NAME   the tables of ~/twinfit-copies/NAME, or of the folder NAME where it is a path
                                            (`.` the repository): with no copy and no runs/queue_at, the repository's own
    python3 tools/met.py M456,M457 --counted   several rows, each read, and exit 0 only where every one is
"""
from __future__ import annotations

import json
import math
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import niceness  # noqa: E402  every process of this work at nice 10, his message 184
COPIES = Path.home() / "twinfit-copies"


def share(hits: float, of: int) -> dict:
    """A share with its standard error, sqrt(p(1-p)/n)."""
    p = hits / of if of else math.nan
    return {"share": p, "standard_error": math.sqrt(p * (1.0 - p) / of) if of else math.nan, "of": of}


def pooled(parts: list) -> dict:
    """Shares of several cases pooled by their draws."""
    of = sum(p["of"] for p in parts)
    return share(sum(p["share"] * p["of"] for p in parts), of)


def hidden_cases(table: dict) -> list:
    """The cases of toy/pointing.json with a term hidden at the size of the twin, not at half of it."""
    return [c for name, c in table.get("cases", {}).items() if c.get("term") and "at 0.5 of the size" not in name]


READ = ("C1", "C2", "C2b", "C3", "C4", "C5", "C5b", "C5c", "C6", "C7", "C8", "C10", "C11")          # the criteria that `readings` reads: one that it does not met nothing, as C5b of row M458 until 2026-10-24
ASKED = ("where", "shape", "what", "changes along")          # the questions of the blind report that a key of toy/scoreboard.py answers. Since the report is blind, "what" is the law of the noise and the shape a question of its own, so C2 read on "what" alone scored one question of noise


def asked_of(cell: dict | None, draws: int = 0) -> list | None:
    """The answers of a group of toy/scoreboard.py pooled over the questions of the blind report: [right, asked, draws],
    the draws those of the group's questions, `group_draws`, since one draw answers several (reviews 97 and 98)."""
    if not cell:
        return None
    return [sum(cell.get(q, [0, 0])[0] for q in ASKED), sum(cell.get(q, [0, 0])[1] for q in ASKED), draws]


def group_draws(table: dict, group: str, weight: str = "law") -> int:
    """The draws of the questions of one group of toy/scoreboard.py at one weight, read from its reports: the most that
    one question was asked on read 25 draws of a bench of 40 (review 98, moderate 1)."""
    names = {q.get("name") for q in table.get("questions", []) if q.get("group") == group}
    return sum(len(r.get("draws") or []) for key, r in (table.get("reports") or {}).items() if key.rsplit("|", 1)[0] in names and key.rsplit("|", 1)[-1] == weight)


def on_draws(asked: list) -> dict:
    """A share of answers whose standard error and count are those of its draws: answers of one draw are not independent."""
    p = asked[0] / asked[1]
    return {"share": p, "standard_error": math.sqrt(p * (1.0 - p) / asked[2]) if asked[2] else math.nan, "of": asked[2]}


def readings(key: str, crit: dict, table: dict) -> dict:
    """The shares that a criterion reads in its table, by name, or a reason why there are none."""
    if key == "C1":
        parts = [c["place"] for c in hidden_cases(table) if c.get("place")]
        return {"the place, pooled over the hidden terms": pooled(parts)} | {f"the place, {c['term']}": c["place"] for c in hidden_cases(table) if c.get("place")} if parts else {}
    if key == "C4":
        parts = [c.get("setting") or c.get("along") for c in hidden_cases(table) if (c.get("setting") or c.get("along")) and ((c.get("expected") or {}).get("setting") or (c.get("expected") or {}).get("along"))]          # `setting` in the blind report from 0.1.0rc16, `along` before: C4 read none of the first table at 0.1.0rc16
        return {"the setting, pooled over the terms that change along one": pooled(parts)} if parts else {}
    if key == "C5":
        parts = [c["something_beyond_the_twin"] for c in table.get("cases", {}).values() if c.get("term") is None and c.get("something_beyond_the_twin")]
        return {"a false alarm on a right model, pooled": pooled(parts)} if parts else {}
    if key == "C5c":          # row M505, toward 0.1.1: the right models each read against the twin of its own fit at the loop's default sets, the false alarm a researcher sees (the own_twin part of toy/pointing.json, row M104)
        parts = [c["something_beyond_the_twin"] for c in (table.get("own_twin") or {}).values() if c.get("term") is None and c.get("something_beyond_the_twin")]
        return {"a false alarm on a right model against its own twin at the default sets, pooled": pooled(parts)} if parts else {}
    if key == "C5b":          # the right models of toy/power.py, on which the package says that something is missing: C5b had no reading, and row M458 could meet nothing
        parts = [c["missing"] for name, c in (table.get("unsuspected") or {}).items() if name.endswith(", right") and isinstance(c.get("missing"), dict)]
        return {"something missing on a right model, pooled": pooled(parts)} if parts else {}
    if key in ("C2", "C2b"):
        cells = (table.get("counted_draws") or {}).get("law") or {}
        bench, held = asked_of(cells.get("bench"), group_draws(table, "bench")), asked_of(cells.get("held out"), group_draws(table, "held out"))
        out = {}
        if bench and bench[1]:
            out["the kind on the bench"] = on_draws(bench)
        if key == "C2b" and held and held[1]:
            out["the kind on the held-out questions"] = on_draws(held)
        return out
    if key == "C3":
        s = table.get("summary") or {}
        return {k: s[k] for k in ("first_class", "three_classes") if isinstance(s.get(k), dict)}
    if key == "C6":
        s, n = table.get("summary") or {}, (table.get("summary") or {}).get("draws") or 0
        return {k: share(s[k] * n, n) for k in ("reached_the_fifth", "refused_the_absent") if isinstance(s.get(k), (int, float)) and n}
    if key == "C7":          # each arm of the table, named: the bar of the collisions held the truth on 0.20 of 20 draws with the noise read from neighbouring samples and on 0.82 of 100 with the toy's law, and the page read the first alone (review 102, moderate 1)
        out = {}
        for s, arm in (((table.get("summary") or {}).get("within_one_bar") or {}, "the noise from neighbouring samples"), ((((table.get("law") or {}).get("summary")) or {}).get("within_one_bar") or {}, "the toy's law of the noise")):
            out.update({f"one bar holds the truth, {p}, {arm}": s[p] if isinstance(s.get(p), dict) else {"share": None, "of": 0} for p in crit.get("parameters", [])})          # a parameter that the table does not score is not met: two of four were (the red team of his message 331, severe 4). An arm not counted reads as not counted, and met nothing alone (review 103, minor 8)
        return out
    if key == "C8":
        s = table.get("summary") or {}
        n = s.get("draws") or 0
        stopped = (s.get("stopped") or {}).get("admitted")
        return {"ended on an admitted model": share(stopped * n, n)} if n and isinstance(stopped, (int, float)) else {}
    if key in ("C10", "C11"):
        s = table.get("summary") or {}
        return {k: v for k, v in s.items() if isinstance(v, dict) and "share" in v}
    return {}


def judged(key: str, crit: dict, got: dict) -> list:
    """For each reading: whether it meets its threshold. The second share of C3 and of C6 has a threshold of its own."""
    out = []
    for name, s in got.items():
        p = s.get("share")
        if key == "C6" and name == "refused_the_absent":
            lo, hi = crit.get("refused_at_least"), None
        elif key == "C3" and name == "three_classes":
            lo, hi = crit.get("three_at_least"), None
        elif key == "C2b" and name == "the kind on the held-out questions":
            bench = got.get("the kind on the bench", {}).get("share")
            lo, hi = (bench - crit["within"] if isinstance(bench, (int, float)) else None), None
        elif key == "C2b":
            continue
        elif key == "C10" and name == "the package, pooled":          # the package against the best usual check of an analyst, at one false alarm, his message 352 (6a)
            best = next((s for k, s in got.items() if k.startswith("the best usual check, pooled")), None)
            gap = 2.0 * math.hypot(s.get("standard_error") or 0.0, (best or {}).get("standard_error") or 0.0)          # two standard errors of the difference, review 98, moderate 4
            lo, hi = (best["share"] + gap if isinstance((best or {}).get("share"), (int, float)) else None), None
        elif key == "C10":
            continue
        elif "around" in crit:          # a band of `errors` standard errors of the draws of the reading itself: set at 20 draws, it read 5.5 of them wide at 100 and could not fail (review 103, moderate 1)
            half = crit["errors"] * math.sqrt(crit["around"] * (1.0 - crit["around"]) / s["of"]) if s.get("of") else None
            lo, hi = (crit["around"] - half, crit["around"] + half) if half is not None else (None, None)
        else:
            lo, hi = crit.get("at_least"), crit.get("at_most")
        ok = isinstance(p, (int, float)) and not math.isnan(p) and (lo is not None or hi is not None) and (lo is None or p >= lo) and (hi is None or p <= hi)          # a criterion with no bound is met by nothing: a held-out kind of 1 in 100 read met beside a bench of no reading (review 92, moderate 9)
        se = s.get("standard_error")
        if se is None and ok is not None and isinstance(p, (int, float)) and s.get("of"):          # a table that gives a share and its draws alone
            se = share(p * s["of"], s["of"])["standard_error"]
        out.append({"criterion": key, "reading": name, "share": p, "standard_error": se, "of": s.get("of"), "at_least": lo, "at_most": hi, "met": ok})
    return out


def met(row: str, copy: Path, thresholds: Path = ROOT / "toy" / "thresholds.toml") -> dict:
    """Every criterion of `row`: its readings and whether each is met, or why it cannot be read. `met` is True only where
    every criterion of the row is read and met."""
    t = tomllib.loads(Path(thresholds).read_text())
    mine = {k: v for k, v in t.items() if isinstance(v, dict) and v.get("row") == row}
    if not mine:
        raise SystemExit(f"REFUSED: no criterion of toy/thresholds.toml names row {row}")
    unread = sorted(k for k in mine if k not in READ)
    if unread:
        raise SystemExit(f"REFUSED: tools/met.py has no reading of {', '.join(unread)}: a criterion that nothing reads meets nothing, whatever its count")
    out = {"row": row, "release": t["release"], "copy": str(copy), "criteria": {}, "met": True, "min_of": {k: int(v.get("min_of", 1)) for k, v in mine.items()}}
    for key, crit in mine.items():
        path = Path(copy) / "toy" / f"{crit['table']}.json"
        if not path.is_file():
            why = f"toy/{crit['table']}.json is not written in {copy}"
        else:
            table = json.loads(path.read_text())
            other = base(table.get("package")) != base(t["release"])          # a table of another version meets nothing, a post release of the same code reading the tables of its release
            mode, asked = mode_of(table), crit.get("mode", "default")
            got = {} if other or mode != asked else readings(key, crit, table)
            why = (f"toy/{crit['table']}.json was made at {table.get('package') or 'no version that it names'}, not at {t['release']}" if other
                   else f"toy/{crit['table']}.json was counted in the {mode} mode, not in the {asked} mode" if mode != asked
                   else None if got else f"toy/{crit['table']}.json holds no reading of {key} yet")
        lines = [] if why else judged(key, crit, got)
        out["criteria"][key] = {"what": crit["what"], "why_not_read": why, "readings": lines}
        out["met"] = out["met"] and not why and bool(lines) and all(x["met"] for x in lines)
    out["counted"] = all(not c["why_not_read"] and c["readings"] and all((x["of"] or 0) >= mine[k].get("min_of", 1) for x in c["readings"])
                         for k, c in out["criteria"].items())          # every reading a share, on the draws that the count asks less a tenth: a parameter not scored, or 98 failed draws of 100, read as counted (review 92, moderate 1)
    return out


BLIND_FROM = 16          # the release candidate of 0.1.0 from which the package is blind by default, the dictionary on its flag alone


def base(version) -> str:
    """The release whose code a version is: 0.1.0.post1 and 0.1.0.post1.dev2 are 0.1.0, since a post release changes its
    documents alone (his messages 441 and 442, the bus test), and a table counted under either reads as one of 0.1.0."""
    return re.sub(r"(?:\.post\d+)?(?:\.dev\d+)?$", "", str(version or ""))


def mode_of(table: dict) -> str:
    """The mode in which a table was counted: the one it names, or, where it names none, "default" from 0.1.0rc16 on,
    whose counts pass no flag, and "flagged" before, since every count before his message 337 ran with the dictionary on
    (the red team of that message, minor 4). No count table of 0.1.0rc16 names its mode, so every one read as flagged."""
    v = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)(?:rc(\d+))?(?:\.post\d+)?(?:\.dev\d+)?", str(table.get("package") or ""))
    blind = bool(v) and (tuple(int(x) for x in v.groups()[:3]) > (0, 1, 0) or (tuple(int(x) for x in v.groups()[:3]) == (0, 1, 0) and (v.group(4) is None or int(v.group(4)) >= BLIND_FROM)))
    return table.get("mode") or ("default" if blind else "flagged")          # 0.1.0 itself and every later version blind (review 97, minor 3)


def mine_of(result: dict, key: str) -> int:
    """The draws that a reading of criterion `key` needs, as the result keeps it."""
    return int(result.get("min_of", {}).get(key, 1))


def said(result: dict) -> str:
    lines = [f"row {result['row']} against the thresholds of the release {result['release']}, the tables of {result['copy']}:"]
    for key, c in result["criteria"].items():
        if c["why_not_read"]:
            lines.append(f"  {key}, {c['what']}: not read, {c['why_not_read']}")
            continue
        for x in c["readings"]:
            if isinstance(x["of"], int) and x["of"] < mine_of(result, key):
                lines.append(f"  {key}, {x['reading']}: read on {x['of']} draws, where the count asks {mine_of(result, key)} at least")
            limit = " and ".join(s for s in ((f"at least {x['at_least']:.2f}" if x["at_least"] is not None else ""), (f"at most {x['at_most']:.2f}" if x["at_most"] is not None else "")) if s)
            value = pm(x["share"], x["standard_error"]) + f" on {x['of']}" if isinstance(x["share"], (int, float)) and isinstance(x["standard_error"], (int, float)) else "not scored by the table"
            lines.append(f"  {key}, {x['reading']}: {value}, {limit}: {'met' if x['met'] else 'NOT MET'}")
    lines.append(("met" if result["met"] else "not met") + (", every criterion read at the release" if result.get("counted") else ", some criterion not read at the release"))
    return "\n".join(lines)


DICTIONARY = ('"a wider pedestal under the curve"', '"a second component of the curve"', '"a polynomial of the level"', '"a tail on one side"', '"a periodic function"')          # the kinds of the report of the dictionary as a key holds them, in quotation marks: his words quoted in a docstring are no key
KEYED = {          # each program that scores the answers of the report: the keys of the dictionary's report that it may no longer hold, what of the blind report it must hold, and whether its count is tried on two draws first (row M472)
    "toy/scoreboard.py": (DICTIONARY + ('"along":',), ('"changes along"', '"shape"'), True),
    "toy/gallery_count.py": (DICTIONARY + ('"along":',), ('"changes along"', '"shape"'), True),
    "toy/pointing.py": (DICTIONARY + ('"along":', '("covariate", "along")'), ('"changes along"', '"shape"'), True),
    "toy/power.py": (DICTIONARY, ('c.family == "trend"',), True),
    "toy/noiseless.py": (DICTIONARY + ('"what", "along")',), ('"changes along"',), False),          # the noiseless run of the folder, which no count makes
}
TRIED = "runs/tried_keys.json"          # what `python toy/tried.py` of the copy keeps: each program's rows on two draws, with the digest of its code


def keyed(copy: Path) -> list:
    """What stands between the counts of `copy` and the blind report, row M472: a copy or a program that is not there,
    a key of the dictionary's report that a program still holds, a question of the blind report that it does not hold,
    and a count not tried on two draws as its code now stands, or whose try failed. A count keyed on the dictionary's
    report scores zero by construction on the blind report, about 45 core hours for nothing (reviews/REDTEAM_343.md,
    S1, and review 95). Empty where nothing stands."""
    import hashlib
    copy = Path(copy)
    if not copy.is_dir():
        return [f"{copy}: no such copy"]
    tried = json.loads((copy / TRIED).read_text()).get("programs", {}) if (copy / TRIED).exists() else {}
    out = []
    for name, (old, blind, counted) in KEYED.items():
        path = copy / name
        if not path.is_file():
            out.append(f"{name}: no such program")
            continue
        text = path.read_text()
        out += [f"{name}: {key} in {text.count(key)} places" for key in old if key in text]
        out += [f"{name}: holds no {key}" for key in blind if key not in text]
        t = tried.get(name)
        if not counted:
            continue
        if t is None:
            out.append(f"{name}: not tried on two draws, `python toy/tried.py` in the copy")
        elif t.get("code") != hashlib.sha256(path.read_bytes()).hexdigest():
            out.append(f"{name}: its code moved since its try")
        elif t.get("failed") or len(t.get("rows", [])) < 2:
            out.append(f"{name}: its try failed on {t.get('failed', 2)} of 2 draws" + (f", {t['error']}" if t.get("error") else ""))
    return out


def pm(share: float, se: float) -> str:
    """A share and its standard error, the error at two significant digits and the share at its decimals (rule 9): 0.440
    +- 0.050 where it printed 0.440 +- 0.05, and the count at none or all of the draws, where the formula gives a bar of
    zero that is not the uncertainty."""
    if se <= 0:
        return f"{share:.3f} +- 0"
    places = max(0, 1 - math.floor(math.log10(se)))
    return f"{share:.{places}f} +- {se:.{places}f}"


def main(argv: list) -> int:
    niceness.to_ten()          # this command and what it starts at nice 10 (his message 184)
    if len(argv) == 3 and argv[1] == "--keys":          # the proof of row M472
        found = keyed(COPIES / argv[2])
        print("\n".join(found) if found else f"the counts of {argv[2]} key on the blind report alone")
        return 1 if found else 0
    words = [w for w in argv[1:] if w != "--counted"]
    if len(words) not in (1, 3) or (len(words) == 3 and words[1] != "--copy"):
        print(__doc__)
        return 2
    named = ROOT / "runs" / "queue_at"
    if len(words) == 3:          # a copy of the author's machine by its name, or any folder by its path, as a reader of the public repository has (his messages 441 and 442)
        copy = Path(words[2]).expanduser() if ("/" in words[2] or words[2] == ".") else COPIES / words[2]
    else:
        copy = Path(named.read_text().strip()) if named.exists() else ROOT          # the repository's own tables where no queue names a copy
    if not (copy / "toy").is_dir():
        print(f"REFUSED: {copy} holds no toy/ folder of tables")
        return 2
    results = [met(row, copy) for row in words[0].split(",")]          # several rows at once, as the deliverable M474 reads its five counts
    print("\n".join(said(result) for result in results))
    return 0 if all(result["counted"] if "--counted" in argv else result["met"] for result in results) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
