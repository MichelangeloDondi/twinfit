"""The marks of a report of the bus test, docs/BUS_TEST.md, read so that none is typed: the mark out of 10 of each
of the six steps, their total out of 60, the stumbles by rating, and whether the run passes by the rule of 60 of 60, no stumble that blocks or misleads, and a verdict that says yes. Two runs in a row must pass,
and two runs in a row whose stumbles do not fall end the row with what stands.

    python tools/bus_score.py REPORT [REPORT ...]     one line of JSON for each report, then whether the rule is met
    python tools/bus_score.py --self-test
"""
import json
import re
import sys
from pathlib import Path

RATINGS = ("BLOCKS", "MISLEADS", "SLOWS", "POLISH")


def score(text: str) -> dict:
    steps = {}
    parts = re.split(r"(?m)^##\s+(?:Step\s+)?(\d)\b", text)
    for k, body in zip(parts[1::2], parts[2::2]):
        marks = re.findall(r"Mark:\s*(\d+)\s*/\s*10", body)
        if marks:
            steps[int(k)] = int(marks[-1])
    counts = {r: len(re.findall(rf"(?m)^\s*-\s+\*\*{r}\*\*", text)) for r in RATINGS}
    tail = text.strip().splitlines()[-12:]
    verdict = next((line for line in reversed(tail) if re.search(r"\bverdict\b|can continue|cannot continue|can not continue", line, re.I)), "")
    yes = bool(re.search(r"\bcan continue\b", verdict, re.I)) and not re.search(r"cannot|can not|not continue", verdict, re.I)
    total = sum(steps.values())
    passed = sorted(steps) == [1, 2, 3, 4, 5, 6] and total == 60 and counts["BLOCKS"] == 0 and counts["MISLEADS"] == 0 and yes
    return {"marks": steps, "total": total, "stumbles": counts, "verdict_yes": yes, "passed": passed}


def rule(scores: list) -> str:
    """Whether the row ends: two passes in a row, or two runs in a row whose stumbles do not fall."""
    if len(scores) >= 2 and scores[-1]["passed"] and scores[-2]["passed"]:
        return "two runs in a row pass: the bus test is passed"
    weight = lambda s: sum(s["stumbles"].values())          # noqa: E731
    if len(scores) >= 3 and weight(scores[-1]) >= weight(scores[-2]) >= weight(scores[-3]):
        return "two runs in a row did not reduce the stumbles: the row ends with what stands"
    return "go on: repair what the last run found and run again"


def self_test() -> int:
    bad = 0
    full = "\n".join(f"## Step {k}. x\n\ndone\n\nMark: 10/10\n" for k in range(1, 7)) + "\n- **POLISH** a word\n\nVerdict: a researcher can continue his work from what PyPI serves, 60 of 60.\n"
    s = score(full)
    for what, ok in (("a full report passes", s["passed"] and s["total"] == 60),
                     ("a mislead fails it", not score(full.replace("**POLISH**", "**MISLEADS**"))["passed"]),
                     ("a step at 9 fails it", not score(full.replace("Mark: 10/10", "Mark: 9/10", 1))["passed"]),
                     ("a verdict that cannot fails it", not score(full.replace("can continue", "cannot continue"))["passed"]),
                     ("two passes end the row", rule([s, s]).startswith("two runs in a row pass")),
                     ("stumbles that do not fall end it", rule([dict(s, passed=False)] * 3).startswith("two runs in a row did not reduce"))):
        print(("ok    " if ok else "FAIL  ") + what)
        bad += not ok
    return bad


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        sys.exit(1 if self_test() else 0)
    got = [score(Path(p).read_text()) for p in sys.argv[1:]]
    for p, s in zip(sys.argv[1:], got):
        print(json.dumps({"report": p} | s))
    print(rule(got))
