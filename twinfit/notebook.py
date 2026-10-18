"""The record of an analysis: its steps in order, its predictions stamped before the runs that they predict, its refusals.

An analysis made with an agent leaves a chat that nobody will read again. What must stay is a folder that a person can
audit with no chat: which rung was opened and with whose words, what was predicted before each run with noise, what
each step returned, and what was refused. This module keeps that folder and holds the order of the work:

    book = Notebook("analysis/")
    book.open_rung("simple profile", by="the researcher, 2026-10-09: start with the simple profile")
    book.step("identify", "simple profile", table)                        # a step and its table, written with the clock
    book.stamp("twin of the simple profile", "the spread of g is 0.9 to 1.1 of its bound ...", rung="simple profile")
    book.step("twin", "simple profile", {"admission": table, "readings": {...}}, stamp="twin of the simple profile")
    book.grade("twin of the simple profile", holds=True, what_came="1.03")
    book.admit("the fit", "simple profile", 1.0)                          # read from the twin step on record

What it refuses, each with a test that tries a careless way around it. Every refusal is written into the record.
  - a step of a rung that nobody opened, a step before the one that it needs, a step with no table or with a table that
    lacks what the step must hold (HOLDS)
  - a step that draws noise (DRAWS) with no prediction stamped for its rung before it, a second stamp under one name, a
    grade with no stamp, a grade before anything was run under the stamp, a second grade
  - a reading of the records that the twin of the rung did not admit at that level of noise: `admit` takes no table from
    its caller. It reads the last twin step that the record of the rung holds at that level, and a table typed by hand
    has no way in. Its line says how many verdicts that level had and how many of them refused, so a level
    that was run again until it passed shows in the record
  - a term added to the working model that was not first read on at least MIN_SETS synthetic sets that hold it, under
    its exact name

What it cannot check, and no page may say that it does: that the numbers of a table are true, since a twin step is
recorded with the table that its caller gives, that a prediction says anything, that it was written before the
computation and not only before its record, and whose words open a rung. It holds the order of its own calls, their
clock and the shape of what it is given. The rest is the word of whoever runs it, which the folder lets a person audit.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

FMT = "%Y-%m-%d %H:%M:%S"
ORDER = ("intake", "identify", "twin", "fit", "read", "propose", "candidate on the twin", "report")
NEEDS = {"identify": "intake", "twin": "identify", "fit": "twin", "read": "fit", "propose": "read", "candidate on the twin": "propose", "report": "propose"}
HOLDS = {"intake": ("records", "samples"),
         "identify": ("parameters", "weakest"),
         "twin": ("admission",),
         "fit": ("shared", "standard_deviations"),
         "read": ("chi_square", "degrees_of_freedom"),
         "propose": ("candidates", "leading_set"),
         "candidate on the twin": ("the_planted_candidate", "read", "share_of_sets_in_which_something_is_seen"),
         "report": ("what_it_does_not_show",)}
DRAWS = ("twin", "candidate on the twin")          # the steps that draw noise: each names the prediction stamped before it
MIN_SETS = 30
FIT = "the fit"                                    # the reading whose admission is the table `admission` of a twin step


class Refused(RuntimeError):
    """A step out of its order or with no table, a rung that nobody opened, a grade with no stamp, a reading not admitted."""


def _count(value) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


class Notebook:
    def __init__(self, folder):
        self.folder = Path(folder)
        (self.folder / "stamps").mkdir(parents=True, exist_ok=True)
        self.log = self.folder / "notebook.jsonl"

    # ---- the record
    def entries(self) -> list:
        if not self.log.exists():
            return []
        return [json.loads(line) for line in self.log.read_text().splitlines() if line.strip()]

    def _write(self, entry: dict) -> dict:
        entry = {"at": time.strftime(FMT), **entry}
        with self.log.open("a") as fh:
            fh.write(json.dumps(entry) + "\n")
        return entry

    def _refuse(self, what: str):
        """A refusal is part of the record: it is written, then raised."""
        self._write({"kind": "refusal", "what": what})
        raise Refused("REFUSED: " + what)

    def refusals(self) -> list:
        return [e["what"] for e in self.entries() if e["kind"] == "refusal"]

    def table_of(self, entry: dict) -> dict:
        return json.loads((self.folder / entry["file"]).read_text())

    # ---- rungs
    def open_rung(self, rung: str, by: str) -> dict:
        """A rung of the ladder is opened with the words of the person who decides. The notebook records the words: it
        cannot check whose they are. Fewer than three words are refused."""
        if len(str(by).split()) < 3:
            self._refuse(f"the rung '{rung}' is opened with the words of the person who decides it, and '{by}' is not a sentence of theirs")
        return self._write({"kind": "rung", "rung": rung, "opened_by": by})

    def open_by_default(self, rung: str, after: str) -> dict:
        """A rung that no person opened. The rung `after` ended, the default rule of the package chose a term, and the
        analysis goes on with it: the researcher can stop or go back at any report. The notebook records the rung as
        opened by default, with the term and what the choice rested on, so that a reader tells the rungs that a person
        opened from those that the rule did.

        Everything is read from the table of the report of `after` on record, under `default_choice`, as
        `criteria.choose` writes it: the term under `chosen`, and what the choice rests on under `rests_on`. Nothing is
        taken from the caller: the first two forms of this method took the term, then its ground and its level, from
        whoever called it (reviews 45 and 46). Refused: a rung `after` with no report, a rung that is already open, a
        report at which the researcher made a choice of their own, a report whose default choice is not a table, a report
        that already opened a rung by default (review 47), a
        report that chose no term or does not say what its choice rests on, a term that was not read on synthetic
        records that hold it, a ground of fewer than three words, and a guess whose level in the family is above 0.05
        or rests on fewer than 19 synthetic sets. What it cannot check is that the table of the report is true: it
        holds what was recorded, and a rung takes one report."""
        if after not in self.rungs() or not self.done("report", after):
            self._refuse(f"the rung '{rung}' is opened by default after '{after}', whose report is not on record: a default goes on from a rung that was closed")
        if rung in self.rungs():
            self._refuse(f"the rung '{rung}' is already open. A rung is opened once, by a person or by the rule")
        if any(e["kind"] == "choice" and e["rung"] == after for e in self.entries()):
            self._refuse(f"the researcher made a choice of their own at the report of '{after}'. The rung that follows it is theirs to open, in their words")
        reports = [e for e in self.entries() if e["kind"] == "step" and e["step"] == "report" and e["rung"] == after]
        choice = self.table_of(reports[-1]).get("default_choice")
        if choice is None:
            choice = {}
        if not isinstance(choice, dict):
            self._refuse(f"the rung '{rung}' is opened by default, and the report of '{after}' on record holds no table under its default choice: "
                         "a default goes on with the term that the report chose and with what that choice rests on")
        used = [e["rung"] for e in self.entries() if e["kind"] == "rung" and e.get("opened_by") == "default" and e.get("after") == after]
        if used:
            self._refuse(f"the report of '{after}' already opened the rung '{used[0]}' by default. A report chooses one term, and opens one rung")
        term, rests = choice.get("chosen"), choice.get("rests_on")
        if not term:
            self._refuse(f"the rung '{rung}' is opened by default, and the report of '{after}' on record chose no term. A default goes on with the term that the report of the rung before chose")
        self.add_term(term, after)
        if not isinstance(rests, dict) or not ("ground" in rests or {"level_in_the_family", "sets"} <= set(rests)):
            self._refuse(f"the report of '{after}' chose '{term}' and does not say what the choice rests on: a ground of the theory, or a level in the family with its number of sets")
        if "ground" in rests:
            if len(str(rests["ground"]).split()) < 3:
                self._refuse(f"the term '{term}' is said to stand on theory, and '{rests['ground']}' is not a sentence of that theory: a term that is added whatever a criterion says rests on its ground alone")
            kept = {"ground": str(rests["ground"])}
        else:
            level, sets = rests["level_in_the_family"], rests["sets"]
            if not (isinstance(level, (int, float)) and 0.0 < float(level) <= 0.05 and _count(sets) >= 19):
                self._refuse(f"the term '{term}' is a guess, and a guess is added by default only at a level in the family of 0.05 or under on 19 synthetic sets or more: "
                             f"the report holds a level of {level} on {sets} sets")
            kept = {"level_in_the_family": float(level), "sets": _count(sets)}
        return self._write({"kind": "rung", "rung": rung, "opened_by": "default", "after": after, "term": term, "rests_on": kept})

    def choose(self, rung: str, term: str, how: str, by: str) -> dict:
        """The researcher's own choice at the report of a rung, which takes the place of the default choice of that report: they can stop or go back at
        any report. `how` is "fitted" or "held", and `by` is their sentence. The record keeps the default and the choice, and `chosen(rung)` gives the
        one in force. Refused: a rung whose report is not on record, a way that is neither fitted nor held, a sentence of fewer than six words, and a
        term that was not read on synthetic records of that rung in the way that is chosen. Fitted: a run that planted the term at a size and fitted
        it. Held: a run that drew the size of the term from its theory and fitted the model with the term held, whose table says `held_at`."""
        if not self.done("report", rung):
            self._refuse(f"a choice is made at the report of a rung, and the rung '{rung}' has no report on record")
        if how not in ("fitted", "held"):
            self._refuse(f"a term is taken fitted or held, and '{how}' is neither")
        if len(str(by).split()) < 6:
            self._refuse(f"a choice of the researcher is recorded with their reason, in six words or more: '{by}' is too short to be one")
        read = [self.table_of(e) for e in self.entries() if e["kind"] == "step" and e["step"] == "candidate on the twin" and e["rung"] == rung]
        read = [r for r in read if r.get("the_planted_candidate") == term and (("held_at" in r) == (how == "held"))]
        if not read:
            self._refuse(f"the term '{term}' was not read {how} on synthetic records of the rung '{rung}'. A term enters the twin, in the way in which it "
                         "will enter the model, before it enters the model")
        return self._write({"kind": "choice", "rung": rung, "term": term, "how": how, "by": str(by).strip()})

    def chosen(self, rung: str):
        """The choice in force after the report of a rung: the researcher's last one if they made any, or else the default of the report. It is
        {"term", "how", "by"} with `by` their sentence or "default". None when the rung has no report, or when its report chose no term."""
        mine = [e for e in self.entries() if e["kind"] == "choice" and e["rung"] == rung]
        if mine:
            return {"term": mine[-1]["term"], "how": mine[-1]["how"], "by": mine[-1]["by"]}
        reports = [e for e in self.entries() if e["kind"] == "step" and e["step"] == "report" and e["rung"] == rung]
        choice = self.table_of(reports[-1]).get("default_choice") if reports else None
        if not isinstance(choice, dict) or not choice.get("chosen"):
            return None
        return {"term": choice["chosen"], "how": choice.get("how"), "by": "default"}

    def note(self, name: str, text: str) -> dict:
        """A dated note beneath a prediction, on the record: what changed after the stamp and before its result, such as a
        run that is no longer made. Refused where no such prediction was stamped, once it is graded, once the run under
        it is on record (review 47), and for a note of fewer than six words. The notebook knows a run to be under a prediction
        by the stamp that its step names: a run recorded with no stamp is not seen here (review 48). The stamp itself is never written again."""
        path = self._stamp_path(name)
        if not path.exists():
            self._refuse(f"the prediction '{name}' was never stamped: there is nothing to write a note beneath")
        if any(e["kind"] == "grade" and e["name"] == name for e in self.entries()):
            self._refuse(f"the prediction '{name}' is graded. A note beneath a prediction is written before its result")
        if any(e["kind"] == "step" and e.get("stamp") == name for e in self.entries()):
            self._refuse(f"the run under the prediction '{name}' is on record. A note beneath a prediction is written before its result")
        if len(str(text).split()) < 6:
            self._refuse(f"a note beneath the prediction '{name}' says what changed and why, in six words or more")
        now = time.strftime(FMT)
        with path.open("a") as fh:
            fh.write(f"\n(Added {now}: {str(text).strip()})\n")
        return self._write({"kind": "note", "name": name, "text": str(text).strip()})

    def opened_by_default(self) -> list:
        """The rungs that the default rule opened, in order: what a person reads first when they come back to an analysis that went on without them."""
        return [e["rung"] for e in self.entries() if e["kind"] == "rung" and e.get("opened_by") == "default"]

    def rungs(self) -> list:
        return [e["rung"] for e in self.entries() if e["kind"] == "rung"]

    # ---- steps
    def done(self, step: str, rung: str) -> bool:
        return any(e["kind"] == "step" and e["step"] == step and e["rung"] == rung for e in self.entries())

    def _stamp_path(self, name: str) -> Path:
        return self.folder / "stamps" / (re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") + ".md")      # a colon or a comma in a file name cannot be checked out on every system

    def step(self, step: str, rung: str, payload: dict | None = None, note: str = "", stamp: str | None = None) -> dict:
        """One step of a rung, with its table. Refused when the rung is not open, when the step before it is not on record,
        when the table lacks what the step must hold, and, for a step that draws noise, when no prediction was stamped for
        this rung."""
        if step not in ORDER:
            self._refuse(f"'{step}' is not a step. The steps are: {', '.join(ORDER)}")
        if rung not in self.rungs():
            self._refuse(f"the rung '{rung}' is not open. Nothing of it is run before the person who decides opens it")
        before = NEEDS.get(step)
        if before and not self.done(before, rung):
            self._refuse(f"the step '{step}' of the rung '{rung}' comes after '{before}', which is not on record")
        if step == "report" and self.done("report", rung):
            self._refuse(f"the rung '{rung}' has its report on record. A rung takes one report: the default that follows it is read from that report")
        lacks = [k for k in HOLDS[step] if not isinstance(payload, dict) or k not in payload]
        if lacks:
            self._refuse(f"the step '{step}' is recorded with its table, and the table given lacks: {', '.join(lacks)}. A step with no table is a label")
        if step in DRAWS:
            if not stamp or not self._stamp_path(stamp).exists():
                self._refuse(f"the step '{step}' draws noise. Its prediction is stamped before it is run: stamp(name, text, rung), then step(..., stamp=name)")
            made = [e for e in self.entries() if e["kind"] == "stamp" and e["name"] == stamp]
            if not made or made[-1].get("rung") != rung:
                self._refuse(f"the prediction '{stamp}' was stamped for another rung. Each rung has its own predictions")
            if any(e["kind"] == "grade" and e["name"] == stamp for e in self.entries()):
                self._refuse(f"the prediction '{stamp}' is already graded. A new run takes a new prediction")
        if step == "twin":
            adm = payload["admission"]
            if not isinstance(adm, dict) or "admitted" not in adm or "level" not in adm or _count(adm.get("sets")) < MIN_SETS:
                self._refuse(f"the table of the twin holds the admission of the fit: whether it is admitted, at which level of noise, on at least {MIN_SETS} sets")
        if step == "candidate on the twin":
            if not str(payload["the_planted_candidate"]).strip() or _count(payload["read"]) < MIN_SETS:
                self._refuse(f"a candidate is read on at least {MIN_SETS} synthetic sets that hold it, under its name")
        if step == "report" and len(str(payload["what_it_does_not_show"]).split()) < 3:
            self._refuse("a report says what it does not show, in a sentence")
        n = sum(1 for e in self.entries() if e["kind"] == "step")
        name = f"{n + 1:03d}_{rung}_{step}".replace(" ", "_")
        (self.folder / f"{name}.json").write_text(json.dumps(payload, indent=1, default=float) + "\n")
        return self._write({"kind": "step", "step": step, "rung": rung, "file": f"{name}.json", "note": note, "stamp": stamp})

    # ---- predictions
    def stamp(self, name: str, text: str, rung: str) -> Path:
        """A prediction of one rung, written with the clock before the run that it predicts. It is written once."""
        path = self._stamp_path(name)
        if path.exists():
            self._refuse(f"the prediction '{name}' is already stamped. A stamp is not written again: a correction goes beneath it, with its own clock")
        if rung not in self.rungs():
            self._refuse(f"the rung '{rung}' is not open: a prediction belongs to a rung")
        if len(str(text).split()) < 8:
            self._refuse("a prediction says what will be seen and in which range, in a sentence")
        path.write_text(f"# {name} (stamped {time.strftime(FMT)}, rung: {rung})\n\n{text.strip()}\n")
        self._write({"kind": "stamp", "name": name, "rung": rung})
        return path

    def grade(self, name: str, holds: bool, what_came: str) -> dict:
        path = self._stamp_path(name)
        if not path.exists():
            self._refuse(f"the prediction '{name}' was never stamped. A result is graded against a prediction written before it")
        seen = self.entries()
        if any(e["kind"] == "grade" and e["name"] == name for e in seen):
            self._refuse(f"the prediction '{name}' is already graded")
        if not any(e["kind"] == "step" and e.get("stamp") == name for e in seen):
            self._refuse(f"nothing was run under the prediction '{name}'. The step that it predicts is recorded first, with stamp=...")
        with path.open("a") as fh:
            fh.write(f"\n**Graded {time.strftime(FMT)}: {'holds' if holds else 'misses'}.** {what_came.strip()}\n")
        return self._write({"kind": "grade", "name": name, "holds": bool(holds), "what_came": what_came})

    # ---- admission of a reading
    def _verdict_on_record(self, reading: str, rung: str, level: float):
        """The verdict of the last twin step of the rung that holds this reading at this level of noise, or None."""
        for e in reversed(self.entries()):
            if e["kind"] == "step" and e["step"] == "twin" and e["rung"] == rung:
                table = self.table_of(e)
                cell = table.get("admission") if reading == FIT else (table.get("readings") or {}).get(reading)
                if isinstance(cell, dict) and "admitted" in cell and "level" in cell and abs(float(cell["level"]) - float(level)) < 1e-9 and _count(cell.get("sets")) >= MIN_SETS:
                    return cell, e["file"]
        return None

    def _verdicts_at(self, reading: str, rung: str, level: float) -> list:
        """Every verdict that the twin steps of the rung hold on this reading at this level, in their order. A level that
        was run again until it passed shows here: the admission takes the last verdict and says how many there were."""
        out = []
        for e in self.entries():
            if e["kind"] == "step" and e["step"] == "twin" and e["rung"] == rung:
                table = self.table_of(e)
                cell = table.get("admission") if reading == FIT else (table.get("readings") or {}).get(reading)
                if isinstance(cell, dict) and "admitted" in cell and "level" in cell and abs(float(cell["level"]) - float(level)) < 1e-9:
                    out.append(bool(cell["admitted"]))
        return out

    def admit(self, reading: str, rung: str, level: float = 1.0) -> dict:
        """The twin's verdict on a reading at one level of noise, read from the twin step that the record of the rung holds.
        It takes no table from its caller. For the fit it is the table `admission` of that step. For another reading it is
        the entry of that name in the table `readings` of that step, which holds `admitted`, `level` and `sets`."""
        found = self._verdict_on_record(reading, rung, level)
        if found is None:
            self._refuse(f"no twin step on record of the rung '{rung}' holds a verdict on '{reading}' at {float(level):g} times the noise, on at least {MIN_SETS} sets. "
                         "Run the twin at that level and record it: the admission is read from that record, never from its caller")
        cell, file = found
        had = self._verdicts_at(reading, rung, level)
        return self._write({"kind": "admission", "reading": reading, "rung": rung, "level": float(level), "passed": bool(cell["admitted"]), "sets": _count(cell.get("sets")), "from": file,
                            "verdicts_at_this_level": len(had), "of_which_refused": sum(1 for v in had if not v)})

    def admitted(self, reading: str, rung: str, level: float = 1.0) -> bool:
        last = [e for e in self.entries() if e["kind"] == "admission" and e["reading"] == reading and e["rung"] == rung and abs(e["level"] - float(level)) < 1e-9]
        return bool(last) and bool(last[-1]["passed"])

    def read_records(self, reading: str, rung: str, level: float = 1.0) -> None:
        """The door to the records for a reading: a refusal unless the twin of the rung admitted it at that level of noise."""
        if not self.admitted(reading, rung, level):
            self._refuse(f"the reading '{reading}' was not admitted by the twin of the rung '{rung}' at {float(level):g} times the noise. "
                         "A reading is used on the records only after synthetic records gave it its law")

    def add_term(self, term: str, rung: str) -> None:
        """The door to a larger model: a refusal unless the term, under its exact name, was first read on synthetic records
        that hold it."""
        if not str(term).strip():
            self._refuse("a term has a name")
        for e in self.entries():
            if e["kind"] == "step" and e["step"] == "candidate on the twin" and e["rung"] == rung and self.table_of(e).get("the_planted_candidate") == term:
                return
        self._refuse(f"the term '{term}' was not read on synthetic records of the rung '{rung}'. A guessed term enters the twin before the model")
