---
name: twinfit
description: Fit recorded curves with a ladder of small models and read each with first order algebra and synthetic records before a term is added. Use when a researcher brings fits with parameters that trade against each other, joint fits across a design, or a model that may lack terms, and asks which term to add next or how far the estimates can be trusted.
---

# twinfit

The package is imported from the folder that holds it, `{{WHERE_THE_PACKAGE_IS}}`, and needs numpy and scipy, which `{{THE_PYTHON}}`, the Python that started this project, has: where `python` stands below and the python of this machine lacks them, type that path in its place. Each
command below is run with `PYTHONPATH={{WHERE_THE_PACKAGE_IS}}` in front of it, or from that folder.

1. Read the guide of the package before anything else: `python -m twinfit guide` prints it, the file `AGENT.md` of the package. It holds the order of the
   work, what to tell the researcher at the start, and the rules that hold at every step. The tutorial, `TUTORIAL.md`
   in the folder of the package, runs the loop from nothing on two worked projects.
2. The project was started with `python -m twinfit init .`: `description.toml` describes the data collection, the
   curves are in `traces/`, and `model.py` is the first model. Fill them with the researcher, or read them as given.
3. Run `python -m twinfit loop .`. It fits the model, reads what it leaves against its twin, and writes
   `loop/model_0.md`: where to search, of which kind, along which setting, and three proposals for the model and three
   for the experiment.
4. Write in `terms/01.py` a term of physical origin for what the report points to, with the reason in its `why`, a
   sentence of six words or more, and run the loop again. `python -m twinfit terms` lists the terms of the package and
   the form of a term file, and `python -m twinfit scan . "second(x)" 0.5 1 1.5` fits a term whose distance, stretch or
   period is chosen at each value. It keeps or refuses the term, and stops at the next file that is not written. The output is
   `loop/best.md`: the best model and five parts of the uncertainty of each parameter.
5. Each iteration follows the ten steps of `AGENT.md`: read the structure of the residuals, search the literature,
   estimate the size of each effect, choose the three most promising, let the scientist choose, refine the model, ask a
   subagent to red team it, run the joint fits again, report the improvement, and present three ways to improve the design of the experiment. Begin with the small model that the
   researcher understands in every part.

What the package cannot do is section 0 of `AGENT.md`. Tell the researcher before the first run. When they are not
there, it goes at the head of your report.
