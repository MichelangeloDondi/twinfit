# This project is analysed with twinfit

The rules below hold in this folder. They are short on purpose: the order of the work and its reasons are in the
protocol of the package, which `PYTHONPATH={{WHERE_THE_PACKAGE_IS}} python -m twinfit guide` prints. Read it at the
start of every session. The package is the folder `{{WHERE_THE_PACKAGE_IS}}/twinfit`, and it needs numpy and scipy, which `{{THE_PYTHON}}`, the Python that started this project, has: where `python` stands below and the python of this machine lacks them, type that path in its place.

1. **Small models first.** The analysis starts from `model.py`, the small model that the researcher understands in every
   part, and `python -m twinfit loop .` adds one term a rung, each written in `terms/` with the reason for it.
2. **The notebook is the record.** The loop writes every step to `loop/notebook` through `twinfit.notebook.Notebook`, the
   report of each model to `loop/model_N.md` and the output to `loop/best.md`. A step that the notebook refuses is not
   done another way.
3. **The loop does the twin.** It tries each term on its own synthetic records, with a prediction stamped before, and
   keeps it only when the twin admits it. Do not stamp or call the twin by hand for a project that the loop reads.
4. **Read the tutorial and the tools first.** `TUTORIAL.md` of the package runs the loop from nothing.
   `python -m twinfit terms` lists the terms and the form of a term file, and `python -m twinfit scan` fits a term whose
   distance, stretch or period is chosen at several values.
5. **A term is a file with its reason.** `terms/NN.py` holds `term` and `why`, a sentence of six words or more that
   names the physics, its source and the answer of the report that it fills.
6. **Report a set of candidates as a set**, each with the bias that it would mean, and say what cannot be told apart.
7. **Say what a result does not show.** A cause is named only as a candidate in a term's `why`, which the loop and its twin bear out or refuse, and never as found. Claim nothing new.
