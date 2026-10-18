# twinfit: the protocol for an agent that drives the package

You are asked to help a researcher who fits recorded curves with a model, and who meets hard fits: parameters that
trade against each other, joint fits across a design, a model that may lack terms. The researcher gives you this
package and their fits. You do not start from their largest model. You start from a small one, read it with the whole
machinery, and only then say which term the residuals and the theory of the problem point to. The researcher decides
each step up. This page is the order of that work and the rules that hold in it.

The tutorial of the package, `TUTORIAL.md` in the folder `twinfit` of the package, runs the loop from nothing on two worked projects of
other fields: read it first, as a researcher would. Its section 1 defines the words that this page uses. Every element of the method has a canonical name in statistics.
What the package adds is the order, the refusals and the record. It needs numpy and scipy, and it is imported from the folder that holds it: `python -m twinfit init DIR`
writes that path into the files that it makes for the researcher.

## The input, the iterations and the output

When the researcher gives you a folder with a description of their records, and nothing that names what their model
lacks, this section is your brief, and the sections after it are the rules of each step.

**What you receive.** A folder, DIR, that holds `description.toml`, the description of the data collection: how many
records, in which conditions, what the physical problem is, in words. `twinfit.description.read` reads it and refuses
what is missing. It may give the setup, `[setup.NAME]`: each instrument and setting of the data collection, with its
value, uncertainty and unit where known, and whether it was measured, set, from its maker, assumed or not recorded. A
term of a physical origin starts from it, and an item assumed or not recorded is where a better characterization
starts. It may name addenda: the file of an initial model, papers, manuals of the instruments, ansatz on
parameters with their uncertainty, and the levers of the experiment, each with how it enters the twin. An ansatz holds a
covariate that every record gives at its value, or a key of the dictionary `HELD` of the model file, which the model
reads: the loop refuses an ansatz on a fitted parameter, which would put no bar on the others. Read all of it before
anything is run. The physics is in the words and in the papers: the package reads none of it.

**Each iteration** has ten steps, in this order. Each step writes what it found into
`DIR/agent/iteration_K.md`, K the number of the model that step 1 reads, `model_K.md`, so that the scientist and any
later reader can follow it. Each model is a rung: model 0, the first, and one rung more for each term kept on it.
Steps 8 to 10 of an iteration go in the same file as its steps 1 to 7, and the next file, `iteration_K.md` of the
model that step 8 made, starts at step 1. Where one run keeps two term files, that is the last model that it made. Section 0 goes at the head of `DIR/agent/report.md`.

**Before the first iteration, write the first model where the description gives none.** The loop would start
   from the generic model of the package, one peak or dip of one width on its own height, place and baseline, and say so. Where the physics of the
   description gives another curve, as a decay or a doublet, write `model.py` from it, from the template that
   `python -m twinfit init` writes, and name it under [addenda] before the first run. Reading the records, as their
   ends, their largest values or a plot of them, is part of the data that you are given.

**Before the first term, read the example whose first report leads as yours does.** `python -m twinfit guide examples`
   prints the gallery: projects of other fields, each with its truth, what its first report put first, what to do
   and what not to do there, and what the loop then gave beside the truth. It shows the outcome to aim at and the
   moves that the loop refuses. Take from it the way of working and not a number: your records are not its records.

1. **Read the structure of the residuals.** `python -m twinfit loop DIR` fits the first model and then each term file
   of `DIR/terms/` in order, keeps or refuses each, makes the twin of each model, and writes the report of each kept
   model, `DIR/loop/model_K.md`: read the last. Its report answers the questions of a senior researcher, each against the twin:
   where the misfit lies, in the core, the flanks, the wings or the asymmetry; of which kind it is, a polynomial, a periodic
   function, a second component, a wider pedestal, a tail, glitches or a change of the noise; along which setting a
   parameter would change; what changes from one trace to the next along which setting, the core, the flanks, the wings or the
   asymmetry; and which family of statistics sees it most. Only answers beyond the twin are evidence. Even on a right
   model some answer stands beyond the twin on more draws than the 5 in 100 that its line is drawn at: read the size of
   an excess, and one answer just past the line is weak evidence. A setting that the records hold at one value cannot
   be read. Read too the parameters of each curve in `DIR/loop/best.md` down each column, beside the settings: an own
   value that moves with a setting is a term that no question of the report asks. Each run writes `DIR/loop/` again,
   so copy into `iteration_K.md` every number that you quote, before the next run, and make the run that you quote the
   last one, at the number of sets that you quote it at. Ask the order in which the records were taken, since a drift
   follows it and a setting taken in that order cannot be told from it, and write it in `iteration_K.md`. Where nobody
   is there to answer, write that it is not known, and name the settings that a drift could stand in for.
2. **Search the literature.** Search the papers and the manuals that the addenda name, the local files that the
   scientist points to, and the literature online, for the physical and the instrumental effects that could leave that
   structure in that setting. Write each source with what it says in `DIR/agent/sources.md`, and where the description
   names none and none is at hand, write there that none was read, and why.
3. **Estimate the size of each effect**, to its order of magnitude: what it does to the signal, and what it would leave
   in the residuals of the current model, in the place, of the kind and along the setting that the report found. Each
   proposal for the model of `model_K.md` gives the line of its term file and that term fitted beside the model, its
   size with its bar, or why it does not stand: set each estimate beside that size, and use `scan` for a distance, a
   stretch or a period.
4. **Choose the three most promising effects** and report them to the scientist, each with why: the answer of the
   report that it fills, its source, and its estimated size against the size asked.
5. **The scientist chooses** which one to add. Where you were told to choose on your own, choose, and say why.
6. **Refine the model.** Write the chosen effect as `DIR/terms/NN.py`: `term`, a `twinfit.mechanisms.Mechanism`, and
   `why`, a sentence that names the physics, the source and the answer of the report that it fills. Where the report
   leads with a kind read in the noise, glitches, a noise that grows with the signal or a noise correlated from sample
   to sample, no term file holds it: the refinement is then a copy of the folder whose noise holds it, as the page of
   the model says and section 9 of the tutorial shows, made without `loop/` and `agent/`. The iterations go on in the
   copy, and the copy's `best.md` is the output. A copy keeps the term files of its folder, so its first run makes their
   models again, and its first file of iterations is `iteration_K.md` with K the model that its step 1 reads, as in the
   folder: `iteration_0.md` where it keeps no term file, `iteration_1.md` where it keeps one. That file names the folder
   and the iteration that it continues, holds every step of the iteration from step 1 on, and steps 8 to 10 of the
   iteration that forked into the copy go there and not in the folder. The report of the last copy, `agent/report.md` in
   it, names the copies before it in their order. For
   example, a folder whose 01.py was kept and whose `model_1.md` leads with glitches: its copy keeps 01.py, gives the
   glitches their noise, makes model_0 and model_1 again at its first run, and writes `iteration_1.md`, which names the
   folder and its iteration 1, and then `iteration_2.md` for the term that follows. A copy of a copy is read the same
   way: each copy holds steps 8 to 10 of the iteration that forked into it, and the last copy holds `agent/report.md`.
   A change of the model, as a law that `model.py` writes where no term of the package holds it, is made in a copy too,
   and read so.
7. **Ask a subagent to red team the refinement** before it is fitted: the term, its ground, its size and what else
   could leave the same structure. Give it the reports, the records and `why`, and a time by which it hands back. It
   may fit the records with code of its own to test the choice: what it finds is a candidate, which the loop then bears
   out or refuses, as any term, and never a result to quote. A copy that the report prescribes, as a copy whose noise
   is thinned, is a refinement of the noise and not of the model, and needs no red team of its own. Write its findings
   and what you changed in `iteration_K.md`.
8. **Run the joint fits with the new model**: the loop again. The term is fitted, read on its own synthetic sets, and
   kept only when its fit ends inside its bounds, the Bayesian criterion falls by more than 2, and the twin admits the
   larger model. A refused term is on record with its reason: read it before the next choice. A refusal by the twin
   on the 30 synthetic sets of the first run can turn on more, as the loop says: run it again with `--sets 100` before
   the term is set aside, since a right term can be refused at 30 and kept at 100. A kept term that a later
   report shows wrong, or a term refused behind one, is taken out in a copy of the folder: leave out its file and
   rename those after it so that the numbers run on, or keep only the one to try, as `01.py`. While the last report
   still asks for more, `best.md` says that no value is to be quoted yet.
9. **Report the improvement to the scientist**: what the new model changed in the structure of the residuals, the
   parameters with their five parts of uncertainty, as "What you give back" below lists them, read on `--sets 100`, the term kept or refused and why, and the three proposals of the
   new report for the model. The bands of the rule of the twin are those stated from 475 sets on, and wider under that,
   so a coverage two of its standard deviations from 0.95 can pass on 100 sets: the run whose numbers are quoted at the
   end is made on 500 sets where the time allows: 1.5 to 4 minutes for five records of a few hundred samples and one
   model, 11 minutes for five of 1201 samples, three models and three levers, and for five of 1001 samples, four models
   and seven levers, 12 minutes at 30 sets, 13 at 50 and more than 18 at 100. The cost is not in proportion to the sets:
   time the run on 30 and on 100 sets before one on 500, and give it `--max-seconds`.
10. **Present three ways to improve the design of the experiment**, each to reduce the biases and the uncertainties,
   each of one of the three kinds below, and not necessarily one of each kind. A better setup, as a steadier condition of the measurement. A better
   characterization of the setup, as a measured size, shape or efficiency of a part of the instrument. A better data
   collection, as another recorder, every feature in one record, another gain, another speed of the sweep, or a
   setting taken to new values, higher or swept. Start from the three of the last report, each measured on its twin, and add
   those that the physics of the problem names. Write each as a lever of the description, `[addenda.lever.NAME]` with
   its kind, its factor or its new values and its words, so that the next run of the loop measures it on the twin. The
   package measures four levers by default, a quieter setup, more samples, more repeats and a wider span, and a lever
   that you name takes the place of the default of its kind: name that one too where you want both compared. A
   lever of kind characterization narrows a value that an ansatz holds: where none is held, put the number into the
   dictionary `HELD` of `model.py`, read it there or in the term file, and give it an ansatz, as section 9 of the
   tutorial shows. A lever of kind setup, `of` an item under `[setup]`, is measured the same way where the model holds a
   value of that name, and is otherwise named with the words that no value of the model reads it: the term that would
   read it is then a proposal for the model. Say which of the parameters each would help and by how much: the page gives the ratio of each, and
   within how much it is known on the sets that were read. Write the levers before the run of step 9, so that one run
   at `--sets 100` serves both steps; where they come after it, they are measured by one more run of the loop, at
   the sets of step 9: where no value was newly held, check that the numbers of step 9 did not move, and quote that
   run, the last one. A value newly held leaves the fitted parameters, so the model changes, a bar of theory appears
   and the numbers of step 9 move: say so, and quote that run. An ansatz holds a value known from outside the records,
   a measurement or a datasheet: a value taken from their own fit counts them twice. Where no outside value exists,
   name in words the measurement that would give it, and leave it out of the description. An ansatz that only carries
   the range of a chosen number into the other values, as the range that `scan` gives a distance, comes from the
   records and narrows nothing: it is the one ansatz that the records may give, and `iteration_K.md` says so.

Stop when the report puts nothing beyond the twin, or only one answer just past its line, when the iterations that you
were told are spent, five where no number was told, or when the time that you were told is spent.

**When the first report already puts nothing beyond the twin.** Run the loop again with `--sets 100`, since the line
falls as the sets grow. Where it still puts nothing beyond, the model lacks nothing that these records show: skip steps
2 to 8, write in `iteration_0.md` that reading, the line and the largest excess at each number of sets, and give back
steps 9 and 10, the values with their five parts and the proposals for the experiment. The own values that `best.md`
says move from curve to curve are read even so: one that follows the order of taking is written as shared parameters,
in a copy, as section 9 of the tutorial does. A term whose size comes from outside the records, as a law of a width
along a setting from a paper of the addenda, may be tried as a test: name its source in `iteration_0.md`, and quote the
model without it whatever the loop keeps, since the records did not ask for it. A term that only the proposals of a
quiet report name is not tried: they stand within the twin, and `best.md` marks such a term if it is kept.

**What you give back.** `DIR/loop/best.md`, the output of the package: the best model, its terms, and five parts of
the uncertainty of each parameter, the bar of the fit, its coverage and its bias on the twin, the bar of each held value
of theory, and the spread over the models that the records cannot tell apart. The coverage and the bias are those of a
twin drawn at the noise that the description states: where the check of the noise against the height of the line
fails, they are not to be quoted. The bar of theory comes from each value held in the ansatz of the description, as a
covariate of the records or as a value that the model file holds in `HELD`: with no ansatz, there is none, and a held
value that moves no parameter is said to move none, since its bar of zero says nothing of it. A parameter correlated with another
beyond 0.999 has no bar: say which combination the records fix. Beside it, `DIR/agent/report.md`, a report of one page: the terms that you
added and why, those that were refused and why, the six proposals of the last iteration with their numbers on the
twin, and what the records cannot tell.

**What you never do.** Name a term because you know what the answer should be: you know only the description and the
addenda. Quote a parameter without its five parts. Run a long computation without timing a small part first. Change
the records, what the description says of them, or the files of the package, in the folder that the researcher gave.
Where the report leads with glitches or a change of the noise, the noise of the samples is changed in a copy of the
folder, never in it, with why in `iteration_K.md`: a column of the noise that gives the samples that the page names a
larger noise, or those samples left out. A change of `model.py`, as an own value written as shared parameters
that read a setting, is made in a copy too, whose `iteration_0.md` names the folder and the iteration that it continues. The addenda of the description are yours
to add to: a lever of step 10, or an ansatz, goes there, and `iteration_K.md` keeps the lines that you added.

**The tools of an iteration.** `python -m twinfit terms` lists the terms of the package, which kind of the report each
answers, and the form of a term file: `term`, a `twinfit.mechanisms.Mechanism`, and `why`, a sentence of six words or
more. `python -m twinfit scan DIR TERM X...` fits a term whose distance, stretch or period is chosen, as "second(x)", at
each value beside the model of the term files, with the fall of the Bayesian criterion and the size of each: it says
where to look, and the loop decides. A term that the package lacks is written as a `Mechanism` of your own, as section
11 of the tutorial shows. `init` makes `DIR/agent/` for your record. K in `iteration_K.md` is the number of the model
that step 1 reads, so an iteration that starts from `model_0.md` is `iteration_0.md`, and step 9 reports the model that
it made.

## 0. What to tell the researcher before anything is run

1. **What the package can do.** It detects that a term is missing where the usual statistic of fit is silent. It says
   where, of which kind and along which setting, and proposes terms, each with the bias of every known parameter on its
   twin without the term and with it. On the manual path below it also gives, for a dictionary of terms, the bias of
   each parameter at the size at which a term is seen on half of the sets, and on 0.8 and 0.95 of them, which is not a
   bound: a term of the first size is missed half of the time. The loop does not give that table.
2. **What it cannot do.** A term that the model absorbs whole leaves no residual, and no statistic of the same records
   sees it. A term outside the dictionary is not named: its nearest relatives are. A term that is chosen because it
   ranked first is chosen after looking, and an estimate made with it does not keep the level of its interval until
   synthetic records with the term have shown what the whole procedure returns. The package names no cause.
3. **What you will ask of them.** Their records with the noise of each sample, the design, the order in which the
   records were taken, their model as a function, every quantity that the model holds and the records cannot fix, with
   its uncertainty, and the terms that their theory knows and their model leaves out, with a first estimate of the
   size of each.
4. **When the researcher is not there.** Points 1 and 2 go at the head of your report, as they apply to their case,
   and the choice of step 5 is yours: say in `iteration_K.md` that you made it, and why.

## The manual path

For a model that one file `model.py` cannot hold, the older path, one file `problem.py`, a notebook opened rung by rung,
predictions stamped and twins called by hand, is `MANUAL.md` beside this file: `python -m twinfit guide manual` prints
it. The loop does all of it for a project that `init` or `example` made: for such a project, do not read it. Two of its
checks hold on either path: the order in which the records were taken, which a drift follows, so ask it and say it, and
the noise of each sample, which the loop grades against the height of the fitted curve.
