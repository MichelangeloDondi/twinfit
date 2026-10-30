# twinfit: the manual path beyond one `model.py`

The protocol of an agent is `AGENT.md` beside this file, and the loop that it runs needs none of this page. The sections
below are the older path: one file `problem.py`, a notebook opened rung by rung, predictions stamped and twins called by
hand. Read them only for a model that one file `model.py` cannot hold. Section 0, what to tell the researcher before
anything is run, is in `AGENT.md`. Not covered here: the loop's own protocol, which `AGENT.md` gives.

## 1. Intake

Write one file, `problem.py`, from `agent/problem_template.py`: the records, the model, the names and the starts of
the shared and of the own parameters, the center of a curve when it has one. Then check, and report each check:

- Every record has a positive standard deviation for every sample. If the researcher has none, estimate the noise from
  a part of each record that holds no signal, say that you did, and say what it leaves out.
- The fit that the researcher brings is a point where the fit of this package stops, to within a tenth of a standard
  deviation of each parameter. If it is not, say so before anything else.
- The order of acquisition is known from the time stamps of the records and not from the order of a list. A design
  axis that rises or falls with the order of acquisition cannot be told from a slow change of the conditions of recording: say which
  axes are in that case.

No call of the package makes these three checks. They are yours: write the result of each into the table of the
intake. In the first trial of the package the one axis of the design rose with the order of acquisition in every
record, and only the agent's own check said so.

**An abscissa read from a file is rounded.** Its steps then differ in their last digit, and the terms that apply a kernel refuse it. Build each record with `twinfit.even(x)`, which returns the same span in equal steps and refuses a difference larger than rounding.

**The abscissa of a record rises.** `Record` refuses an abscissa that does not rise strictly from its first value
to its last: the widths, the periods and the bounds of a center are counted from its ends and its steps. Reverse a
record that was taken from high to low, with its values and their standard deviations, before it is read.

**The package holds no list of terms.** The candidates that `propose` ranks are yours, each a `Candidate` or made
from a term of `mechanisms.py` by `from_mechanism`, written from the theory of the problem.

**The covariates.** Every key of the design whose value is a number with more than two values is read as a covariate:
the number of a repeat and the place in the order of acquisition too. `covariates_of(problem)` names them. Keep an index
out with `covariates_of(problem, covariates=[...])`, and say in the report which keys were read.

## 2. The ladder

Ask the researcher which model is the simple case that they understand in every part. That is the first rung. Write
beside it the terms that their theory knows, in the order in which they would add them, one term a rung. A rung is
opened by the researcher, and their words go into the notebook: `book.open_rung(name, by="their sentence")`. The
notebook refuses a rung that nobody opened.

**Who opens a rung, said once.** While the researcher answers you, they open
each rung: propose it and wait. When they are not there, because they said to go on without them or because nobody
can answer, the first rung is opened with the words that they left and each later one by default:
`book.open_by_default(name, after=rung)`.

The call takes nothing from you. It reads the term and what the choice rests on
from the report of the rung before, so record that report with the choice as `decide` gives it, under
`default_choice`. The notebook marks such a rung as opened by default, and `book.opened_by_default()` lists them: the
researcher sees on coming back which rungs they opened and which the rule did, and can stop or go back at any report.

A choice of theirs at a report, where it is not the default, is recorded in their words with
`book.choose(rung, term, how, by)`, and `book.chosen(rung)` returns the choice in force with who made it. The rung
after such a choice is theirs to open: the notebook refuses to open it by default. A rung that goes back restates
the model of an earlier rung, and the rung that was left stays on the record with what refused it.

A rung takes one report. If something changes after a prediction is stamped and before its result, such as a run
that will not be made, write it beneath the stamp with `book.note(name, text)`: it becomes an entry of the record.

## 3. The cycle of a rung

Each step is a call, a line in the notebook and a short account to the researcher. The notebook refuses a step while
the one before it is not on record, a step with no table, and a step that draws noise with no prediction stamped
for its rung before it. It writes every refusal into the folder. What it cannot check is that a table is true, since
it records the table that you give it: that stays your word, which the folder lets a person audit.

Every step is recorded with its table: `book.step(step, rung, table)`. A prediction belongs to a rung:
`book.stamp(name, text, rung)`. The two steps that draw noise name their prediction:
`book.step("twin", rung, {"admission": table, "readings": {...}}, stamp=name)`, where `table` is what `admission`
returned and `readings` holds, for every other reading that you will use on the records, whether the twin admits
it, at which level and on how many sets.

A reading is then opened for the records by `book.admit(reading, rung,
level)`, which takes no table from you: it reads the twin step that the folder holds, at that level, and refuses when
there is none. Then `book.read_records(reading, rung, level)`. A term is opened for the model by
`book.add_term(name, rung)`, which refuses a name that no step "candidate on the twin" holds.

| step | call | what it gives | what you say |
|---|---|---|---|
| intake | `fit(problem)` | whether the search ended, the statistic of fit | the checks of section 1 |
| identify | `identify(Linear(problem, result), block=...)` | the first order bound of each parameter, the weakest, the block that each rests on | what the design cannot fix, before any residual is read |
| twin | `noiseless(problem, result)` first, then `book.stamp(...)`, then `admission(problem, result, sets=..., level=...)` at each level of noise, then `planted(problem, result, [[None] * len(problem.records)], [0.0], candidates=...)`, then `book.admit(...)` | whether the fit returns the truth, and how often the table of candidates cries wolf | the prediction, then its grade |
| fit | `book.read_records("the fit", rung)`, then `book.step("fit", rung, table)` | the estimates, admitted | values with their bars |
| read | `slices(lin)`, your own candidates | the candidates to be ranked | which terms of the theory you wrote as candidates |
| propose | `propose(lin, candidates)`, `lines(table, show=[names of the terms of the theory])` | the leading set, the range of the bias, what cannot be excluded | section 4 |
| candidate on the twin | `book.stamp(...)`, then `planted(problem, result, columns, amplitude, name=..., candidates=...)`, recorded with `stamp=` | what the procedure returns when the term is really there | whether the readings of the records are those of synthetic records with the term |
| report | `book.step("report", rung, {"what_it_does_not_show": ...})` | the rung, closed | section 5 |

**Several workers.** `admission`, `planted`, `null_gains` and `on_the_twin` take `workers=`: the sets are read in
that many processes and the tables are the same. Where one fit takes ten seconds, 200 sets with one worker are half
an hour. Set `VECLIB_MAXIMUM_THREADS=1` before the run, so that each process keeps to one core.

A worker can be lost: on some platforms a forked process cannot use the threads of the library of linear algebra and
ends at its first large matrix. The package asks each worker for one thread where it can. If a worker is lost all
the same, the call closes its pool, computes the missing rows in the calling process and warns
(`RuntimeWarning: twinfit.pool`). The tables are the same and the time is that of one worker: say so in the report of
the stage.

A job that raises in a worker is run again in the calling process, so that its exception comes out as it
is. If it does not raise there, the fault was the worker's alone, and the call warns as for a lost worker. The interpreter may warn at each fork that a process with several threads can deadlock in its child: it names
the same hazard, and the package does not hide the warning.

**A search with bounds that stops is tried once with none.** If the second search ends inside the bounds at a lower chi square, its point is the fit, the result holds `finished_without_bounds`, and the report of the stage says so. It finishes some stalls and explains none: say in the report when it happened.

**A search that did not end.** The limit is 80 steps, and at the noise of the records a model with records that
hold almost no signal can need more. On the sets where this was looked at, what still moved was the center of a
record whose peak is under its noise, after the shared estimates had settled. Ask for more,
`fit(problem, max_iter=400)`, and the same `max_iter` in `admission` and `planted`. Say how many steps it took, and
say so when a search still has not ended: its shared estimates may be right and its report must not be read as if
it had.

A search can also stop before its limit, when no step lowers the chi square any more: `develop_one` says
which of the two it was under `why_not`. More steps do not help the second. Report the term as not fitted, with
that sentence.

**A parameter held at a bound.** It has no bar. `admission` does not admit a fit in which a parameter is held on
more than one synthetic set in twenty, and its table says why under `why_not`: move the bound away from the
estimate, or report that parameter as one that the bound sets. Under that share the cell says on how many sets it
was held.

**A term that bends the abscissa, with a kernel in the model.** It can be added outside the kernel: the kernel is
then read on equal steps that do not move with the bend, and taken back to the bent points by a cubic. The first
repair made the curve right and left a kink in the chi square at each point of its grid, where a search along
the bend stopped.

**The prediction comes first.** Before any replica is drawn, compute what first order gives, write the prediction with
its range through `book.stamp`, and say for each clause whether it could miss, repeats a number already seen, or holds
by construction. Grade it afterwards with `book.grade`. A miss is recorded beside the stamp and is a result.

**The ladder of noise.** Run the twin with no noise first, `noiseless(problem, result)`: started at the starts of the
problem, the fit must end and return what was put in to many digits, and with `planted=` it gives what a term does
by itself. It cannot show an error of the model function, which the records and the fit share.

Then run `admission` at a tenth, a third, one, three and ten times the recorded noise. Its searches start at the starts of the problem too, as
the researcher's fit does. Say which starts you used: a fit with several minima is admitted or refused by them. A reading is used on the records only at a level at which
the twin admitted it: `book.read_records` refuses otherwise.

**How the ladder is recorded.** One twin step for each level of noise, all under the one prediction, which is graded
after the last. `book.admit(reading, rung, level)` reads the last twin step of the rung at that level, and its line in
the notebook says how many verdicts the level had and how many of them refused.

**The bands of the rule, and a level that was refused.** A fit is admitted when its searches end, its bias is under
0.25 of its spread, its spread is 0.85 to 1.15 of the bound and its interval holds the truth on 0.92 to 0.98 of the
sets. With few sets a band is made as wide as three standard errors of its statistic, and the table says so: the
three are held as stated from 475 sets on. Use 500 sets or more where the fit is fast.

A level that the twin refused is run again once, on at least three times the sets, under a new prediction that gives
the arithmetic of the first refusal. Report both verdicts. A third run is a search for a pass.

**The list that the twin reads.** `candidates=` is a function of the linearized fit of a replica that returns the
list of candidates that you read on the records: the slices and the terms of the theory. A rate of false alarms holds
for the list that it was counted on. A planted term whose name the list does not hold is refused, and so is a twin with
no list: the package holds none.

**Nuisances.** A quantity that the records cannot fix and that is known with an uncertainty is drawn from its
distribution in every replica and held at its central value in the fit: pass `nuisances=` to the twin. Say whether one
value is drawn for each record or one for the whole set: the first averages over the records and the second does not.

## 4. Guessing the next term from residuals and theory

1. **Describe before you name.** From the table: is something seen at all. Is it even or odd about the center. Where
   along the abscissa does it sit. Along which axis of the design does it grow. Does it sit in one block.
2. **Read the leading set as a set.** Report every member with the bias that it would mean, and the range of that
   bias for each parameter. The members can give biases of opposite sign. Never report the first member alone.
3. **Bring the theory.** For each mechanism that the theory of the problem knows and the model leaves out, write what it
   adds to the mean of a record as a function, with its expected size, and add it to your list as a `Candidate`.
   Rank again. A mechanism whose candidate is absorbed whole cannot be tested by these records: say so.

   **A term with a solid ground in the theory of the problem is added to the model whatever a ranking or an
   information criterion says of it.** A criterion weighs a gain of fit against a count of parameters and knows
   nothing of the system. Ask the researcher which mechanisms stand on such ground and take those first, one at a
   time: each still goes into synthetic records before it goes into the model, as point 4 says, each in a rung of its
   own, opened as section 2 says, and each with one free multiplier of the size that the theory gives. Report what the records
   make of that multiplier and what the term costs the other parameters in precision.

4. **Try it on the twin before the model.** Plant the candidate in synthetic records at the size that the residuals ask
   for. Compare what comes back with the records: the score, the leading set, the readings by slices, the estimates.
   A candidate that reproduces the readings is not excluded. That is all that can be said of it.
5. **Propose, and wait if the researcher is there** (section 2 says who opens a rung). Give them the candidate or the
   set, the evidence, the range of the bias, what the records cannot tell apart, and what the next rung would cost each
   parameter in precision. The same proposal goes into the report in every case.
   `book.add_term` refuses a term that was not tried on the twin.

## 5. The report of a rung

Short, in plain words, the answer first. Three columns in every comparison: the canonical method at its strongest, the
readings of this package, and the two together. Say where the readings add nothing. Give each parameter with two
uncertainties: the one that the noise and the nuisances give, from the twin, and the one that the model may hide, from
the range over the leading set and from the line of what cannot be excluded. End with what the rung does not show and
with what is the researcher's to decide.

## 6. Rules for every step, each from a past fault

- **"Is not excluded" and "fits" are two statements.** Write the second only when a computation tested it, and name
  the file.
- **Nothing is chosen after looking and then reported as found.** A record or a block that stands apart was chosen by
  its distance: say so, and give the statistic over all of them.
- **Each case with its sign and its size.** "Cannot be told apart" is written with what each case would mean.
- **Two methods are compared on the same samples**, and the canonical one at its best setting.
- **A number keeps its condition and names its file.** An uncertainty has two significant digits.
- **One change at a time**, with no noise first.
- **A refusal is code and has a test that tries the careless way around it.** A sentence that asks is not a refusal.
  A refusal that reads a label and not the thing is not one either: that fault was made three times here.
- **A correction is made where the sentence stands.**
- **When the researcher keeps to models that are symmetric about the center**, call `propose(..., symmetric_only=True)`:
  the odd candidates are alarms and none is proposed as a term.
- **A table that sees nothing has excluded nothing.** On synthetic records, two hidden terms of five were seen on one
  set in ten or fewer, and one of them moved an estimate by 4.7 standard deviations. Report, for the term that the
  model absorbs most, the bias at the size at which it would be seen on half of the sets and on 95 per cent of them,
  and say that it is no bound and holds for that list alone. Its size depends on how many candidates the list
  holds, since the threshold of the family does: in the first trial one term had 10.5 standard deviations in a list
  of 77 and 6.1 tested alone. Name the list with the number.
- **The range of bias over the leading set is not an interval.** On synthetic records, for the three hidden terms
  that were seen on most sets, it held the true bias on 22 to 100 per cent of the detections. Report it as a
  disagreement between candidates, never as an uncertainty.
- **The first of the table is not the cause.** For the three hidden terms that were seen on most sets, a candidate of the
  hidden kind was first on 43 to 73 per cent of the sets. Every candidate goes into synthetic records before it is added to the model.
- **A verdict names its yardstick.** On one dataset the departure of a block had a level of 0.095 against its own
  scatter, 0.013 against a pooled scatter, 0.4 for the paired step from the block before, 0.0000016 against the noise
  alone, and 0.004 once the member chosen as the farthest was dropped, which is a selection. The first three are fair
  and the last two are not. Give each with its name.
- **No cause is named, and no claim of novelty is made.**
- **The notebook is the record.** A person must be able to audit the analysis from its folder with no chat.

## 7. Files and their contents

| file | what it holds |
|---|---|
| `problem.py` | records, model, parameters |
| `description.py` | the description of a problem in one file, its addenda and the levers of its experiment, read with their refusals |
| `where.py` | where to search a new term: the place, the class, the covariate and the family of statistics, each against a twin of the fit |
| `best.py` | the output: the best model after the iterations, its terms, and five parts of the uncertainty of each parameter |
| `proposals.py` | three proposals for the model and three for the experiment, each measured on the twin |
| `loop.py` | the loop without an agent, `python -m twinfit loop DIR`, whose term files you write |
| `fit.py` | the canonical fit of every record at once. A problem may give bounds (`bounds_shared`, `bounds_own`): give them whenever some records hold almost no signal, an area that is not negative and a center inside the abscissa. Without them the own parameters of such a record wander. Read three numbers of every fit before its estimates: `held`, the own parameters at a bound, `undetermined`, the own parameters that their record does not fix, and `at_bound`, the shared parameters at a bound, whose bars are not bars |
| `shapes.py` | a curve from its characteristic function: for a model that is a convolution of several kernels, or a mixture of such convolutions, and has no formula |
| `linear.py`, `ledger.py` | first order at a fit: the absorbed and the visible part of any discrepancy |
| `identify.py` | what the design can fix |
| `families.py` | families of statistics of a curve, each as a rule of weights: generalized moments, windows along the abscissa and along the value, the two sides of a curve, the abscissa against the value, a transformation of the value. For a family at a fit: what it fixes by itself (`bounds`), what it reads in the residuals (`reading`), its bias under a missing term (`shift`), its own estimate (`estimate`), how much of each own parameter it sees (`seen`), and its estimate set against the fit's (`contrast`, with `draws` for a level that can be trusted). Every family of the catalog at a fit: `table`, and `planted_table` with each named term put into records of no noise. A combination of estimates in which the first order bias of named terms cancels: `free_of`. One weighted by the covariance of the estimates on synthetic records of the fit, with its standard deviation read on further sets: `together`. Which of the cells of objects and techniques exist: `cells`. What the module refuses, each with its reason in the error: a fit that did not end, a hold of an own parameter that the family sees at over 0.05 of its column, an estimate from a family that leaves a direction of the shared parameters free or that is at 100 times the bound of the fit, and a contrast whose search did not end. A first order level of a contrast comes with `first_order_verified`, which is False outside the range where that law was found to hold: read the level with `draws` there. An own parameter that a family sees by rounding alone is not estimated by it, and `blind` names it. The choice is made record by record: `free_at` says where each is estimated |
| `losses.py` | chi functions: the fit under a cost of the residual that is not its square (`chi_fit`): a power of the residual, Huber's cost, a logarithmic cost, the likelihood of a Student law. Each with the factor of its first order standard deviation over the one that the weighted fit prints, and `converged`: a search that did not end has no standard deviation. A power under 1.5 or over 3 is refused. The residuals of a fit counted on the plane of position and value, as a test of independence (`grid`), and that reading under each cost (`crossed`) |
| `standard.py` | the checks of residuals that an analyst runs with no method of this package, to stand beside its readings in a table: the correlation of neighboring residuals (`serial`), the runs test (`runs`), the chi square of binned residuals (`binned`), and the score test against a named term (`against`), which at first order is the likelihood ratio of the model with that term added |
| `matched.py` | statistics that are formed on the noisy values of a record, the crossings of levels, the counts in bands of value and the logarithm, each matched to its mean on synthetic records of the model (`match`) |
| `candidates.py`, `readings.py` | the form of a candidate, and the pairs of slices |
| `propose.py` | the table: scores, the leading set, the bias of each member, what cannot be excluded |
| `twin.py` | replicas, admission, a planted term |
| `mechanisms.py` | a guessed term as the working model with one parameter more (`extend`), or held at the size that theory gives (`hold`): a kernel with exponential wings, a displacement of linear density, a parameter that changes along a covariate, a pedestal, a bent abscissa |
| `refine.py` | the candidates of a stage fitted (`develop`), the gains of the same candidates when nothing is missing (`null_gains`), the default choice (`decide`), a chosen candidate put into synthetic records (`on_the_twin`), a term on theory held with its bar carried (`carried`) and read on synthetic records that draw its size (`held_on_the_twin`) |
| `criteria.py` | AIC, AICc, BIC with both counts, Hannan and Quinn, the level of a gain in the family of candidates, and the rule that a term on theory is chosen whatever they say |
| `report.py` | the report of a stage in five parts: results, structure of the residuals, the fitted candidates, the default choice, what the stage does not show |
| `noise.py` | the weights from a law of the noise when the researcher has no noise for each sample (`flat`, `settle`), and whether the residuals agree with the law (`check`, `fitted_law`) |
| `notebook.py` | the steps, the stamps, the refusals |
| `examples/walkthrough.py` | one rung walked on synthetic records with a hidden term: run it first |

The five files from `mechanisms.py` to `noise.py` are new, and sections 2 to 4 of this page do not walk them yet: until they do, a stage ends as section 4 says, and these files are tools that you may call and must report as such.

**What the default choice of a stage does, so that you can say it to the researcher.**

- **A term on theory carries the range that its theory allows.** A term that can only add to the curve is written with its bound, as in `pedestal(8.0, ground="...", lower=0.0)`. Left free, a fit can give it the opposite sign, and it would be added on a ground that forbids that sign. A term that ends at its bound has no gain and is not added on its ground alone.
- **A term on theory comes first.** Give it its ground, `ground="..."`, in the words of the researcher, and the size
  that the theory gives, `theory=(centre, bar)`. If it can be fitted it is chosen whatever its gain. Of several, one
  is added in a stage and the report names the others as the next. If its own parameter cannot be fitted on the
  records, it is added held at the center of the theory: pass the problem and its fit to `decide`, and report the
  systematic bar that each estimate gains beside the bar of the fit. With no size given it is named, not held.

  Before a held term enters, read it on the twin as you would a fitted one:
  `held_on_the_twin(problem, result, choice["held"], table_of=...)`. Each synthetic set holds the term at a size drawn
  from the law of the theory, and is fitted with the term held at the center. The scatter of each estimate over those
  sets holds the noise and the bar of the theory together: where it is longer than the bar of the fit, it is the bar
  to quote. It cannot show that the center of the theory is right.

- **A guess needs the law of the gains.** `null_gains` fits the same candidates on synthetic sets of the working
  model. The guess with the largest gain is chosen when the best candidate gained as much on at most 5 per cent of
  those sets, counted as (count + 1) / (sets + 1). That count cannot be under 0.05 with fewer than 19 sets: ask for 30
  at least, and for 200 or more where a fit is fast, since a level of 0.05 is known to 0.015 from 211 sets.

  With fewer than 39 sets a guess is chosen only when its gain beats every set, so the rule is stricter than 0.05 and
  misses more terms that are there. Counted where nothing is missing, on the worked example with three
  candidates, a term was added on 19 of 320 sets of records, where 1 in 25 to 1 in 31 was expected. Say that rate. With no
  law, or too few sets, no guess is chosen and the stage ends: say so, and say that the researcher may still choose.

- **A law holds for its list.** `decide` refuses a list of candidates other than the one that the law was simulated
  for. Report how many sets were asked, how many ended and how many were dropped.
- **A parameter at its bound** is a term that the records put at zero, or as large as the code lets it be. Report it
  as such, with no bar.
- **The chi square, when the weights come from a law of the noise.** The law has bars, and one bar of it can move
  the mean of the chi square by more than a missing term does. Give `noise.allowance` the law with each of its numbers
  raised by one bar and pass what it returns to the report: it says how far the chi square may sit from its mean for
  no other reason. Inside that, the chi square says nothing of the model: read the structure of the residuals.
- **The table of the residuals by position** gives each sum over its standard deviation after the fit, so 2 is a
  fair limit for one number. With twelve numbers one beyond 2 is no surprise, and a table with none does not show
  that nothing is missing.

This path by hand is the older one. Its table is first order, its fit is least squares, and its twin draws normal
noise of the stated size. The loop of the card reads by windows and by moments, with the modules of the table above,
and reads the shape of the noise too. The other losses of `losses.py` are this path's: the loop fits with none of them.
