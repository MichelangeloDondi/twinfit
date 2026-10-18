"""The version of the package: its number, and for each version the digest of its code and what it changed.

A result names the version that made it. Two results under one number were made by one code: `check` raises when a file
of the package is not the one that its number was given to. So a change of the package is a new number here, with the
digest of the changed code and a line for each thing that it changed.

    python -m twinfit.version        the number, and the digest of the code as it stands

Nothing here is run to make a result: a new number alone moves none.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

__version__ = "0.0.7"

VERSIONS = {
    "0.0.2.dev0": {"digest": "fd5bae3ecf4e99de29c84a71083d6f0dd1463b96536e93887f74abc4af644c27", "changes": ("the package as it was taken from the repository twin-moments at its commit 2a55db6",)},
    "0.0.3": {"digest": "b374c06e2d7ead04d50e8897b7cbddfecf57a8ccb323cab3386a59457ee52008", "changes": (
        "the version in a file of its own, version.py, with the digest of the code of each version and what it changed: check() tells a code that was changed with no new number",
        "a model with no shared parameter is read as any other: the ledger takes a design with no column, a family of statistics has then nothing to estimate and is a reading, and the weighted fit stands for the square cost. Another chi function is refused in plain words",
        "unsuspected.py, what can be said before the missing term has a name. level_for_all and missing: one level for all the readings together, counted on synthetic sets of the fit, and one answer at it",
        "classes: the candidates that the records cannot tell apart, as classes, with the range of the bias of each shared parameter over the members of a class",
        "known_parameters and read: five estimates of each shared parameter, by the fit, with the generic term of the largest gain added, fitted and at first order, and by the combination of the families that is free of that term, or of the leading class",
        "licence: whether first order may be used on the problem at hand, counted on its own synthetic sets beside the fit made again",
    )},
    "0.0.4.dev1": {"digest": "f41aa7be16461c59fe5364e77ecd610028bd679a8c8c3aad2b7a78c713ca9438", "changes": (
        "description.py, the input of the package, his message 140: the description of the data collection in one file, the traces, the conditions of each, which covariates vary, the noise, the physical problem in words, and optional addenda, read with its refusals and made into records",
    )},
    "0.0.4.dev2": {"digest": "0a60c375b3f961a04750a8f84649dd1e68fbdc6ddc83e268f741d0e25ec712b2", "changes": (
        "where.py, where to search a new term, his message 139: the place, the class, the covariate and the family of statistics that the residuals of a fit point to, each with its level and its excess on synthetic sets of the fit, which a count reads once for many records of one design, the classes that the records cannot tell apart, a threshold counted on the twin over which an answer stands beyond it, and the covariates with one value named as what cannot be read",
        "best.py, the output, his message 140: the best model after the iterations by the Bayesian criterion, its terms, the models that the records cannot tell apart from it, what was tried and refused, and five parts of the uncertainty of each parameter: the bar of the fit, the coverage and the bias on the twin, the bar of each held value of theory, and the spread over the models not told apart. A parameter correlated with another beyond 0.999 has no bar",
        "proposals.py, the six proposals of his message 142: three terms for the model, each the shape that the residuals give the class that the report leads, fitted and read on its twin with and without it, and three levers of the experiment, each entering the twin as its kind says, measured with the working model and with the first term held in the model",
        "description.py reads the levers of the experiment in its addenda, each with its kind, its factor or its values, the conditions it starts from and its words, with their refusals",
        "loop.py and `python -m twinfit loop DIR`, the loop without an agent: each model is a rung of the notebook with its twin, its report of where to search and its six proposals, each term file of the researcher is fitted, read on its twin and kept or refused with its reason, and the output of best.py ends it. Run with report=False it reads only what keeps or refuses a term, for a count, and a first model may have no shared parameter",
    )},
    "0.0.4.dev3": {"digest": "5669799cae049582e7d279de8e6a50691d433cbe8144ff4a277457f9909de57a", "changes": (
        "agent/AGENT.md begins with the brief of a fresh chat, his message 140: the description and its addenda as the input, each term of physical origin written as a term file of the loop, and the output of best.py with its five parts and the six proposals",
        "the repairs of review 21: the loop and its output say that the twin takes the stated noise as right, and check the residuals against the height of the line; where.py says its threshold as counted, names a covariate that it cannot read, and gives the first answer of a question that stands beyond the twin; a proposal says whether the report found it; too few sets are refused in words",
    )},
    "0.0.4.dev4": {"digest": "2f6e329b1d881a1efea2508053de20d68787606e130f75886d7575b292a57024", "changes": (
        "loop.py: a value that the description holds in its ansatz reaches the bar of theory of the output, as a covariate that every record gives at its value or as a value that the model file holds in its dictionary HELD. An ansatz that would carry nothing, or that disagrees with its records, is refused in words. The import of numpy that nothing used is gone",
        "best.py: theory_moves computes the bar of theory once, for the output and for the proposals, and puts each held value back at its centre",
        "proposals.py and description.py: a lever of kind characterization narrows the uncertainty of a held value, and every lever is ranked on the error of each known parameter with its bar of theory added in squares, so that a better characterization of the setup and a better instrument stand on one scale",
    )},
    "0.0.4.dev5": {"digest": "451c7dc819c37fea5296483251888eef1260d9d2245d72bfb3b599e1377ac779", "changes": (
        "agent/AGENT.md, the second reading of review 21: a right model passes the line of the twin more often than 5 in 100, so the size of an excess is read and one answer just past the line is weak evidence; the coverage and the bias of the output are not quoted where the check of the noise fails; and the bar of theory comes from the values held in the ansatz of the description",
    )},
    "0.0.4.dev6": {"digest": "3ffeabd9667f37437d71df7fad93b5385227fba4936a70af25593f23ce2d468a", "changes": (
        "loop.py, review 22: the output says when a held value moves no parameter of the best model, beside a failed check of the noise it says that the coverage and the bias on the twin are not to be read, a loop that refused nothing says so, the page of a model gives the bands of the twin where few sets widen them, and the values that a model file holds are found from its module, whatever the model is",
        "description.py, agent/AGENT.md and agent/problem_template.py: an ansatz holds a covariate that every record gives at its value or a key of HELD, and the loop refuses any other, said where the input is documented",
        "proposals.py: a lever of characterization records the value that it narrows, and the bar of theory is said to be that of the design of the records",
    )},
    "0.0.4.dev7": {"digest": "39da86bec27677ac45e1abbb64df420ba8ab345102bd286884b7b9aced62463f", "changes": (
        "fit.py, his messages 156 and 164: a second test of an end. On records much cleaner than their noise the step that the first test asks of a model that lacks a term shrinks with the noise and the arithmetic does not, so a fit that had reached its minimum said that it had not ended. A search that stops short of the first test has ended where a whole step promises a fall of chi square under ARITHMETIC of it, and the fit says so in ended_by_the_arithmetic. At the noise of the records the first test ends a fit as before",
        "refine.py: a refit that the model refuses, as a line that overflows does, is a fit that did not end, in each of the six places where a candidate or a set of a twin is fitted, as the twin of the report already took it. On noiseless traces the loop stopped on its fifth model after 57 minutes",
        "loop.py, rule 10 and his messages 133 and 135: a term file that says hold = True and gives theory = (centre, bar) is held at the centre and not fitted, read first on synthetic records that draw its size from its theory, and the output carries what one bar of it moves each other value. Until now the loop fitted every term free, and on noiseless traces a free coefficient of theory took up the widths of the terms still missing",
    )},
    "0.0.4.dev8": {"digest": "c1f88eb509aab429263e5598327f4433ced1d55dd2ad61bd69c4634c565c3f9d", "changes": (
        "where.py, his messages 142 and 165, the scoreboard at a thousandth of the noise: every statistic of the report grows as the square of the residuals, glitches as the largest square and the two of the law of the noise as standard normal deviates on the twin, so that on records much cleaner than their twin no answer leads by the power of the residuals that it reads. On the clean ripple of the toy those two had stood at 1e12 and 2.5e10 against 2.3e6 for the best shape",
        "where.py: periodic_best scans every period from the span to four steps, half a Fourier step apart, all at once and scored as propose scores a candidate, and the period that it finds enters the periodic class. The five octaves of the width of the dictionary, which the loop keeps, named a ripple of 1.3 widths a second component",
        "candidates.py: the curvature and the cubic term of the level are offered each its own in each record too, as a background that each record carries at its own size. The bench of the toy plants a parabola that grows with the height of each trace, and the one for all records did not name it",
        "proposals.py, his message 142 and the red team of his message 167: where the description holds fewer than three levers, the package adds those that every design has, a quieter setup, more samples, more repeats and a wider span, and a better measured value for each that it holds. Each lever says which of his three ways it is, a better setup, a better characterization of the setup or a better data collection. Where the working model shares no parameter, a lever is measured by the size of the term that the report proposes first: the noiseless run of an agent gave no proposal for the experiment on its first model",
    )},
    "0.0.4.dev9": {"digest": "8f4cef5e09e7062ba2ab9cbf31cdd24efdc8d7ae45272d20e4ce17491dc82e1c", "changes": (
        "README.md at the root and pyproject.toml, his message 169: the package is installed with numpy and scipy alone, and carries every file that is not Python. A Python under 3.11, or one without numpy or scipy, is refused in one sentence",
        "examples/fields.py and `python -m twinfit example FIELD DIR`: worked projects of chromatography and of nuclear magnetic resonance, each drawn from a truth held by name and given in the units of the fit",
        "report.by_position_and_record and the question changes along of where.py, his message 170: the structure of the residuals curve by curve, of each family, and how it changes along each setting, against the twin",
        "`init` writes a folder that `loop` reads: the description and the first model to fill, traces/ and terms/, and the rules and the skill of an agent, which name the loop. AGENT.md gives an iteration of the agentic version in the nine steps of his message 171",
        "mechanisms.second, tail, level and periodic: a term of one parameter for each kind that the report names. The tail takes out the baseline of its wide copy and rises over one sample",
        "the report of each model defines its words before it uses them, gives the first answer of each question, says when none is beyond the twin, names the kind of each proposal for the experiment, and the size of the term where the model shares no parameter. The levers are measured with the first proposal for the model only where it stands beyond the twin",
        "best.md defines the words of its table, with the precision of the coverage, and every value is printed at the decimals of its bar. A refused term says whether the criterion rose or fell too little. The docstring of best.py no longer says that the bias on the twin is taken out of the value, which no version did",
        "`python -m twinfit terms` lists the terms and the form of a term file, and `python -m twinfit scan DIR TERM X...` fits a term whose distance, stretch or period is chosen at each value, with the range of one standard deviation of a profile",
        "tutorial.py and TUTORIAL.md: the tutorial, written from a run of every command in an empty folder, each sentence that reads a result refused when the files of the run do not hold it, and `--check`",
        "the repairs of round 1 of the cold tests, reviews/COLD_TUTORIAL_1.md: a term whose twin ends on too few sets is refused with its reason, and a run that stops on a fault puts the last output back. Each refusal names its fault: a term file, a setting or a parameter of along, a column of a curve, a CENTRE, a missing description, a word of help, a word without quotes",
        "the twin says what it measured, and a refusal by it names the statistic that failed. The check of the noise says how many bands it read and how to read a failure while the report finds something",
        "each proposal for the model gives the line of its term file, a distance in the unit of the abscissa, or says that no term of the package holds it. Kinds that tie are said to tie. The dictionary has a slope of the level. best.md gives the own parameters of each curve with bars marginal over the shared ones, best.own_bars",
        "the chance of an answer, on the page, where the level was a chance and a baseline at once. The arguments of each term are said. The worked chromatograms are at positive times",
        "AGENT.md gives an iteration in the ten steps of the owner's messages 171 and 172: the tenth presents three ways to improve the design of the experiment, each as a lever that the twin measures",
        "the repairs of round 2 of the cold tests, reviews/COLD_TUTORIAL_2.md: the proposals for the model come by the gain that each brings, a term that the model holds is not proposed again, the area of a tail is bounded at zero, every lever measured is said, and a number is taken as the size of widened and the reach of displaced",
        "the page and best.md say what a width is in the names of the proposals, how the two questions along a setting differ, what a search and seen by are, in which units and on which sets each bias is read, a negative area, and a size that holds at a chosen number. The command line says when one answer stands just past its line, a term file named without two digits, and an option before the folder",
    )},
    "0.0.4.dev10": {"digest": "fdbf2365eb0ba79c7b053fa369d25656b8a9089ccda7475815dc1638d9e561d3", "changes": (
        "the repair of round 3 of the cold tests that the counts read, reviews/COLD_TUTORIAL_3.md: second and pedestal take out of their copy the straight line that the model draws far from its curve, as the tail did, so that they add no share of the offset and the slope of a record. The own offset read 8 per cent low under a second component of 0.08, and on curves on a level of 1 the search of the size stopped. level takes its middle and its half span from all the records together, so that one size is one slope on records of any span",
    )},
    "0.0.4.dev11": {"digest": "93f443f4ee04f82dcb747e78b3401cc189a9861702afbb861d6cb783c779a32f", "changes": (
        "the repairs of round 3 of the cold tests that the counts do not read, reviews/COLD_TUTORIAL_3.md: each proposal for the model fits the term of its line beside the model and says its size and bar, or that it ends at its bound and does not stand, and those that stand come first. The page of a model gives the fit in labelled lines with its chi square per degree of freedom read, names the kinds that tie in the row of the class, the parameter in the row along a setting, and what to do when the first kind is one of the noise. Each lever gives its ratio for each known parameter, and within how much it is known where it draws noise of its own",
        "the loop refuses a term that builds on a refused term and goes on, marks the output that a refused run puts back, gives the size and bar of a refused term, names the term files past --n, and says the unit of a level. The scan refuses a term that the model holds and skips a term file that cannot be fitted. A fault that the model raises is a refusal in words",
        "the tutorial defines its words before it uses them, gives every key of the description in a block that its builder reads, says what to quote when the noise fails and why a coverage at 30 sets moves, the unit of a scan and its x, and that a term of theory is held with hold = True. The template of the model starts a peak or a dip. The guide drops its contradictions on the description and on a table that the loop does not give",
    )},
    "0.0.4.dev12": {"digest": "6b09d71a5bf7b6da3a455c4dca67a916ba3c3b98039ae4f43bf30f91322b9c23", "changes": (
        "the repairs of round 4 of the cold tests, reviews/COLD_TUTORIAL_4.md, none in the code that the counts read: best.md says what to quote of each parameter from its five parts, or why not yet, marks a table of fewer than 100 sets as provisional, lists the terms held at their theory, and sets the own values of each curve beside its settings. The loop stops without asking a next file where only families of statistics stand past their line, and the first table of a page follows the flags of its answers",
        "a term file reads the dictionary HELD of the model file by its name, and `scan DIR terms/NN.py X...` scans a term of one's own that defines term_at(x), beside the model of the other term files",
        "the guide holds the loop's protocol alone, and the older manual path is MANUAL.md, `python -m twinfit guide manual`. The tutorial defines the loop and a kernel, says what kept means, how to give a noise that is not known, the two widths of the chromatograms, the forms of a scan, and works a ripple fixed to the abscissa as a term of one's own. The docstrings of mechanisms.py say the frame of periodic and the forms of the size of widened",
    )},
    "0.0.4.dev13": {"digest": "fb7d579986f23841ac9811b313c5e4f5238f02dd3d6999f08258ddd16f4eadf2", "changes": (
        "the repairs of round 5 of the cold tests, reviews/COLD_TUTORIAL_5.md, none in the code that the counts read: `abscissa_offset` in the description, taken from the abscissa of every trace, and every refusal of a run on curves far from zero against their span names it, since the fit steps the place of a curve by a share of its size and the first fit of lines near 587.6 nm did not end",
        "best.md gives a value beside a term along a setting at the mean of that setting, and its growth per unit of the setting with the bar that the covariance gives. A term held at its theory stays in the model that the output chooses, and a bar of theory of the arithmetic alone is zero. The residuals of each model are in loop/residuals_K.csv, the check of the noise says which way a bar is wrong, and a setting that the model reads by another name is named",
        "`scan` gives the range of x from a parabola where the range is under its step, and how far that range moves the other values of the model. The tutorial gives a recipe for a noise not known, a dependence written in model.py, the full form of a Mechanism, and two forms of one effect compared in two folders",
    )},
    "0.0.4.dev14": {"digest": "9106fa39d63f29fcd4ca384297320a3768432be3c4db2201d8ecda50bff60307", "changes": (
        "the repairs of round 6 of the cold tests, reviews/COLD_TUTORIAL_6.md, none in the code that the counts of the scoreboard, of power and of pointing read: held.py, a term built again at each call from the dictionary HELD, so that a number of a term of the package held by an ansatz puts its bar of theory on the other values, which read 0 when the term took the number as its file was read",
        "the scan writes the top of its parabola as the value to write, inside the range that it gives, passes a whole number as one, so that a copy of a kept level is refused, and finds a term file from the folder or from where it runs. A background of each record's own size is proposed with the term of one size for all too. best.md gives the value of a term along a setting at each value of the setting with its bar, says which held value moves a parameter, writes a small bar in plain decimals, and says when a held value is not read",
        "the first line of a report no longer says that the noise agrees where its check fails, each lever gives each parameter its ratio and its spread, a negative size of a term of the package is said, and the refusal of a lever names the ansatz that there are. The tutorial works a distance held by an ansatz, says that the offset holds everywhere, that a ripple of unknown phase is a cosine and a sine, and why the line beyond the twin sits above 1. The guide says what K is, and that the run that is quoted is the last",
    )},
    "0.0.4.dev15": {"digest": "2d065289d798fbc50e2faa7d9b60bbfa2dbf7dd03b6ca88c19f889d9fe35da3f", "changes": (
        "the repairs of round 7 of the cold tests, reviews/COLD_TUTORIAL_7.md, none in the code that the counts of the scoreboard, of power and of pointing read: where glitches lead, the page of a model names the samples beyond 5 in whitened residual and sends them to the noise of a copy of the folder, its proposals for the model say that they are read on residuals that glitches lead, and the loop says no term file holds them in the place of write the next one",
        "no value of best.md is to be quoted while the last report still has an answer beyond its twin, said once above the values, and the kept line says when the noise of the best model fails its check",
        "the twin line of a proposal says what the term would cost if it were real, and not a sign that it is real. A refused term behind kept ones is followed by how to take a kept term out, in a copy",
        "the width of the proposals is said as the code reads it, at half the height above the median. A wrap that does not take three arguments, a term file with no why and a setting that the curves lack are refused in words. Rows of a curve are sorted before the noise of neighbours is read",
        "the tutorial runs nothing and writes nothing with no flag, --write writes it. Its words say how to take a kept term out, what a quiet chi square beside a term found means, and how a term multiplies. The guide gives spikes their lawful way in a copy and five iterations where none is told",
    )},
    "0.0.4.dev16": {"digest": "315d080d33c154129edd095252052cb38259b99e944b5560f0680bede3a7e575", "changes": (
        "the repairs of round 8 of the cold tests, reviews/COLD_TUTORIAL_8.md, none in the code that the counts of the scoreboard, of power and of pointing read: the page of the best model no longer ends in a trace of the code on a parameter with no bar, and a parameter whose bar is not finite is said to be at a bound or tied, not quoted with nan",
        "a proposal of a term that the model holds is left out by the parameter of its line, a term refused before says so, and a background of each record's own size says which own value of the model its sizes follow, read before a term along a setting",
        "a refusal by the twin on fewer than 100 sets says that it can turn, the line of the command on a failing noise names a missing term first, the scan skips only a value that a file holds, refuses an option in words and prints its range to the digits of its half width",
        "best.md sets its notes apart and says which term files it did not read, and a report defines the known parameters and the rule of the twin",
        "the tutorial gives the noise of a counter by the slope of the spread against the fitted value, the noise held in a copy, the second 01.py and HELD, and a term that reads the parameters of the model. The guide says that a newly held value moves the numbers of step 9 and that an ansatz comes from outside the records",
    )},
    "0.0.4.dev17": {"digest": "974025b9dac0454d5ac99f75da84ae1eb35dc8a0921798aebf1d1ed27cdac549", "changes": (
        "the repairs of round 9 of the cold tests, reviews/COLD_TUTORIAL_9.md, none in the code that the counts of the scoreboard, of power and of pointing read: where a noise correlated from sample to sample leads, the page gives the number of samples over which the correlation falls, and the rule to run several first samples",
        "best.md reads each own value of the model against each setting and names one that follows a setting along a line, with its slope and its bar, and the command line says so: a shift of a Raman line with the power, absorbed by the own place of each curve, was left to the reader",
        "a term that duplicates the freedom of an own parameter is refused with that cause, one answer just past its line gives one message on the command line and in best.md, a chi square a hundred times too large names the starts first, and a kept term that falls within two bars of none is named",
        "a size along a setting has its unit on its line to quote, and the value at zero of a setting is given where zero lies near the values taken. The bars of the fit, of theory and over the models are said to add in squares where one is asked",
        "the report prints its line of beyond the twin beside its table, defines the level apart, and names no module of the code. No share of an own value of both signs is read, since it is a place and not a size, and walk refuses a hidden term that it does not know in words",
        "the tutorial says that the first proposal can be wrong and a kept term is not looked at again, how an own value that follows a setting is written, how a correlated noise is thinned, the noise curve by curve, the columns of the residuals and the drift of a setting taken in the order of time",
    )},
    "0.0.4.dev18": {"digest": "14d63926c1a375a0ae09f94af8a97faa457a3c80231f2a941cf9bce218ba6f47", "changes": (
        "the repairs of round 10 of the cold tests, reviews/COLD_TUTORIAL_10.md, none in the code that the counts of the scoreboard, of power and of pointing read: where the asymmetry leads, the page names shapes that no term of the package writes, an exponential tail of the line itself, a dispersive share and a derivative, for a Mechanism of one's own",
        "best.md says when the best model still leaves more than ten of chi square a degree of freedom, lists each kept term with its reason, and says the order in which the curves were taken, or that the description does not say it",
        "the proposals are ranked by what the term of the package gains once fitted, a term file whose curve is not a number is named in the refusal, and a term that replaces the line after kept terms is refused with that cause",
        "a description that holds a key or a table that the package reads nowhere is refused, and it takes the order of the curves under [problem] as order. The line of one answer just past its line names the questions that it counts",
        "the tutorial says that the first model is required, how a centre is kept where the place follows a setting, why the own values of the worked spectra stay own, how a line is replaced and that the replacing term goes first, and to read the noise again where a tail reaches the wings",
    )},
    "0.0.4.dev19": {"digest": "9ceb819feddcc78de95cb2aae3da172b3b4f51eee6242dd4bd72e42e69fa91cb", "changes": (
        "the repairs of round 11 of the cold tests, reviews/COLD_TUTORIAL_11.md, none in the code that the counts of the scoreboard, of power and of pointing read: the page of a model no longer ends on a division by zero where a twin's scatter is zero, and the terms of the best model after a held one are all listed",
        "the offset of the abscissa is named from three spans from zero where the first fit does not end, an own value that steps from curve to curve along no line is named, a chosen number has its line to quote, and a negative gain at a rough number says to scan it",
        "the proposal of the kind that leads the row of what comes first, the scan asks for a finer scan after a parabola of three coarse falls, and a refusal for too few sets names the sets and a term that jumps",
        "a term file whose wrap reads past the parameters of the model is refused in words, and the tutorial says how a centre is kept, a term written smooth, a tail written with erfcx, a proposal within the twin read, and an ansatz from the scan's range carried without narrowing",
    )},
    "0.0.4.dev20": {"digest": "81cd13ebe170211b0687081bcfb4458ff64e3a34044209410ba0923d780c7d4a", "changes": (
        "the repairs of round 12 of the cold tests, reviews/COLD_TUTORIAL_12.md, none in the code that the counts of the scoreboard, of power and of pointing read: the start of the template takes the first sample from the left that reaches four fifths of the largest distance from the baseline, so that every curve of a doublet starts on the same line, and the page of a model whose places differ by more than a width while the fit leaves much says that the starts may sit on different lines",
        "a term kept where the report before it asked for none is marked on the command line and in a section of best.md, the report defines the size of a term, a chosen number may be written as a shared parameter, and the middle of a doublet is the place plus half its distance",
        "the scan says that the top of its parabola is the value to write where the values about it are spaced within its range, refuses a scan with no values and a term with too few numbers in words, and the guide says what to do when the first report puts nothing beyond the twin",
    )},
    "0.0.4.dev21": {"digest": "b86c72cee79336e5fa1968b0990f29026234ed2517dbd6e30acee2b65a969974", "changes": (
        "init writes into the rules and the skill of an agent the Python that ran it, which has numpy and scipy, where the plain python of a machine may lack them, and says when a project was already started",
        "the repairs of round 13 of the cold tests, reviews/COLD_TUTORIAL_13.md: a value at zero of a setting by the line is said to be an extrapolation whose bar carries no form, best.md says that the spread over the models reads one folder, an own value along a curve is told to be written as the physics gives it, an order that is not known is said so, and a misspelt setting is refused with no sentence on values that are not numbers",
        "the guide says that a refusal on 30 sets is read again on 100, that the bands are those stated from 475 sets so that the quoted run is made on 500, how K counts in a copy, and that a proposal for the experiment need not be one of each kind, and the tutorial gives the glitches a noise and says that shared and own are arrays",
        "a refusal by the twin whose failed statistics each miss their band by less than their own standard error is said to be a draw as much as a verdict, and the command asks for 500 sets before the next term file; a model whose starts or bounds do not fit its parameters is refused with its file and its counts; a proposal near a kept term is left out; tied tails and second components point to a component of its own width written in model.py; a run on 100 sets is not told to read again on 100; a kept term of one's own with a number in its name has the note of a chosen number; a noise stated too large is said so; and the loop says each step as it starts",
        "the tutorial gives a component of its own width at its own distance in model.py, the whole file of a shared place, the letters of the tailed curve, the copy command, a bar and the five parts in section 1, and the width of the proposals before their list",
    )},
    "0.0.4.dev22": {"digest": "8d11944fa6a4f2b7771d5213c376aa8b1e7775290a83aa4b461ff04fc30440ae", "changes": (
        "the repairs of round 14 of the cold tests, reviews/COLD_TUTORIAL_14.md: a wider span gives each new sample the noise of the nearest sample of the records, where it gave every sample their median, so that a curve of counts, whose noise grows with its height, was forecast 0.26 of its error where a real window twice as long gave 0.98",
        "a term file reads the names of model.py by their names, the scan says when its term does not read x, an along proposal says that it is a polynomial and not a law, the report defines a second component, a kept term of negative size is flagged on the command line, a start that gives no number or is of the wrong scale is named in the refusal, the listing says the side of displaced, and the tutorial shows its outputs in the order of a terminal, reads curves that are not peaks, and warns of the side of a distance",
        "a periodic proposal takes the period of the residuals, read from a periodogram summed over the records, where it took a multiple of the width of the curves, and the page says whether the ripple stands fixed to the abscissa or about each centre; best.md refuses to quote a coverage outside the band of the rule of the twin at the sets read, one rule with the admission; the tutorial says what loop.before is, and the guide how a copy is numbered, what to write where no source is at hand, and that a named lever replaces the default of its kind",
    )},
    "0.0.4.dev23": {"digest": "4d02f28a55c669e923b986775b8fe6bb58db79cd474f5af98003c940478bef0c", "changes": (
        "the repairs of round 15 of the cold tests, reviews/COLD_TUTORIAL_15.md: a first fit that does not end names a start at an end of the abscissa and a curve too narrow for its unit, the start of the template reads a dip against a line through the ends of its curve, a refusal on thirty sets says in best.md that more sets can turn it, and not read is defined",
        "under a noise that leads the report the levers are measured on the model as it is, with a word that more samples buy less under a correlated noise, the check of the noise says when the residuals correlate, a quiet report with families beyond the twin says that a noise read from the records moves with them, a fitted share of model.py is not said chosen by its name, the default levers say what they are for counts, the size of a term names its denominator, and the command thin writes three copies of a folder that keep one sample in every N",
        "the tutorial gives a tailed curve that holds at every distance, with its side, its share and the order of terms, the thin command and a refusal on thirty sets in section 7, and the guide one name of the first file of iterations of a copy",
    )},
    "0.0.4.dev24": {"digest": "e638ed170f03bdc90f63a0454b221f17a14b8ee322b6d827b00023acc04e94fe", "changes": (
        "the repairs of round 16 of the cold tests, reviews/COLD_TUTORIAL_16.md: a term refused by the Bayesian criterion is fitted in the place of the kept term before it, and a tie is said, held in the spread over the models and named in the command, with the parameters that a tied model lacks",
        "a term under which the curve reads no longer an own parameter that the model read is refused, the criterion of each model is on the best page, a quote said with a caution is no longer said not to be quoted, and the default levers of the kinds not named are always measured",
        "row M205: the description alone starts from the generic model of the package, whose width is read from the curves, a description that reads TO FILL is refused, and every term file is read before the first fit",
        "row M207: no word of the work of the repository in the pages of the package, and the glossary of section 1 of the tutorial defines each word of its list before its use, which python -m twinfit.tutorial --check-words checks",
        "section 11 of the tutorial in five steps, with a curve that replaces the line written out whole, its baseline added again, and the guide's step 0",
        "row M213: a curvature exactly singular, which raised an error on Linux, is inverted by its pseudo-inverse, so the fit ends and gives the parameters that move the curves alike a correlation of 1: the records cannot tell them apart",
    )},
    "0.0.4.dev25": {"digest": "47e76137eae7f1f26d44002bd2c038222f8903193c9add05dab500830df6ea56", "changes": (
        "the repairs of round 17 of the cold tests, reviews/COLD_TUTORIAL_17.md: a parameter along a setting through stacked terms is given at each value and at zero by the whole chain, p (1 + a z (1 + b z)), with the bar of the covariance of all of them, where the first term alone gave a width at 0 K of 0.711 for 1.154",
        "a noise given as a column beside families of statistics read beyond the twin is named on the best page and in the command as a column that may come from the values of the records",
        "the red team of the package, his message 181: one run of the loop at a time in a folder, by a lock, the last whole run kept through a kill, a budget of seconds, --max-seconds, a file of the run's progress, loop/STATUS.json",
        "every line to quote carries its caution, the lines of a growth along a setting too, best.json gives a verdict of each parameter, quotable, and is marked stale beside a refused run, an interrupt is said as one, and the correlation of the parameters is read on those whose deviation is a number",
        "a term refused by the twin behind a kept one is fitted in its place, as one refused by the criterion is, and a noise correlated from sample to sample beyond the twin but not first gets the span of its correlation and the way to thin",
        "the top of a scan is read from a parabola through three neighbouring values whose fits ended and lies between them, a best value beside one whose fit did not end is said as the edge of what the fits reach, init writes no model beside a description that names none, and the guide's first step and its times are said again",
    )},
    "0.0.4.dev26": {"digest": "0d9aaaf55ed12a2cfe2ddc5cbe4010b07effbe486238e87b9a8d3e3349ce6f35", "changes": (
        "the repairs of round 18 of the cold tests: a shared parameter that a step leaves moving no curve, as the spacing of a second line whose share the step took to zero, is held for that step and read again at the next, where the whole fit was refused and said to move no curve at its start, and one that moves no curve where the search ends has no bar",
        "past a chi square of 2 for each degree of freedom a correlation of neighbouring residuals is said to be the misfit of the mean, and no thinning is advised, init writes model_template.py beside a description that names no model, the paragraph of the guide before its first iteration is named as it is, and a negative growth along a setting is said as a fall",
        "a value is not to be quoted where the twin does not admit the fit of its model, as the report of the model says",
        "the budget of seconds of the loop, --max-seconds, in the usage of the command and in the list of the README",
        "python -m twinfit counts DIR writes a copy whose counted curves carry the noise of a counter from the fitted values of the best model, which three testers of round 18 wrote by hand",
        "To quote names the parameters of each curve, where a model shares none and a researcher was left with nothing to quote, and says that the twin reads no coverage and no bias of them",
        "thin refuses a copy of fewer than 10 samples a curve and names the largest N that keeps them, the report does not advise such a thinning, and a copy of thin named wrong names the copies that are there",
        "terms lists each term by the name that a report prints beside the call of a term file that writes it",
        "README.md of the package alone, which the pyproject names and the cold tests hand out",
        "best.md says how many searches of its own synthetic sets ended, of how many, and how far its coverage scatters from that of the report of a model",
        "the setup in the description, his message 184 and row M216: [setup.NAME] gives each instrument and setting with its value, uncertainty and unit where known and how it is known, best.md gives it back, and a lever of kind setup is measured on the twin where the model holds a value of its name and is otherwise named as read by no value of the model",
        "the flanks, his messages 184 and 187: three places of where to search, the core within half a width of the centre of a curve, the flanks from there to two widths and the wings beyond, each band of the report in one of them",
        "the gallery, his message 188: six projects of other fields with their truth, written by python -m twinfit example NAME DIR and printed by python -m twinfit guide examples, and twinfit/EXAMPLES.md drawn from runs of the loop on each, what to do and what not to do beside the truth",
    )},
    "0.0.4.dev27": {"digest": "edc23ab825453107e39b37ca6dd396d94da824930fa8e85fb190923aff311b35", "changes": (
        "the figure of a report larger, and its words at 13 points where they were at 9 and 10, his message 191: \"make the fonts bigger to make the figures and the slides more clear\"",
        "the gallery of nine projects says which six a fresh agent worked in the cold tests and which three were written for it and worked by no agent, where its head said that six were, each worked by an agent",
    )},
    "0.0.5rc1": {"digest": "1a7a624f6b6ad8692dffb462dfa3f502e4ccc02366f3b72aebd8abe1da3ba598", "changes": (
        "the first candidate of release 0.0.5, his message 188: the engine of version dev27, whose code the counts of the release read, with the gallery named in every guide and in the usage of the command without a count of its projects written by hand, where seven places said six of nine, and the tutorial written again",
    )},    "0.0.5rc2": {"digest": "c11850e921f22a96058534a4dddb23faeb154a0a40108beaa50b14fc4985eb8f", "changes": (
        "repairs: the acceptance of release 0.0.5, round 19 of the tests of the tutorial: a first model far from its curves was refused as a fit that did not end, since one search takes 80 steps and the minimum of a resonance pulled by more than a width took 161. The loop's first fit is searched again from where it stopped, up to four times",
        "repairs: best.md names the twin that refused a fit, its report's or its own on its own sets, where it said \"as its report says\" beside a report that said the twin admits the fit; the command names the source of a noise stated too large, where it named neighbouring samples beside a noise given as a value; best.json keeps whether a term was refused by a hair",
        "repairs: the gallery acts as a researcher does: a scan writes the top that its command prints, with finer scans until the top stands inside them, each move is said with its command, and each refusal with its sets and whether it was by a hair. Three of its pages drew their numbers from values that no command gave, and one called refused a term that 100 sets keep",
    )},    "0.0.5rc3": {"digest": "83e3bd837fb86b617a3f25a9d4532d1cc6d0c42fe7d02b9fa6db1b9bbdcc2efb", "changes": (
        "his message 194, which moves the candidate of release 0.0.5 to this version: each project of the gallery is written at its own noise, at three times it and with none drawn, and at two kinds of noise, its own and a second, white, correlated over neighbouring samples or that of a counter, by python -m twinfit example NAME DIR --noise L --kind K, and the gallery runs the first model and what to do at each of the four other noises, with the output's verdict on quoting beside each value. At its own noise each project draws the very numbers that it drew before",
    )},    "0.0.5rc4": {"digest": "1f84a5154d38fd164cee06c08465cf6afeeed6af4955300d9b91cb4a255f6364", "changes": (
        "repairs: the blocking finding of round 20 of the tests of the tutorial, the acceptance of 0.0.5rc3. A term on its twin is read on synthetic sets each fitted with the working model first, and that fit had one search of 80 steps where the first fit of the records has had up to 4 more since 0.0.5rc2: a right term, the resonance of a cantilever moving with the pressure, was refused on 0 of 30 sets and 0 of 100, and each of those sets ends in one search of 400 steps. The twin of a term now gets the steps of the records' searches, from twinfit/loop.py",
        "repairs: the misleading findings of round 20 in the loop. A refusal on none of its sets no longer sends to more sets. A refusal of a first fit names the steps of a search and where the last one stopped, and a start of a place outside its abscissa, as starts left in the recorded frame beside an offset. Best.md gives a place with the offset of the abscissa added back and its unit, and the next step where the twin does not admit a fit. The name of a model parameter with a digit, f0, is no longer read as a chosen number. The cost of a proposal on its twin gives the bias in the unit of the parameter, and says where the records stand far from that twin",
        "repairs: the misleading findings of round 20 in the gallery, the command, the tutorial and the guide. The gallery says what the noise changes in its tables, the rows that share one draw, the counter at three times, the term written in each row and the range of its scan, and the excess of each answer. The command names every example in its refusal and says its options plainly, and the README gives them. The tutorial says that starts and bounds are written less the offset, the endings of the loop, that along cannot move a place near zero and that second takes one distance, the lever of the setup, and that the range of a scan is that of white noise. The guide says what the red team of step 7 may do and where the steps of a chain of copies go",
    )},    "0.0.5rc5": {"digest": "1386dde6991a2096fc52653314b87dacc5a2c17b78ab1a0f9c43db158bb5b0c3", "changes": (
        "repairs: rule 16, nothing of his line in the package, which a whole run of the tests found broken since 0.0.5rc1: a comment of where.py named a term of the toy, the setup of the example satellite named the levels of its gas by a word of his line, and the usage of the command named an example by a word of his line. No count reads a comment or an example, so no count is made stale",
    )},    "0.0.5rc6": {"digest": "3654b692d3bfb90ebf589c383d075b4ac31e3aee97d7f9905af4f7f5908e3cce", "changes": (
        "repairs: the README of the package, which 0.0.5rc4 said gave the options of an example, --noise and --kind, as round 20 asked, where only the README of the repository gave them. No count reads it",
    )},    "0.0.5rc7": {"digest": "45d54e9fe4011413cf1dccc4c00bbf538c779af26816888a84ceb52ad84e7721", "changes": (
        "repairs: the README of the package, which lacked the option --floor of counts that the usage of the command gives: found by a test that holds the README to every command of the usage, his message 196's mechanism for a repair said and not made. No count reads it",
    )},    "0.0.5": {"digest": "12bf2884930573dacb136b8bd09831b420daed25761fc74bc85f9a7e024aa120", "changes": (
        "the release of 0.0.5rc7, row M209: its number without rc, and the tutorial and the README of the package naming PyPI first, `pip install twinfit`, for the group meeting of his message 199. The tutorial prints the number, and so is written again",
        "the README of the package opens on five minutes with the atomic clock of the gallery, from `pip install twinfit` to the light shift written in model.py, for the group of Nicola Poli (his message 199)",
        "the findings of round 21 of the tests of the tutorial, `reviews/COLD_TUTORIAL_21.md`: the clock holds the 2 per cent calibration of its power by an ansatz, so its light shift carries a bar of theory, which the gallery keeps and its page shows beside the model that what to do writes, the README gives that model whole with the times of a run, and how one's own curves and the noise of a fraction of atoms are written, and `guide tutorial` and `guide readme` print those pages with their path, where an unknown topic printed the protocol of an agent.",
        "his messages 200 and 203: the release carries 16 of his recorded traces as plots alone, PNG figures under examples/rubidium/ with the page RUBIDIUM.md that `guide rubidium` prints, written by toy/release_traces.py, and the MIT licence and the author of version 0.0.1. The digest of 0.0.5 was written again before its release each time its files moved: no count, kept fit or reading was made under the earlier ones, ed6730aa, 0c9ca6e3 and 59d60701",
        "the summary of `python -m twinfit loop` built outside its f-string: a line break inside a field of one is the grammar of Python 3.12, and pyproject.toml admits 3.11, on which every command of the package would have stopped. The release refuses such a line (2026-10-17)",
    )},    "0.0.6": {"digest": "f5b04d81fdfad2ee8495a848712d69b24ccad5beb3509748538928c819ce6cc8", "changes": (
        "repairs: review 30 of candidate 30, before the upload: 0.0.5 was committed and never uploaded, and a committed version keeps its digest (rule 16), so its repairs are this number",
        "repairs of review 30, moderate 9, before the upload: `guide rubidium` prints the folder of its figures, as the README and RUBIDIUM.md say it does, RUBIDIUM.md says that the plots print the values fitted on his traces, where it said that no number of them is in the package, and its caption of figure 2 names a change along the power and not its cause",
    )},    "0.0.7": {"digest": "2bc416d2103e99650cf0dad239ced7d91416af35cb6266053bc69b475cf7ef8a", "changes": (
        "the opening of the README, which PyPI shows, written as one page that presents the package, how it works, what it gives and what it cannot do, for a reader who meets it first, his message 232. 0.0.6 was never uploaded, so the release is this number",
        "the floors of the versions derived from the speed of numpy, as the other repository of this work derives them, his messages 233 and 234: numpy 2.5, which convolves about ten times faster than numpy 1.26 on the convolution of the term widened, and so Python 3.12 and scipy 1.16, each run by a clean leg, Python 3.12 with numpy 2.5.0 and scipy 1.16.0 and Python 3.14 with numpy 2.5.2 and scipy 1.18.0. pip refuses a version under a floor. The import of the package keeps its check of Python 3.11, numpy and scipy, since every count of the toy reads the file that holds it",
        "the two figures of the README drawn on a blind run of an agent on eight synthetic curves of a two photon line, his messages 235, 237 and 238: every sample, and the sums of what a fit leaves in pairs of windows cut at six scales, each a number that is noise about zero when nothing is missing, read against the twin of the fit, turn by turn",
    )},
}


def digest(root=None) -> str:
    """The digest of the package as it stands in a folder: each of its files by name and text, but this one."""
    home = Path(root) if root else Path(__file__).resolve().parent
    h = hashlib.sha256()
    for p in sorted(home.rglob("*")):
        name = p.relative_to(home).as_posix()
        if p.is_file() and "__pycache__" not in p.parts and p.suffix in (".py", ".md", ".json", ".toml", ".txt", ".csv") and name != "version.py":
            h.update(name.encode() + b"\0" + p.read_bytes() + b"\0")
    return h.hexdigest()


def check(root=None) -> str:
    """The number of the package, when its code is that of its number. Raises when a file was changed and the number
    was not: the change then needs a new number, its digest and what it changed, in `VERSIONS`."""
    found, held = digest(root), VERSIONS[__version__]["digest"]
    if found != held:
        raise ValueError(f"the code of the package is not that of its version {__version__}: its digest is {found} and the version holds {held}. "
                         "A change of the package is a new number in twinfit/version.py, with this digest and what it changed")
    return __version__


def stamp(root=None, last: str | None = None) -> str:
    """Writes the digest of the package as it stands into the entry of its current number, the one step that makes a
    version, and returns it. Refused where the current number is the one that was stamped last and the code moved since:
    a stamp never moves a version that was handed to a tester or committed, and a change asks a new number. `last`: the
    number that was stamped last. The step lived in the scratchpad of one session, which a reboot empties (the red team
    of his message 176)."""
    import re
    home = Path(root) if root else Path(__file__).resolve().parent
    found = digest(home)
    path = home / "version.py"
    text = path.read_text()
    entry = re.compile(r'("' + re.escape(__version__) + r'": \{"digest": ")([0-9A-Za-z]+)(")')
    held = entry.findall(text)
    if len(held) != 1:
        raise ValueError(f"version.py holds no single entry of {__version__} to stamp")
    if last == __version__ and held[0][1] not in (found, "PENDING"):
        raise ValueError(f"version {__version__} was stamped already with {held[0][1]}, and its code moved: give the change a new number")
    path.write_text(entry.sub(lambda m: m.group(1) + found + m.group(3), text))
    return found


if __name__ == "__main__":
    import sys
    if sys.argv[1:2] == ["--stamp"]:          # python -m twinfit.version --stamp [LAST]
        print(__version__, stamp(last=sys.argv[2] if len(sys.argv) > 2 else None))
    else:
        print(__version__, digest())
