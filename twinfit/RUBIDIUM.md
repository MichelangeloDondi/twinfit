# A real case: the two photon line of rubidium

These figures show recorded traces of the two photon line of rubidium from 5S to 6S, by the author of the package, as plots alone: their samples are not in the package, only pictures of them and, printed on the figures, the values that the fits read from them. There are 16: the eight conditions of one peak, four temperatures at one power and five powers at one temperature, then at one of them the three other peaks, four more repetitions and a ruler.

![The eight conditions under a model of five terms](examples/rubidium/1_eight_conditions.png)

*Figure 1. The eight traces of the peak at 993.4192 nm, fitted together under a model of five terms: the natural width, collisions between the atoms, a gas that permeated the cell, the noise of the laser, and the passage of the atoms through the beam. Below each trace, what the fit leaves over its noise.*

![What the package reads on them](examples/rubidium/2_where_to_search.png)

*Figure 2. The report of where to search on those traces, beside one synthetic draw of the same model. On the recorded traces its first answer that points to a term is a change of the core along the power, at an excess of 3.59, where a synthetic set passes 2.84 5 times in 100. Above it in the figure stand a family of statistics that reads the misfit by the square root of the value, at an excess of 10.7, and a noise skewed to one side at 4.81: they say how the misfit is read and what the noise is like, and point to no term.*

That answer stands just past its line, 1.26 times it, on 30 sets, as a model that lacks nothing shows on some of its synthetic sets: the package reads it again on 100 sets before a term is written. Where it holds, something changes with the power that the five terms do not hold. The package names the change and not its cause: in this line the light shift of the probe, a broadening by the power and saturation each change with it, and a term of each, tried by the loop against its twin, tells them apart.

![The other peaks, the repetitions and a ruler](examples/rubidium/3_peaks_repetitions_ruler.png)

*Figure 3. At 130 C and 225 mW, the four peaks over their height, five repetitions of one peak as the oscilloscope wrote them, and a ruler, whose sidebands of a known frequency turn the time of a sweep into a frequency.*

The twin of those traces, a model of the five terms whose values follow their joint fit, three slides of the author's report on it:

![The height and the laser noise of each trace](examples/rubidium/4_twin_each_trace.png)

*Figure 4. The height and the laser noise with which the twin draws each trace.*

![What the fit of five terms leaves](examples/rubidium/5_twin_what_the_fit_leaves.png)

*Figure 5. What the fit of the five terms leaves on the recorded traces.*

![What the report reads, beside the twin](examples/rubidium/6_twin_what_the_report_reads.png)

*Figure 6. The statistics of the report that pass the line of the twin, the line that a right model passes once in 20.*

`python -m twinfit example rubidium DIR` writes this twin as a project: its eight synthetic traces, their description and a first model of one Lorentzian line, so that the loop runs on it as on the recorded traces. The values that drew them are in `examples/rubidium_twin/truth.json`, outside the project.

These are readings of the package on one set of recorded traces, and not measurements. Figures 1 and 2 are pages of the folder of the author's toy model, whose foot names the files that made them. `python -m twinfit guide rubidium` prints this page with the folder of its figures.
