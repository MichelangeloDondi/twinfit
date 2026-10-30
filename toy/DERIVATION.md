# Details 1 to 3 from the physics alone
written by a Sonnet 5.5 subagent that was given the experiment and the three details, and did not read the code of the model.

## 1. The natural width

Let Gamma = 1/tau = 2.1944e7 per s be the decay rate of the 6S population, with tau = 45.57 ns. In a weak field the 5P states are far off resonance and drop out, leaving one coupling of Rabi frequency Omega2 between 5S and 6S. The coherence obeys d(rho_eg)/dt = -(gamma_c + i Delta) rho_eg - i Omega2/2, with Delta = 2 pi (nu_0 - 2 nu_L). Its steady state gives the 6S population

`rho_ee = (Omega2^2 / (2 Gamma)) gamma_c / (gamma_c^2 + Delta^2)`.

The ground state is stable, so the coherence decays at half the population rate, gamma_c = Gamma/2, and the line is a Lorentzian of FWHM Gamma in angular frequency:

- On the axis of the atomic transition frequency (2 nu_L): w_N = 1/(2 pi tau) = 3.4925 MHz. The uncertainty is 3.4925 x 0.17/45.57 = 0.013 MHz, so w_N = 3.493 ± 0.013 MHz.
- On the axis of the laser frequency: nu_atom = 2 nu_L halves every width, giving 1.7463 ± 0.0065 MHz.

Time domain: the line is the Fourier transform of the coherence, L(nu) = integral of C(t) exp(2 pi i (nu - nu_0) t) dt, with

`C(t) = exp(-|t| / (2 tau)) = exp(-pi w_N |t|)`,

a decay rate of 10.972 per microsecond.

## 2. Collisions between Rb atoms

A (6S, 5S) pair has no resonant 1/R^3 coupling, because the 5S-6S dipole matrix element vanishes by parity. It feels the van der Waals potential -C6/R^6, and C6 is far larger than for (5S, 5S) because 6S is much more polarizable. During an encounter the transition energy shifts by Delta V(R(t)) and the coherence gains a random phase phi = (1/hbar) integral of Delta V dt. Encounters last about 1e-11 s, far less than the microseconds between them, so they are independent and Poisson distributed (impact limit). Averaging exp(-i sum of phi) over encounters at the rate nu_c = n sigma vbar (density times cross section times mean relative speed) gives exp(-nu_c (1 - <exp(-i phi)>) t). The real part is an extra coherence decay rate gamma_col = n sigma_eff vbar and the imaginary part a small pressure shift. The rate is linear in n because each atom meets perturbers independently and only binary encounters matter.

Density: n = p / (kB T), with p = 101325 Pa x 10^(4.312 - 4040/T).

| T (C) | p (atm) | n (cm^-3) |
|---|---|---|
| 70 | 3.457e-8 | 7.39e11 |
| 90 | 1.539e-7 | 3.11e12 |
| 110 | 5.859e-7 | 1.12e13 |
| 130 | 1.954e-6 | 3.56e13 |

So n rises 48.1-fold from 70 to 130 C. The extra decay multiplies C(t):

`C(t) = exp(-pi (w_N + k n) |t|)`, with k n = gamma_col / pi, so k = sigma_eff vbar / pi.

The line stays one Lorentzian and the widths add. The unit of k is MHz cm^3 (MHz per atom per cm^3). The total width w_N + k n(T) rises with T through n(T) alone and does not depend on power. Order of magnitude: with vbar near 450 m/s and a cross section sigma_eff of 1e-13 to 1e-12 cm^2 (the van der Waals phase reaches one radian at a few nm), k is near 1e-14 MHz cm^3, uncertain by a factor of ten either way, so k n(130 C) lies between 0.04 and 4 MHz.

## 3. Gas that permeated into the cell

Helium broadens the line by the same impact mechanism, with the helium density n_g as perturber density: w_g = n_g sigma_g vbar_g / pi. A closed cell keeps a fixed number of helium atoms over the days of a measurement, because permeation takes years, so n_g is constant: it does not follow the Rb vapour pressure and does not change with T. Only the collision speed and cross section vary: for a van der Waals potential sigma is proportional to v^(-2/5), so sigma v is proportional to T^(3/10), a rise of 5 % from 70 to 130 C. At fixed gas content w_g is a constant of the cell, independent of power and, to about 5 %, of T. It adds a factor exp(-pi w_g |t|) to C(t), so the line stays a Lorentzian and the widths add.

## 4. The line shape and its width

The three decays are independent, so their factors multiply and C(t) = exp(-pi w |t|). The line is one Lorentzian, `L(nu) = (w / (2 pi)) / ((nu - nu_c)^2 + w^2 / 4)`, centred at nu_c, which carries the pressure shifts. Its width is

`w(T) = 1/(2 pi tau) + k n(T) + w_g = 3.4925 MHz + k n(T) + w_g`,

with n(T) from section 2.

## 5. Area and height

Take a signal proportional to the 6S population. The two-photon Rabi frequency Omega2 = Omega_a Omega_b / (2 Delta_int), summed over the 5P intermediate states, is a product of two fields, so it is proportional to the intensity and hence to P. Integrating the weak-field population of section 1 over Delta gives pi Omega2^2 / (2 Gamma) per atom, which does not contain gamma_c: collisions move area from the peak to the shoulders and leave it unchanged. The number of atoms in the beam is proportional to n(T). So

`area ∝ P^2 n(T)`, `height ∝ P^2 n(T) / (w_N + k n(T) + w_g)`.

From 25 to 225 mW the area grows by 9^2 = 81. From 70 to 130 C it grows by 48.1 and the height by slightly less.

## 6. Identifiability

Trace i has width w_i = w_N + k n(T_i) + w_g, with its own height, centre and baseline.

(a) The power sweep alone cannot separate w_N and k: its five traces share one density, so only w_N + k n(130 C) is determined and the two estimates are perfectly anticorrelated. The temperature sweep can, if k n(130 C) exceeds the width error, because four densities spanning a factor 48.1 make the columns 1 and n(T_i) independent: the 70 C trace has 2 % of the 130 C collisional width and fixes w_N, then the 130 C trace fixes k. All eight traces can and do best, since the power traces repeat the top density and cut the variance of k by about 2.4 for equal weights.

(b) With w_N, k and w_g free, w_N and w_g multiply the same column of ones, up to the 5 % factor T^(3/10), so the fit is degenerate: the likelihood is flat along w_N + s, w_g - s, and even keeping the T^(3/10) factor the variance inflation is near 3e3. The sum w_0 = w_N + w_g is determined, and k from the temperature sweep. The lifetime breaks the tie: w_N = 3.493 ± 0.013 MHz gives w_g = w_0 - w_N. Without a model of transit and laser noise, w_g absorbs every constant excess width and only bounds the gas width from above.

(c) With one constant width w_1 and no collisions, w_1 averages the true widths with weights equal to the squared peak height, proportional to (P^2 n)^2 for equal noise. The 130 C traces at 225, 175 and 125 mW carry 92 % of the weight, so w_1 = w_N + w_g + 0.95 k n(130 C): about 3.5 MHz plus a few tenths of MHz for k near 1e-14 MHz cm^3, above the natural width by far more than its 0.4 % uncertainty. For a relative width error eps the misfit is an even residual of eps/2 of the peak height at the line centre, with opposite sign on the shoulders. It is negligible at 130 C and about 5 % of the peak on the 70 and 90 C traces, where the model is too wide, but visible only where the peak signal to noise exceeds a few tens.

## 7. What the model leaves out

In order of expected importance for these traces:

1. Transit-time broadening: atoms cross a beam focused to tens of micrometres in about 0.2 microsecond, comparable to the 6S lifetime, giving a non-Lorentzian width of 1 to 4 MHz that scales as sqrt(T) / waist.
2. The AC Stark shift follows the local intensity, so over the beam profile it gives a shift, skew and broadening that grow linearly with power.
3. Saturation and depletion at the highest powers: when Omega2 nears Gamma the area grows slower than P^2 and the line broadens.
4. Laser frequency noise and scan imperfections (nonlinearity, hysteresis, drift): laser jitter counts twice on the atomic axis.
5. The Doppler pedestal of same-direction photon pairs (about 0.9 GHz wide at 130 C) and neighbouring hyperfine and isotope lines: a slope under the line that a straight baseline only partly absorbs.
6. The detection chain: filter bandwidth and the noise law.
