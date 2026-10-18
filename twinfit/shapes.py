"""A curve from its characteristic function: a density that is a convolution or a mixture and has no formula.

The model of a recorded curve is often a density that is easier to write by its characteristic function than by itself.
The characteristic function of a convolution is a product, and that of a mixture is a weighted mean, so a curve that
mixes over a latent variable the convolution of several kernels, with the kernels depending on that variable, costs a
few lines. This module turns such a function into the density on the abscissa of a record.

    phi(t) = integral of f(u) exp(2 pi i u t) du             the convention: t counts cycles for each unit of u
    f(u)   = 2 Re integral over t > 0 of phi(t) exp(-2 pi i u t) dt

The integral is taken by Gauss and Legendre quadrature on nodes of t that do not move with the parameters, so the
density is a smooth function of whatever phi depends on, and a finite difference of it is a derivative. The nodes are
laid in panels that begin at t = 0: a factor with a kink there, exp(-pi w |t|), is smooth inside every panel.

    tr = Transform(x, reach)                 x the abscissa. Beyond t = reach the function phi must be negligible
    tr.t                                     the nodes, where the caller gives the values of phi
    tr.density(values, shift=c)              the density at x - c
    reach_of(fwhm, tol)                      a reach for any phi that holds a Cauchy factor of at least that full width
    cauchy(t, fwhm), normal(t, sd), laplace(t, scale), ramp(t, depth)        the functions of four kernels
    voigt(x, fwhm, sd)                       the convolution of the first two by its formula, to check a Transform against

What it costs: one matrix of the size of x by the number of nodes for each abscissa, about 400 nodes for a curve of
width 4 read to 20 each side, and then one product for each density. Nothing here fits or draws anything.
"""
from __future__ import annotations

import numpy as np
from scipy.special import voigt_profile

ORDER = 16                # nodes a panel
CYCLES = 2.5              # at most this many cycles of exp(-2 pi i u t) a panel, at the largest |u|: the rule is then exact to 1e-15


def reach_of(fwhm: float, tol: float = 1e-12) -> float:
    """The t beyond which a Cauchy factor of full width `fwhm` is under `tol`: a reach for any product that holds it."""
    if fwhm <= 0:
        raise ValueError("a reach needs a Cauchy factor of positive width: without one, give the reach yourself")
    return float(-np.log(tol) / (np.pi * fwhm))


class Transform:
    """The density on a fixed abscissa from the values of a characteristic function at fixed nodes."""

    def __init__(self, x, reach: float, margin: float = 3.0, order: int = ORDER):
        x = np.asarray(x, dtype=float)
        if x.ndim != 1 or reach <= 0:
            raise ValueError("an abscissa of one dimension and a positive reach")
        span = float(np.max(np.abs(x))) + margin                 # how far from zero a density is asked, with room for a shift
        panels = max(4, int(np.ceil(span * reach / CYCLES)))
        nodes, weights = np.polynomial.legendre.leggauss(order)
        edges = np.linspace(0.0, reach, panels + 1)
        half = 0.5 * np.diff(edges)
        self.t = (0.5 * (edges[:-1] + edges[1:])[:, None] + half[:, None] * nodes[None, :]).ravel()
        self.w = (half[:, None] * weights[None, :]).ravel()
        self.x, self.reach, self.margin = x, float(reach), float(margin)
        self._kernel = 2.0 * self.w[None, :] * np.exp(-2j * np.pi * x[:, None] * self.t[None, :])

    def density(self, values, shift: float = 0.0) -> np.ndarray:
        """The density at x - shift, for the values of the characteristic function at the nodes `t`."""
        values = np.asarray(values, dtype=complex)
        if values.shape != self.t.shape:
            raise ValueError(f"one value for each of the {self.t.size} nodes")
        if abs(shift) > self.margin:
            raise ValueError(f"a shift of {shift:g} is beyond the margin of {self.margin:g} that the nodes were laid for")
        if shift:
            values = values * np.exp(2j * np.pi * shift * self.t)
        return (self._kernel @ values).real


def cauchy(t, fwhm: float):
    """The Cauchy density of full width at half maximum `fwhm`."""
    return np.exp(-np.pi * fwhm * np.abs(t))


def normal(t, sd: float):
    """The normal density of standard deviation `sd`."""
    return np.exp(-2.0 * (np.pi * sd * np.asarray(t)) ** 2)


def laplace(t, scale: float):
    """The two sided exponential density exp(-|u| / scale) / (2 scale): its variance is 2 scale^2."""
    return 1.0 / (1.0 + (2.0 * np.pi * scale * np.asarray(t)) ** 2)


def ramp(t, depth: float):
    """A displacement depth * v, with v between 0 and 1 of density 2 v: mean 2 depth / 3, variance depth^2 / 18."""
    th = 2.0 * np.pi * depth * np.asarray(t, dtype=float)
    out = np.empty(th.shape, dtype=complex)
    small = np.abs(th) < 1e-2
    s = th[small]
    out[small] = 1.0 + 2j * s / 3.0 - s * s / 4.0 - 1j * s ** 3 / 15.0 + s ** 4 / 72.0     # the series: the formula below loses digits near zero
    b = th[~small]
    out[~small] = 2.0 * (np.exp(1j * b) * (1.0 - 1j * b) - 1.0) / (b * b)
    return out


def voigt(x, fwhm: float, sd: float):
    """The convolution of a Cauchy density of full width `fwhm` and a normal density of standard deviation `sd`."""
    return voigt_profile(np.asarray(x, dtype=float), sd, 0.5 * fwhm)
