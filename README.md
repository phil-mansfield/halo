Halo
====

This package allows for the calculation of the physical properties of
dark-matter halos.

`Halo` instances are created with the function `halo.New` and then halo
properties are calculated through that halo's methods. The primary
concern of this library is correcting for biases in the x-ray mass
estimates of mass in large galaxy clusters, but it will work in
other contexts as well.

Documentation can be found at
<http://godoc.org/bitbucket.org/phil-mansfield/halo>.
Documentation for a conveient ASCII table interface can be found at
<http://godoc.org/bitbucket.org/phil-mansfield/table>.

*****
Cosmo
-----

The `halo/cosmo` subpackage defines a number of functions which calculate global
cosmological parameters.

***
Num
---

The `halo/num` subpackage is a small numerical computation library. It
currently supports differentiation, several types of one-dimensional
integration, and various searching functions.

*************
Table-Scripts
-------------

The `halo/table-scipts` subpackage is a collection of quick analysis scripts.
These can serve as usage examples. A real programmer wouldn't include these in a
public repository, but I'm a physicist and not a programmer.

**************************************
The Theoretical Basis For This Package
======================================

The following section is a work in progress.

The central purpose of this library is to correct cluster mass measurements
derived from x-ray measurments. X-ray mass measurements are made by assuming
that the clusters are in hydrostatic equilibrium, so the equation

    dP_th/dr + rho_gas * GM/r^2 = 0

is valid. Since `dP_th/dr` and `rho_gas` can be measured directly
cluster mass can be found at an arbitrary radius. However, from numerical
simulations we know that there is a non-insignificant bulk flow in most
clusters, resulting in an effective pressure gradient which is larger than
the empirically measured thermal pressure gradient.

By dividing `M_true` (using `dP_eff/dr`) by `M_bias` (using `fP_thermal/dr`) we
arrive at the following relation:

    M_true(r)/M_bias(r) = b(r) = 1/f_th(r) (1 - beta(r)/alpha(r))

where

    beta(r) = d ln(f_th)/dr

and

    alpha(r) = d ln(P_th)/dr.

Here `f_th = P_th / P_eff`. `f_th` cannot be measured directly and must come
from average relations in hydrodynamic N-body simulations
(e.g. _Battaglia, et al, 2012_ and _Nelson, et al, 2012_ ).

It is straight-forward to insert fitted formulae into beta and alpha to arrive
at a closed-form expression for b(r), but special care must be taken. Any
reasonable fit for `f_th` will depend on at least the mass or the radius of
the halo (for the purposes of this overview, a halo's radius and mass correspond
to its 500c overdensity sphere, but this library is fully capable of working
with alternate mass definitions). Since `f_th` is measured from simulations
the fit will be made with respect to the *true* mass of the halo, but any
any empirical fits for `dP_th/dr` (e.g. the *Planck Pressure Profile*,
or Arnaud et al's *Universal Pressure Profile* ) will be made with respect to
the halo's *biased* mass. For this reason it is strongly suggested that you
use pressure profiles derived from simulations. For maximal internal
consistency the pressure profile used in the calculate of `f_th` should be
used.

It is still possible to compute a halo's bias from empirical pressure profiles,
but if one starts with the true mass of the halo and wishes to compute the
biased mass, one must solve

     0 = R_bias - R_overdenisty(rho_500, R_bias)

numerically. R_overdensity cannot be found trivially because biased halos
do not follow an NFW profile (see below), so the `dP_th/dr` fit and its
dependance on R_bias must be used.

Regardless of what type of pressure profile is assumed, calculating the
bias of a halo based off of its biased mass requires numerically solving

    0 = b - 1/fth(R_bias) (1 - alpha(R_Bias, M_true = b(R_bias)*M_bias) / beta(R_bias, M = b(R_bias)*M_bias*)

for b. This is neccesary because although haloes are traditionally assumed
 to follow NFW profiles, these profiles are measured using the true mass
profiles found in simulations. There is no a priori garuantee that
`M_true(r) / b(r)` can also be fit to the same functional form (in fact,
for the profiles used here it most certainly cannot be). Thus, in order to
 calculate `M_bias`, onne must compute `M_true(r)/b(r)` and to calculate
`rho_bias(r)` one must calculate `(d/dr (M_true(r)/b(r)))/(4 pi r^2)`.

******************
Additional Caveats
------------------

There are three quantities which are very easy to confuse, but are
quantitatively fairly different. These are `b(R_bias)`, `b(R_true)`, and
`M_true/M_bias`. Since the radius of a halo is related to its mass by the
relation `r = (m/(rho pi 4/3))^(1/3)`, a change in halo mass corresponds
to a a change of radius meaning that `b(R_bias)` and `b(R_true)` are
significantly different quantities. Furthermore, `b(R_true)` must not be
confued with `M_true/M_bias`. `M_bias(R_true)` is neccesarily larger than
`M_bias(R_bias)`, meaning that to correct the measured x-ray masses it is
insufficient to simply compute `b(r)`. Instead equivelent haloes must actually
be constructed.