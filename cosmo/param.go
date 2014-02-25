package cosmo

import (
	"math"
)

// HubbleFrac calculates h(z) = H(z)/H0. Here H(z) is from Hubble's Law,
// H(z)**2 + k (c/a)**2 = H0**2 h100**2 (OmegaR a**-4 + OmegaM a**-3 + OmegaL).
// The hubble's constant in const.go is H0 = H(z = 0). An alternate
// formulation is h(a) = da/dt / (a H0). Assumes k = 0.
func HubbleFrac(z float64) float64 {
	return math.Sqrt(OmegaR*math.Pow(1.0+z, 4.0) +
		OmegaM*math.Pow(1.0+z, 3.0) + OmegaL)
}

func rhoCriticalMks(z float64) float64 {
	H := HubbleFrac(z) * H0Mks * H100
	return 3.0 * H * H / (8.0 * math.Pi * GMks)
}

// RhoCritical calculates the critical density of the universe. This shows
// up (among other places) in halo definitions and in the definitions of
// the omages (OmegaFoo = pFoo / pCritical).  The returned value is in
// comsological units.
func RhoCritical(z float64) float64 {
	return rhoCriticalMks(z) * MpcMks * MpcMks * MpcMks / MSunMks
}

// RhoAverage calculates the average density of matter in the universe. The
// returned value is in cosmological units.
func RhoAverage(z float64) float64 {
	return RhoCritical(z) * OmegaM * math.Pow(1+z, 3.0)
}

// RhoBaryonAverage calculates the average density of baryons in the
// universe. The returned values is in cosmological units.
func RhoBaryonAverage(z float64) float64 {
	return RhoCritical(z) * OmegaB * math.Pow(1+z, 3.0)
}
