package halo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

func mNFW(x float64) float64 { return math.Log(1.0+x) - x/(1.0+x) }

// MassEnclosed calculates the mass enclosed within a shell centered on the
// halo's center with radius r.
//
// Implementation note: Calculations are done with m200c and c200c, so
// be sure to calculate 200c DensityInfo before constructing 200a
// DensityInfo and 500c DensityInfo with this function.
func (h *Halo) MassEnclosed(bt BiasType, r float64) float64 {
	switch bt {
	case Biased:
		return h.MassEnclosed(Corrected, r) / h.BFrac(r)
	case Corrected:
		switch h.mpt {
		case NFW:
			x := r / h.Rs 
			return h.C200.M * (mNFW(x) / mNFW(h.C200.C))
		}
		panic("Given unrecognized MassProfileType.")
	}
	panic("Given unrecognized BiasType.")
}

// OverdensityRadius calculates the radius at which h has an average density
// of rho.
func (h *Halo) OverdensityRadius(bt BiasType, rho float64) float64 {
	radiusToRho := func(r float64) float64 {
		m := h.MassEnclosed(bt, r)
		return haloDensity(r, m)
	}

	return num.FindEqualConst(radiusToRho, rho, h.Rs, h.Rs)
}

// Acceleration computes the acceleration due to gravity of point charge
// placed a ditance r from the center of the halo. Returned value is in mks.
func acceleration(h *Halo, bt BiasType, r float64) float64 {
	mEnclosed := h.MassEnclosed(bt, r)
	dist := r * cosmo.MpcMks
	return cosmo.GMks * (mEnclosed * cosmo.MSunMks) / (dist * dist)
}

// Rho calcualtes the density of the halo's dark matter at radius r. The
// returned value is in cosmological units.
func (h *Halo) Rho(bt BiasType, r float64) float64 {
	switch bt {
	case Biased:
		m := func(r float64) float64 { return h.MassEnclosed(bt, r) }
		return num.Derivative(m, r)(r) / (4 * math.Pi * r * r)
	case Corrected:
		rH, mH, c := h.C500.R, h.C500.M, h.C500.C
		ampl := mH * c * c * c / (4.0 * math.Pi * rH * rH * rH * mNFW(c))
		x := r / h.Rs
		return ampl / (x * (1.0 + x) * (1.0 + x))
	}
	panic("Unrecognized BiasType.")
}

// RhoGas returns the gas density at a given distance from the center of
// the halo. The returned value is in cosmological units.
func (h *Halo) RhoGas(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, r float64) float64 {
	rho := -h.DPdr(pbt, ppt, AllPressure, r) / acceleration(h, bt, r)
	return rho / cosmo.MSunMks * math.Pow(cosmo.MpcMks, 3.0)
}

// GasEnclosed returns the mass of all the gas in a halo enclosed within the
// given radius. The returned value is given in MSolar.
func (h *Halo) GasEnclosed(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, r float64) float64 {
	rho := func(r float64) float64 { return h.RhoGas(bt, pbt, ppt, r) }
	return num.Integral(rho, h.MinR(), 0.1, num.Log, num.Spherical)(r)
}
