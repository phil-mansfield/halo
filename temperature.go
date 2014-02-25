package halo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

const (
	kelvinToKeV       float64 = cosmo.KBMks / (1000.0 * cosmo.EVMks)
	densityCosmoToMks         = cosmo.MSunMks / (cosmo.MpcMks * cosmo.MpcMks * cosmo.MpcMks)
	coolingConst              = -1.0 // This gets canceled out.
)

func coolingLambda(temp float64) float64 {
	return coolingConst * math.Sqrt(temp)
}

// Switch temp calculation from using an arbitrary profile to a thermal profile.
// Get rid of pE factor.

func createNumFunc(h *Halo, pbt PressureBiasType, ppt PressureProfileType, bt BiasType) num.Func1D {
	return func(r float64) float64 {
		density := h.RhoGas(bt, pbt, ppt, r) * densityCosmoToMks
		pE := h.Pressure(pbt, ppt, ElectronPressure, r)
		temp := pE * cosmo.ElectronMu * cosmo.MHyMks * kelvinToKeV /
			(cosmo.KBMks * density)
		return pE * density * density * temp * coolingLambda(temp)
	}
}

func createDenFunc(h *Halo, pbt PressureBiasType, ppt PressureProfileType, bt BiasType) num.Func1D {
	return func(r float64) float64 {
		density := h.RhoGas(bt, pbt, ppt, r) * densityCosmoToMks
		pE := h.Pressure(pbt, ppt, ElectronPressure, r)
		temp := pE * cosmo.ElectronMu * cosmo.MHyMks * kelvinToKeV /
			(cosmo.KBMks * density)
		return pE * density * density * coolingLambda(temp)
	}
}

// EWTemperature computes the emission-wieghted temperature of a halo
// integrated out to the specified radius.  Here, the emission-weighted
// temperature is:
//
// int dr r**2 rhoG**2 * Lambda * T / int dr r**2 * rhoG**2 * Lambda
//
// Here, Lambda is a cooling coefficient taken to be proportional to 
// sqrt(T) and the temperature is calculated from the pressure profile which
// is specified by pbt and ppt. bt is necceary because we're using PV = nRT
// and thus need ot figure out the 
//
// Return value is in keV.
func (h *Halo) EWTemperature(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, rMax float64) float64 {
	rMin := h.MinR()
	scale := math.Log10(rMax) - math.Log10(rMin)

	numFunc := createNumFunc(h, pbt, ppt, bt)
	numInt := num.Integral(numFunc, rMin, scale, num.Log, num.Spherical)

	denFunc := createDenFunc(h, pbt, ppt, bt)
	denInt := num.Integral(denFunc, rMin, scale, num.Log, num.Spherical)

	return numInt(rMax) / denInt(rMax)
}
