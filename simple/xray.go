package simple

import (
	"math"

	"github.com/phil-mansfield/halo/cosmo"
	"github.com/phil-mansfield/num"
)

type XRayCorrectionType uint8

const (
	Uncorrected XRayCorrectionType = iota
	FlatCorrection

	FlatCorrectionSize = 1.2
)

const (
	kelvinToKeV       float64 = cosmo.KBMks / (1000.0 * cosmo.EVMks)
	densityCosmoToMks         = cosmo.MSunMks / (cosmo.MpcMks * cosmo.MpcMks * cosmo.MpcMks)
	coolingConst              = -1.0 // This gets canceled out.
)

func acceleration(h *Halo, r float64) float64 {
	mEnclosed := h.MassEnclosed(r)
	dist := r * cosmo.MpcMks
	return cosmo.GMks * (mEnclosed * cosmo.MSunMks) / (dist * dist)
}

func (h *Halo) RhoGas(xct XRayCorrectionType, ppt PressureProfileType, r float64) float64 {
	var rho float64

	switch xct {
	case Uncorrected:
		rho = -h.DPThermalDr(ppt, AllPressure, r) / acceleration(h, r)
	case FlatCorrection:
		rho = -h.DPThermalDr(ppt, AllPressure, r) / acceleration(h, r) * 
			FlatCorrectionSize
	default:
		panic("Unrecognized XRayCorrectionType")
	}
	return rho / cosmo.MSunMks * math.Pow(cosmo.MpcMks, 3.0)
}

func (h *Halo) GasEnclosed(xct XRayCorrectionType, ppt PressureProfileType, r float64) float64 {
	rho := func(r float64) float64 { return h.RhoGas(xct, ppt, r) }
	return num.Integral(rho, h.MinR(), 0.1, num.Log, num.Spherical)(r)
}

func coolingLambda(temp float64) float64 {
	return coolingConst * math.Sqrt(temp)
}

// Switch temp calculation from using an arbitrary profile to a thermal profile.
// Get rid of pE factor.

func createNumFunc(h *Halo, xct XRayCorrectionType, ppt PressureProfileType) num.Func1D {
	return func(r float64) float64 {
		density := h.RhoGas(xct, ppt, r) * densityCosmoToMks

		var pE float64
		switch xct {
		case Uncorrected:
			pE = h.ThermalPressure(ppt, ElectronPressure, r)
		case FlatCorrection:
			pE = h.ThermalPressure(ppt, ElectronPressure, r) * FlatCorrectionSize
		}

		temp := pE * cosmo.ElectronMu * cosmo.MHyMks * kelvinToKeV /
			(cosmo.KBMks * density)
		return pE * density * density * temp * coolingLambda(temp)
	}
}

func createDenFunc(h *Halo, xct XRayCorrectionType, ppt PressureProfileType) num.Func1D {
	return func(r float64) float64 {
		density := h.RhoGas(xct, ppt, r) * densityCosmoToMks

		var pE float64
		switch xct {
		case Uncorrected:
			pE = h.ThermalPressure(ppt, ElectronPressure, r)
		case FlatCorrection:
			pE = h.ThermalPressure(ppt, ElectronPressure, r) * FlatCorrectionSize
		}

		temp := pE * cosmo.ElectronMu * cosmo.MHyMks * kelvinToKeV /
			(cosmo.KBMks * density)
		return pE * density * density * coolingLambda(temp)
	}
}

func (h *Halo) EWTemperature(xct XRayCorrectionType, ppt PressureProfileType, rMax float64) float64 {
	rMin := h.MinR()
	scale := math.Log10(rMin) - math.Log10(rMin)

	numFunc := createNumFunc(h, xct, ppt)
	numInt := num.Integral(numFunc, rMin, scale, num.Log, num.Spherical)

	denFunc := createDenFunc(h, xct, ppt)
	denInt := num.Integral(denFunc, rMin, scale, num.Log, num.Spherical)

	return numInt(rMax) / denInt(rMax)
}
