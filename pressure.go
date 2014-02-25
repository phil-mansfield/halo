package halo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

// PressureProfileType is a flag corresponds to the fit used to derive a
// radial pressure function. When calcualting h.Pressure(), the profile may
// either require that h.M500cBias and h.R500cBias are initialized or that
// h.C500.M and h.C500.R are initialized. The required type can be found
// with the ppt.RequiresBias() method.
type PressureProfileType int

const (
	Planck2012 PressureProfileType = iota
	Arnaud2009
	BattagliaAGN2012
	BattagliaShockHeating2012
)

// PressurePopulationType is a flag that corresponds to the population of
// particles for which the pressure is being calculated.
type PressurePopulationType int

const (
	AllPressure PressurePopulationType = iota
	ElectronPressure
)

// PressureBiasType is a flag that corresponds to whether the calculated
// pressure is the thermal pressure or the effective pressure. Here,
// "effective" pressure means that the pressure is the pressure required
// for the halo to be in hydrostatic equilibrium. "Naive" thermal
// pressure is calculated using the unbiased masses and radii in the place
// of biased masses and radii.
type PressureBiasType int

const (
	EffectivePressure PressureBiasType = iota
	ThermalPressure
	NaiveThermalPressure
)

const (
	kevToPascal float64 = 1.602177e-10
)

const (
	planckPivotM500H float64 = 3e14
	planckA0kev             = 1.65e-3

	planckP0   = 6.41
	planckC500 = 1.81

	planckAlpha = 1.33
	planckBeta  = 4.13
	planckGamma = 0.31
	planckDelta = 0.12

	arnaudPivotM500H = 3e14
	arnaudA0kev = 1.65e-3
	arnaudAP = 0.12

	arnaudP0   = 8.403
	arnaudC500 = 1.177

	arnaudAlpha  = 1.051
	arnaudBeta = 5.4905
	arnaudGamma = 0.3081

	battagliaPMassPivot = 1e14

	battagliaPGamma = -0.3
	battagliaPAlpha = 1.0

	battagliaP0AGN = 7.49
	battagliaPmAGN = 0.226
	battagliaPzAGN = -0.957

	battagliaX0AGN = 0.710
	battagliaXmAGN = -0.0833
	battagliaXzAGN = 0.853

	battagliaB0AGN = 4.19
	battagliaBmAGN = 0.048
	battagliaBzAGN = 0.615

	battagliaP0ShockHeating = 20.7
	battagliaPmShockHeating = -0.074
	battagliaPzShockHeating = -0.743

	battagliaX0ShockHeating = 0.438
	battagliaXmShockHeating = 0.011
	battagliaXzShockHeating = 1.01

	battagliaB0ShockHeating = 3.82
	battagliaBmShockHeating = 0.0375
	battagliaBzShockHeating = 0.535
)

var (
	// This needs to be multiplied by rhoCrit * m500c / r500c.
	pDeltaBattagliaPre = cosmo.GMks * 250 * (cosmo.OmegaB / cosmo.OmegaM) *
		math.Pow(cosmo.MSunMks, 2.0) / math.Pow(cosmo.MpcMks, 4.0)
)

func (ppt PressureProfileType) RequiresBiasedMass() bool {
	switch ppt {
	case Planck2012:
		return true
	case Arnaud2009:
		return true
	case BattagliaAGN2012:
		return false
	case BattagliaShockHeating2012:
		return false
	}
	panic("Given unrecognized PressureProfileType")
}

func pThermal(h *Halo, ppt PressureProfileType, pt PressurePopulationType, r float64) float64 {

	type battagliaParam func(m500c, z float64) float64
	makeBattagliaParam := func(a0, am, az float64) battagliaParam {
		return func(m500c, z float64) float64 {
			return a0 * math.Pow(1 + z, az) *
				math.Pow(m500c/battagliaPMassPivot, am)
		}
	}

	pDeltaBattaglia := pDeltaBattagliaPre * cosmo.RhoCritical(h.Z) *
		h.C500.M / h.C500.R

	switch pt {
	case AllPressure:
		muFrac := cosmo.Mu / cosmo.ElectronMu
		return pThermal(h, ppt, ElectronPressure, r) * muFrac
	case ElectronPressure:
		switch ppt {
		case Planck2012:
			// This is a measured fit, so we need to use the biased mass.
			mFrac := h.M500cBias / (planckPivotM500H / cosmo.H70)
			x := r / h.R500cBias
			y := planckC500 * x

			scaledPressure := planckP0 / (math.Pow(y, planckGamma) *
				math.Pow(1.0+math.Pow(y, planckAlpha),
					(planckBeta - planckGamma)/planckAlpha))

			P500 := (planckA0kev * math.Pow(mFrac, 2.0/3.0 + 0.12) *
				math.Pow(cosmo.HubbleFrac(h.Z), 8.0/3.0) *
				(cosmo.H70 * cosmo.H70))

			PkeV := P500 * scaledPressure
			return PkeV * kevToPascal

		case Arnaud2009:

			mFrac := h.M500cBias / (arnaudPivotM500H / cosmo.H70)
			x := r / h.R500cBias
			y := arnaudC500 * x

			app := 0.1 - (arnaudAP + 0.1) * math.Pow(x/2, 3.0) /
				(1 + math.Pow(x/2, 3.0))

			scaledPressure := + math.Pow(mFrac, arnaudAP + app) *
				math.Pow(cosmo.H70, -1.5) *
				arnaudP0 / (math.Pow(y, arnaudGamma) *
				math.Pow(1 + math.Pow(y, arnaudAlpha),
				(arnaudBeta - arnaudGamma)/arnaudAlpha))

			P500 := (arnaudA0kev * math.Pow(mFrac, 2.0/3.0) *
				math.Pow(cosmo.HubbleFrac(h.Z), 8.0/3.0) *
				(cosmo.H70 * cosmo.H70))

			PkeV := P500 * scaledPressure
			return PkeV * kevToPascal


		case BattagliaAGN2012:
			x := r / h.C500.R

			p0 := makeBattagliaParam(battagliaP0AGN,
				battagliaPmAGN,
				battagliaPzAGN)
			xc := makeBattagliaParam(battagliaX0AGN,
				battagliaXmAGN,
				battagliaXzAGN)
			beta := makeBattagliaParam(battagliaB0AGN,
				battagliaBmAGN,
				battagliaBzAGN)

			xFrac := x / xc(h.C500.M, h.Z)

			muFrac := cosmo.Mu / cosmo.ElectronMu

			return p0(h.C500.M, h.Z) * math.Pow(xFrac, battagliaPGamma) *
				math.Pow(1.0 + math.Pow(xFrac, battagliaPAlpha),
				-beta(h.C500.M, h.Z)) * pDeltaBattaglia / muFrac

		case BattagliaShockHeating2012:
			x := r / h.C500.R

			p0 := makeBattagliaParam(battagliaP0ShockHeating,
				battagliaPmShockHeating,
				battagliaPzShockHeating)
			xc := makeBattagliaParam(battagliaX0ShockHeating,
				battagliaXmShockHeating,
				battagliaXzShockHeating)
			beta := makeBattagliaParam(battagliaB0ShockHeating,
				battagliaBmShockHeating,
				battagliaBzShockHeating)

			xFrac := x / xc(h.C500.M, h.Z)

			muFrac := cosmo.Mu / cosmo.ElectronMu

			return p0(h.C500.M, h.Z) * math.Pow(xFrac, battagliaPGamma) *
				math.Pow(1.0 + math.Pow(xFrac, battagliaPAlpha),
				-beta(h.C500.M, h.Z)) * pDeltaBattaglia / muFrac
		}
		panic("Given unrecognized PressureProfileType.")
	}
	panic("Given unrecognized PressureType.")
}

// P returns either the electron pressure or effective pressure at a
// distance r from the center of a halo.  The type of pressure profile
// can also be specified, although only one has been coded up so far.
//
// The returned pressure is in MKS units.
func (h *Halo) Pressure(pbt PressureBiasType, ppt PressureProfileType, pt PressurePopulationType, r float64) float64 {
	switch pbt {
	case NaiveThermalPressure:
		tmpR := h.R500cBias
		tmpM := h.M500cBias

		h.R500cBias = h.C500.R
		h.M500cBias = h.C500.M

		p := h.Pressure(ThermalPressure, ppt, pt, r)

		h.R500cBias = tmpR
		h.M500cBias = tmpM

		return p
	case ThermalPressure:
		return pThermal(h, ppt, pt, r)
	case EffectivePressure:
		return pThermal(h, ppt, pt, r) / h.FThermal(r)
	}
	panic("Given unrecognized PressureBiasType.")
}

// DPdr calculates the derivative the halo's pressure profile with respect
// to radius. The returned quantity is in MKS units.
func (h *Halo) DPdr(pbt PressureBiasType, ppt PressureProfileType, pt PressurePopulationType, r float64) float64 {
	p := func(r float64) float64 { return h.Pressure(pbt, ppt, pt, r) }
	return num.Derivative(p, r)(r) / cosmo.MpcMks
}

// ThomsonY calculates the spherical Thompson Y out to the specified radius.
func (h *Halo) ThompsonY(pbt PressureBiasType, ppt PressureProfileType, r float64) float64 {
	P := func(r float64) float64 {
		return h.Pressure(pbt, ppt, ElectronPressure, r)
	}

	intTerm :=  num.Integral(P, h.MinR(), h.Rs, num.Log, num.Spherical)(r)
	constTerm := cosmo.SigmaTMks / (cosmo.MeMks * cosmo.CMks * cosmo.CMks)
	return intTerm * constTerm * cosmo.MpcMks
}
