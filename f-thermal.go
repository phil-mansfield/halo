package halo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

type FThermalType int

const (
	Battaglia2012 FThermalType = iota
	Battaglia2013
	Nelson2012
)

type FThermalCurveType int

const (
	MeanCurve FThermalCurveType = iota
	PlusSigmaCurve
	MinusSigmaCurve
)

const (
	alphaNelson = 0.14
	alphaPNelson = 0.14 + 0.09
	alphaMNelson = 0.14 - 0.09

	nNelson = 1.16
	nPNelson = 1.16 + 0.17
	nMNelson = 1.16 - 0.17

	a0AmpBattaglia   = 0.18
	betaAmpBattaglia = 0.5
	nrAmpBattaglia   = 0.8
	nmAmpBattaglia   = 0.2

	pivotMassBattaglia = 2e14

	a0Battaglia  = 0.927852
	a0mBattaglia = 0.838161
	a0pBattaglia = 0.966221

	xc0Battaglia  = 1.92058
	xc0mBattaglia = 1.95477
	xc0pBattaglia = -2.38706

	alpha0Battaglia  = 0.345055
	alpha0mBattaglia = 0.396751
	alpha0pBattaglia = 0.386167

	xcNzBattaglia    = -0.44
	xcNmBattaglia    = -0.32
	alphaNzBattaglia = 0.33
	alphaNmBattaglia = -0.18
	aNzBattaglia = -0.14
	aNmBattaglia = -0.027
)

type RadialFuncType func(h *Halo, r float64) float64

// FThermalFunc creates a function which calculates fThermal at some radius
// r.  The return function takes the parameters f(r, r500, m500, z).  The
// calculation is done as if r, r500, and m500 are the true mass and radii
// of the halo.
func FThermalFunc(ftt FThermalType, ftct FThermalCurveType) RadialFuncType {
	switch ftt {
	case Nelson2012:
		switch ftct {
		case MeanCurve:
			return func(h * Halo, r float64) float64 {
				return 1 - alphaNelson * math.Pow(r / h.C500.R, nNelson)
			}
		case PlusSigmaCurve:
			return func(h * Halo, r float64) float64 {
				return 1 - alphaPNelson * math.Pow(r / h.C500.R, nPNelson)
			}
		case MinusSigmaCurve:
			return func(h * Halo, r float64) float64 {
				return 1 - alphaMNelson * math.Pow(r / h.C500.R, nMNelson)
			}
		default:
			panic("Unrecognized ftct")
		}

	case Battaglia2012:
		switch ftct {
		case MeanCurve:
			// This is measured from simulations, so the pivot radius is
			// the corrected r500.
			return func(h *Halo, r float64) float64 {
				x := r / h.C500.R
				return 1.0 - a0AmpBattaglia*
					(math.Pow(1.0+h.Z, betaAmpBattaglia)*
						math.Pow(h.C500.M/pivotMassBattaglia, nmAmpBattaglia)*
						math.Pow(x, nrAmpBattaglia))
			}
		default:
			panic("Only MeanCurve is accepted for Battaglia2012")
		}
	case Battaglia2013:
		makeFTh := func(A0, xc0, alpha0 float64) RadialFuncType {
			return func(h *Halo, r float64) float64 {
				x := r / h.C500.R
				mFrac := h.C500.M / pivotMassBattaglia

				xc := xc0 * math.Pow(1.0+h.Z, xcNzBattaglia) *
					math.Pow(mFrac, xcNmBattaglia)
				alpha := alpha0 * math.Pow(1.0+h.Z, alphaNzBattaglia) *
					math.Pow(mFrac, alphaNmBattaglia)
				A := A0 * math.Pow(1.0+h.Z, aNzBattaglia) *
					math.Pow(mFrac, aNmBattaglia)

				return A * math.Pow(1.0+math.Pow(x/xc, 2.0), -alpha)
			}
		}

		switch ftct {
		case MeanCurve:
			return makeFTh(a0Battaglia, xc0Battaglia, alpha0Battaglia)
		case PlusSigmaCurve:
			return makeFTh(a0pBattaglia, xc0pBattaglia, alpha0pBattaglia)
		case MinusSigmaCurve:
			return makeFTh(a0mBattaglia, xc0mBattaglia, alpha0mBattaglia)
		}

		panic("Unrecognized FThermalCurveType for Battaglia2013")
	}
	panic("Unrecognized FThermalType")
}

func BetaBiasFunc(fTh RadialFuncType) RadialFuncType {
	return func(h *Halo, r float64) float64 {
		fThR := func(r float64) float64 { return fTh(h, r) }
		return num.Derivative(fThR, r)(r) * r / fThR(r)
	}
}

func AlphaBiasFunc(fTh RadialFuncType) RadialFuncType {
	return func(h *Halo, r float64) float64 {
		return h.DPdr(ThermalPressure, h.ppt, ElectronPressure, r) * r /
			h.Pressure(ThermalPressure, h.ppt, ElectronPressure, r) *
			cosmo.MpcMks
	}
}

// BFracFunc returns a function which calculates b(r) for a given halo.
// Here b(r) = mTrue(< r) / mBias(< r).
func BFracFunc(fTh RadialFuncType, corr num.CurveCorrectionType) RadialFuncType {
	return func(h *Halo, r float64) float64 {
		fThR := func(r float64) float64 {
			return fTh(h, r)
		}

		switch corr {
		case num.FirstOrder:
			return 1.0 / fThR(r)
		case num.SecondOrder:
			return (1.0 - h.BetaBias(r)/h.AlphaBias(r)) / fThR(r)
		}
		panic("Recieved unknown CurveCorrectionType")
	}
}
