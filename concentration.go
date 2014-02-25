package halo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

type ConcentrationType uint32

const (
	Duffy2011 ConcentrationType = iota
	Prada2011
	Bhattacharya2013

	concentrationTypeCount
)

const (
	bhattacharyaA     = 1.12
	bhattacharyaB     = 0.53
	bhattacharyaC     = 5.9
	bhattacharyaAlpha = 0.3
	bhattacharyaBeta  = 0.54
	bhattacharyaGamma = -0.35

	bhattacharyaPivotMassH = 5e13

	duffyA float64 = 5.71
	duffyB         = -0.084
	duffyC         = -0.47

	duffyPivotMassH = 1e12

	pradaC0    = 3.681
	pradaC1    = 5.033
	pradaAlpha = 6.940
	pradaX0    = 0.424

	pradaSigma0 = 1.047
	pradaSigma1 = 1.646
	pradaBeta   = 7.386
	pradaX1     = 0.526

	pradaA = 2.881
	pradaB = 1.257
	pradaC = 1.022
	pradaD = 0.060
)

func minFunc(c0, c1, mult, x0 float64) num.Func1D {
	return func(x float64) float64 {
		return c0 + (c1-c0)*(math.Atan(mult*(x-x0))/math.Pi+0.5)
	}
}

func bhattacharyaNu(d, m200cH float64) float64 {
	return (1.12 * math.Pow(m200cH/bhattacharyaPivotMassH, 0.3) + 0.53) / d
}

func pradaX(a float64) float64 {
	return a * math.Pow(cosmo.OmegaL/cosmo.OmegaM, 1.0/3.0)
}

// ConcentrationFunc returns a function which transforms a 200c mass into
// a concentration via the specified paper's fit for all halos (not just
// relaxed ones) at the specified redshift.
//
// A halo's concentration, c is equal to r_s / r_200c, where r_s is a
// parameter in the halo's NFW density profile specifying the distance at
// which the slope is -2 on a logarithmic scale.
func ConcentrationFunc(cType ConcentrationType, z float64) num.Func1D {
	switch cType {
	case Duffy2011:
		firstTerm := duffyA * math.Pow(1+z, duffyC)
		return func(m200c float64) float64 {
			m200cH := m200c * cosmo.H100
			return firstTerm * math.Pow(m200cH/duffyPivotMassH, duffyB)
		}

	case Prada2011:
		x := pradaX(1.0 / (1.0 + z))
		sigma := cosmo.SigmaFunc(cosmo.MultiDark2010, z)
		cMin := minFunc(pradaC0, pradaC1, pradaAlpha, pradaX0)
		sigmaMin := minFunc(pradaSigma0, pradaSigma1, pradaBeta, pradaX1)

		return func(m200c float64) float64 {
			B0 := cMin(x) / cMin(1.393)
			B1 := sigmaMin(x) / sigmaMin(1.393)
			sigmaP := B1 * sigma(m200c)
			C := (pradaA * (math.Pow(sigmaP/pradaB, pradaC) + 1.0) *
				math.Exp(pradaD/(sigmaP*sigmaP)))
			return B0 * C
		}

	case Bhattacharya2013:
		d := cosmo.DFluctuation(1.0 / (1.0 + z))
		return func(m200c float64) float64 {
			nu := bhattacharyaNu(d, m200c*cosmo.H100)
			return (math.Pow(d, 0.54) * 5.9 * math.Pow(nu, -0.35))
		}
	}

	panic("Unrecognized ConcentrationType")
}
