package cosmo

import (
	"math"

	"bitbucket.org/phil-mansfield/halo/num"
)

type SigmaType uint32

const (
	MultiDark2010 SigmaType = iota

	sigmaTypeCount
)

const (
	mdA     = 16.9
	mdB     = 1.102
	mdC     = 6.22
	mdAlpha = 0.41
	mdBeta  = 0.2
	mdGamma = 0.333

	mdPivotMassH = 1e12
)

func pradaX(a float64) float64 {
	return a * math.Pow(OmegaL/OmegaM, 1.0/3.0)
}

// DFluctuaction gives the linear growth rate of density fluctuations at
// a = 1 / (1 + z).
func DFluctuation(a float64) float64 {
	x := pradaX(a)

	term1 := (5.0 / 2.0) * math.Pow(OmegaM/OmegaL, 1.0/3.0)
	term2 := math.Sqrt((1.0 + x*x*x) / (x * x * x))

	dimlessIntegral := num.Integral(func(xp float64) float64 {
		return xp / math.Pow(1.0+xp*xp*xp, 1.5)
	}, 0, x, num.Linear, num.Flat)

	return term1 * term2 * dimlessIntegral(x)
}

// SigmaFunc returns a function which transforms a 200c mass into the
// cosmology-independant mass-proxy, sigma.
func SigmaFunc(sType SigmaType, z float64) num.Func1D {
	switch sType {
	case MultiDark2010:
		a := 1.0 / (1.0 + z)
		D := DFluctuation(a)
		return func(m200c float64) float64 {
			m200cH := m200c * H100
			y := mdPivotMassH / m200cH

			return (D * mdA * math.Pow(y, mdAlpha) /
				(1 + mdB*math.Pow(y, mdBeta) +
					mdC*math.Pow(y, mdGamma)))
		}
	}

	panic("Given unrecognized SigmaType")
}
