package simple

import (	
	"github.com/phil-mansfield/num"

	"github.com/phil-mansfield/halo"
	"github.com/phil-mansfield/halo/cosmo"
)

type haloInterface interface {
	MassEnclosed(r float64) float64
	OverdensityRadius(rho float64) float64
	ThermalPressure(ppt PressureProfileType, pt PressurePopulationType, r float64) float64
	DPThermalDr(ppt PressureProfileType, pt PressurePopulationType, r float64) float64
	MinR() float64 
}

type DensityInfo struct {
    M, R, C float64
}

type Halo struct {
    A200, C200, C500 DensityInfo
	Rs               float64
    Z                float64

    FThermal num.Func1D
    BFrac    num.Func1D
    AlphaBias num.Func1D
    BetaBias num.Func1D

	ppt PressureProfileType
}

// Typechecking
var _ haloInterface = new(Halo)

func New(fTh RadialFuncType, ppt PressureProfileType, cType halo.ConcentrationType, m200c, z float64) *Halo {
	h := new(Halo)
	h.Z = z
	h.ppt = ppt

	cFunc := halo.ConcentrationFunc(cType, z)

	h.C200.R = haloRadius(m200c, 200 * cosmo.RhoCritical(z))
	h.C200.M = m200c
	h.C200.C = cFunc(m200c)

	h.Rs = h.C200.C * h.C200.R

	h.C500.R = h.OverdensityRadius(500 * cosmo.RhoCritical(z))
	h.C500.M = haloMass(h.C500.R, 500 * cosmo.RhoCritical(z))
	h.C500.C = h.Rs * h.C500.R

	h.A200.R = h.OverdensityRadius(200 * cosmo.RhoAverage(z))
	h.A200.M = haloMass(h.A200.R, 200 * cosmo.RhoAverage(z))
	h.A200.C = h.Rs * h.A200.R

	return h
}
