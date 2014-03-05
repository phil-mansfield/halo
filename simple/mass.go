package simple

import (
	"math"

	"github.com/phil-mansfield/halo/cosmo"

	"github.com/phil-mansfield/num"
)

const (
    MinHaloMass = 1e10
    MaxHaloMass = 1e16
)

// Utility functions

func haloRadius(m, density float64) float64 {
    return math.Pow(3.0/(4.0*math.Pi)*m/density, 1.0/3.0)
}

func haloMass(r, density float64) float64 {
    return density * 4.0 * math.Pi / 3.0 * r * r * r
}

func haloDensity(r, m float64) float64 {
    return m / (4.0 * math.Pi / 3.0 * r * r * r)
}

func mNFW(x float64) float64 { return math.Log(1.0+x) - x/(1.0+x) }
 

// Public functions

func (h *Halo) MinR() float64 {
	return haloRadius(MinHaloMass, cosmo.RhoCritical(h.Z) * 500) / 100
}
 
func (h *Halo) MassEnclosed(r float64) float64 {
	x := r / h.Rs
	return h.C200.M * (mNFW(x) / mNFW(h.C200.C))
}


func (h *Halo) OverdensityRadius(rho float64) float64 {
    radiusToRho := func(r float64) float64 {
        m := h.MassEnclosed(r)
        return haloDensity(r, m)
    }

    return num.FindEqualConst(radiusToRho, rho, h.Rs, h.Rs)
}
