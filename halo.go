/*
package halo-physics provides basic (and not so basic) methods
describing the physical state of dark matter halos. */
package halo

import (
	"fmt"
	"math"

	"bitbucket.org/phil-mansfield/halo/cosmo"
	"bitbucket.org/phil-mansfield/halo/num"
)

const (
	MinHaloMass = 1e10 // Halos smaller than MinHaloMass may not be created.
	MaxHaloMass = 1e16 // Halos larger than MaxHaloMass may not be created.
)

type BiasType int

const (
	Biased BiasType = iota
	Corrected
)

type DensityType int

const (
	A200 DensityType = iota
	C200
	C500
)

type MassProfileType int

const (
	NFW MassProfileType = iota
)

type DensityInfo struct {
	M, R, C float64
}

type halo interface {
	DensityInfo(d DensityType) *DensityInfo
	Redshift() float64
	PivotRadius() float64
	MinR() float64
	MaxR() float64

	MassEnclosed(bt BiasType, r float64) float64
	Rho(bt BiasType, r float64) float64
	OverdensityRadius(bt BiasType, rho float64) float64

	GasEnclosed(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, r float64) float64
	RhoGas(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, r float64) float64

	Pressure(pbt PressureBiasType, ppt PressureProfileType, pt PressurePopulationType, r float64) float64
	DPdr(pbt PressureBiasType, ppt PressureProfileType, pt PressurePopulationType, r float64) float64
	ThompsonY(pbt PressureBiasType, ppt PressureProfileType, r float64) float64

	EWTemperature(bt BiasType, pbt PressureBiasType, ppt PressureProfileType, rMax float64) float64
}

type Halo struct {
	A200, C200, C500 DensityInfo
	Z                float64
	Rs               float64

	M500cBias, R500cBias float64

	FThermal num.Func1D
	BFrac    num.Func1D
	AlphaBias num.Func1D
	BetaBias num.Func1D

	ppt PressureProfileType
	mpt MassProfileType
}

// typechecking
var (
	_ halo = new(Halo)
)

func haloRadius(m, density float64) float64 {
	return math.Pow(3.0/(4.0*math.Pi)*m/density, 1.0/3.0)
}

func haloMass(r, density float64) float64 {
	return density * 4.0 * math.Pi / 3.0 * r * r * r
}

func haloDensity(r, m float64) float64 {
	return m / (4.0 * math.Pi / 3.0 * r * r * r)
}

// DensityInfo returns the DensityInfo struct associated with the given
// DensityType.
func (h *Halo) DensityInfo(d DensityType) *DensityInfo {
	switch d {
	case A200:
		return &h.A200
	case C200:
		return &h.C200
	case C500:
		return &h.C500
	}
	panic("DensityInfo given invalid DensityType")
}

// Redshift returns the redshift of h.
func (h *Halo) Redshift() float64 { return h.Z }

// PivotRadius returns Rs for h.
func (h *Halo) PivotRadius() float64 { return h.Rs }

// MinR gives the minimum radius of the halo for the purposes of consistent
// numerical integration. Values less than MinR are not garuanteed to be
func (h *Halo) MinR() float64 { return 0.0001 * h.Rs }

// MaxR
func (h *Halo) MaxR() float64 { return 2 * h.C200.R }

func initSearching(h *Halo, cFunc num.Func1D, m, r float64) {
	rho200c := cosmo.RhoCritical(h.Z) * 200

	m200cToR := func(m200c float64) float64 {
		h.C200.R = haloRadius(m200c, rho200c)
		h.C200.C = cFunc(m200c)
		h.C200.M = m200c

		h.Rs = h.C200.R / h.C200.C

		return h.MassEnclosed(Corrected, r)
	}

	m200c := num.FindEqualConst(m200cToR, m, m, m)

	h.C200.R = haloRadius(m200c, rho200c)
	h.C200.C = cFunc(m200c)
	h.C200.M = m200c

	h.Rs = h.C200.R / h.C200.C
}

// initDensityInfo modifies h so that its various DensityInfo fields
// correspond to a halo with the given m500c. Also updates h.Rs and h.rhoS
func initDensityInfo(h *Halo, cFunc num.Func1D, m, r float64) {
	rho500c := cosmo.RhoCritical(h.Z) * 500
	rho200a := cosmo.RhoAverage(h.Z) * 200

	// This sets c200 for us.
	initSearching(h, cFunc, m, r)

	h.A200.R = h.OverdensityRadius(Corrected, rho200a)
	h.A200.C = h.A200.R / h.Rs
	h.A200.M = haloMass(h.A200.R, rho200a)

	h.C500.R = h.OverdensityRadius(Corrected, rho500c)
	h.C500.C = h.C500.R / h.Rs
	h.C500.M = haloMass(h.C500.R, rho500c)

	// This will already be set if the biased mass is required.
	if !h.ppt.RequiresBiasedMass() {
		h.R500cBias = h.OverdensityRadius(Biased, rho500c)
		h.M500cBias = haloMass(h.R500cBias, rho500c)
	}
}

func initBiasedHalo(h *Halo, cFunc num.Func1D, m500cBias float64) {
	rho500c := 500 * cosmo.RhoCritical(h.Z)

	h.R500cBias = haloRadius(m500cBias, rho500c)
	h.M500cBias = m500cBias

	bFracLhs := func(b float64) float64 { return b }
	bFracRhs := func(b float64) float64 {
		initSearching(h, cFunc, b * h.M500cBias, h.R500cBias)
		h.C500.R = h.OverdensityRadius(Corrected, rho500c)
		h.C500.M = haloMass(h.C500.R, rho500c)

		return h.BFrac(h.R500cBias)
	}

	// h.BFrac evaluated at r500cBias
	b := num.FindEqual(bFracLhs, bFracRhs, 1.0, 2.0)

	mR500cBias := m500cBias * b
	initDensityInfo(h, cFunc, mR500cBias, h.R500cBias)
}

func initCorrectedHalo(h *Halo, cFunc num.Func1D, m500cCorrect float64) {
	rho500c := cosmo.RhoCritical(h.Z) * 500

	h.C500.M = m500cCorrect
	h.C500.R = haloRadius(h.C500.M, rho500c)
	initSearching(h, cFunc, h.C500.M, h.C500.R)

	// initDensityInfo may require that r500cBias and m500cBias be initialized.
	if h.ppt.RequiresBiasedMass() {
		// initial guess
		h.M500cBias = h.C500.M
		h.R500cBias = h.C500.R

		r500cBiasLhs := func(r500cBias float64) float64 { return r500cBias }
		r500cBiasRhs := func(r500cBias float64) float64 {
			h.R500cBias = r500cBias
			h.M500cBias = haloMass(r500cBias, rho500c)
			return h.OverdensityRadius(Biased, rho500c)
		}

		h.R500cBias = num.FindEqual(r500cBiasLhs, r500cBiasRhs,
			h.C500.R/2, h.C500.R)
		h.M500cBias = haloRadius(h.R500cBias, rho500c)
	}

	initDensityInfo(h, cFunc, h.C500.M, h.C500.R)

	h.R500cBias = h.OverdensityRadius(Biased, rho500c)
	h.M500cBias = haloMass(h.R500cBias, rho500c)
}

// New creates a new Halo instance using the given parameters. If the given
// mass is the m500c of a biased halo, bt should be set to Biased. If the
// mass is the true m500c of the halo.
func New(fTh RadialFuncType, ppt PressureProfileType, cFunc num.Func1D, bt BiasType, m500c, z float64) *Halo {
	if m500c < MinHaloMass || m500c > MaxHaloMass {
		panic(fmt.Sprintf("The %.5g is not within halo bounds [%.5g, %.5g]",
			m500c, MinHaloMass, MaxHaloMass))
	}

	h := new(Halo)
	h.Z = z
	h.ppt = ppt
	h.mpt = NFW

	bFrac := BFracFunc(fTh, num.SecondOrder)
	alphaBias := AlphaBiasFunc(fTh)
	betaBias := BetaBiasFunc(fTh)

	h.FThermal = func(r float64) float64 { return fTh(h, r) }
	h.BFrac = func(r float64) float64 { return bFrac(h, r) }
	h.AlphaBias = func(r float64) float64 { return alphaBias(h, r) }
	h.BetaBias = func(r float64) float64 { return betaBias(h, r) }

	switch bt {
	case Biased:
		initBiasedHalo(h, cFunc, m500c)
	case Corrected:
		initCorrectedHalo(h, cFunc, m500c)
	default:
		panic("Unrecognized BiasType")
	}

	return h
}
