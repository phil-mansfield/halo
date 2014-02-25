package main

import (
	"os"
	"path"
	"math"

	"bitbucket.org/phil-mansfield/halo"
	"bitbucket.org/phil-mansfield/table"
)

const (
	steps = 200
	simFth = halo.Battaglia2013
	simPpt = halo.BattagliaAGN2012
	valPpt = halo.Planck2012

	tempBt = halo.Corrected

	f0 = 0.0764
	alpha = 0.037
	m0 = 1e15
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)
	fThP = halo.FThermalFunc(simFth, halo.PlusSigmaCurve)
	fThM = halo.FThermalFunc(simFth, halo.MinusSigmaCurve)

	cFunc0 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)

	colNames = []string {
		"m500",
		"m500-biased",
		"fgas500c-corrected-z=0",
		"fgas500c-naive-z=0",
		"mgas500-corrected",
		"mgas500-naive",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]
	outTable := table.NewOutTable(colNames...)

	fGas := func(h *halo.Halo, bt halo.BiasType, pbt halo.PressureBiasType) float64 {
		if bt != halo.Biased {
			return h.GasEnclosed(bt, pbt, valPpt, h.C500.R) / 
				h.MassEnclosed(bt, h.C500.R)
		} else {
			return h.GasEnclosed(bt, pbt, valPpt, h.R500cBias) / 
				h.MassEnclosed(bt, h.R500cBias)
		}	
	}

	minMassLog, maxMassLog := math.Log10(1e13), math.Log10(1e15)
	logWidth := (maxMassLog - minMassLog) / steps

	for massLog := minMassLog; massLog <= maxMassLog; massLog += logWidth {
		mass := math.Pow(10, massLog)
		
		h := halo.New(fTh, simPpt, cFunc0, halo.Corrected, mass, 0.0)

		fGasCorrected := fGas(h, halo.Corrected, halo.EffectivePressure)
		fGasNaive := fGas(h, halo.Corrected, halo.NaiveThermalPressure)
		outTable.AddRow(mass, h.M500cBias, fGasCorrected, fGasNaive,
			fGasCorrected * mass, fGasNaive * h.M500cBias)
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "mass-fgas.table"))
}
