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
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)
	fThP = halo.FThermalFunc(simFth, halo.PlusSigmaCurve)
	fThM = halo.FThermalFunc(simFth, halo.MinusSigmaCurve)

	cFunc0 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)

	colNames = []string {
		"m500",
		"m500-biased",
		"mgas-corrected-500",
		"mgas-corrected-500th",
		"mgas-biased-500th",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]
	outTable := table.NewOutTable(colNames...)

	minMassLog, maxMassLog := math.Log10(1e13), math.Log10(1e15)
	logWidth := (maxMassLog - minMassLog) / steps

	for massLog := minMassLog; massLog <= maxMassLog; massLog += logWidth {
		mass := math.Pow(10, massLog)
		
		h := halo.New(fTh, simPpt, cFunc0, halo.Corrected, mass, 0.0)

		mGas500Corrected := h.GasEnclosed(halo.Biased,
			halo.EffectivePressure, valPpt, h.C500.R)
		mGas500ThCorrected := h.GasEnclosed(halo.Biased,
			halo.EffectivePressure, valPpt, h.R500cBias)
		mGas500ThNaive := h.GasEnclosed(halo.Biased,
			halo.ThermalPressure, valPpt, h.R500cBias)

		outTable.AddRow(mass, h.M500cBias, mGas500Corrected,
			mGas500ThCorrected, mGas500ThNaive)
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "plot-type-comp.table"))
}
