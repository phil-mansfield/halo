package main

// mass-bias.go constructs a table computing Q_500 for a range of unbiased
// masses at different redshifts and different std-dev bins.

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
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)
	fThP = halo.FThermalFunc(simFth, halo.PlusSigmaCurve)
	fThM = halo.FThermalFunc(simFth, halo.MinusSigmaCurve)

	cFunc0 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)
	cFunc2 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.2)

	colNames = []string {
		"m500c-bias",
		"m500c-true/m500c-bias-z=0",
		"m500c-true/m500c-bias-z=0+",
		"m500c-true/m500c-bias-z=0-",
		"m500c-true/m500c-bias-z=0.2",
		"m500c-true/m500c-bias-z=0.2+",
		"m500c-true/m500c-bias-z=0.2-",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("Must provite a target directory.")
	}

	outDir := os.Args[1]
	outTable := table.NewOutTable(colNames...)

	minMassLog, maxMassLog := math.Log10(1e13), math.Log10(1e15)
	logWidth := (maxMassLog - minMassLog) / steps

	for massLog := minMassLog; massLog <= maxMassLog; massLog += logWidth {
		mass := math.Pow(10, massLog)

		h0 := halo.New(fTh, simPpt, cFunc0, halo.Biased, mass, 0.0)
		h0p := halo.New(fThP, simPpt, cFunc0, halo.Biased, mass, 0.0)
		h0m := halo.New(fThM, simPpt, cFunc0, halo.Biased, mass, 0.0)

		h2 := halo.New(fTh, simPpt, cFunc2, halo.Biased, mass, 0.2)
		h2p := halo.New(fThP, simPpt, cFunc2, halo.Biased, mass, 0.2)
		h2m := halo.New(fThM, simPpt, cFunc2, halo.Biased, mass, 0.2)

		outTable.AddRow(mass,
			h0.C500.M / h0.M500cBias,
			h0p.C500.M / h0p.M500cBias,
			h0m.C500.M / h0m.M500cBias,
			h2.C500.M / h2.M500cBias,
			h2p.C500.M / h2p.M500cBias,
			h2m.C500.M / h2m.M500cBias)
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "mass-bias.table"))
}
