package main

// mass-bias-false.go generates a table comparing various types of
// incorrect or partial bias measurements.

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
	cFunc = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)

	colNames = []string {
		"m500true",
		"q500",
		"fth500",
		"fth500bias",
		"bfrac500",
		"bfrac500bias",
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
		h := halo.New(fTh, simPpt, cFunc, halo.Biased, mass, 0.0)

		outTable.AddRow(mass,
			h.C500.M / h.M500cBias,
			1.0 / h.FThermal(h.C500.R),
			1.0 / h.FThermal(h.R500cBias),
			h.BFrac(h.C500.R),
			h.BFrac(h.R500cBias))
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "mass-false-bias.table"))
}
