package main

import (
	"os"
	"path"
	"math"

    "github.com/phil-mansfield/halo"
    "github.com/phil-mansfield/halo/simple"

    "github.com/phil-mansfield/table"
)

const (
	steps = 200
	simFth = simple.Battaglia2013
	simPpt = simple.BattagliaAGN2012
	valPpt = simple.Planck2012
	cType = halo.Bhattacharya2013
)

var (
	fTh = simple.FThermalFunc(simFth, simple.MeanCurve)
	fThP = simple.FThermalFunc(simFth, simple.PlusSigmaCurve)
	fThM = simple.FThermalFunc(simFth, simple.MinusSigmaCurve)

	colNames = []string {
		"m500",
		"mgas500-uncorrected",
		"mgas500-flat-correction",
		"fgas500-uncorrected",
		"fgas500-flat-correction",
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
		
		h := simple.New(fTh, simPpt, cType, mass, 0.0)
		mGasUnc := h.GasEnclosed(simple.Uncorrected, simPpt, h.C500.R)
		mGasFC := h.GasEnclosed(simple.FlatCorrection, simPpt, h.C500.R)
		fGasUnc := mGasUnc / h.C500.M
		fGasFC := mGasFC / h.C500.M
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "mass-fgas-simple.table"))
}
