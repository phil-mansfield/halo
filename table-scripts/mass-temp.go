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
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)
	fThP = halo.FThermalFunc(simFth, halo.PlusSigmaCurve)
	fThM = halo.FThermalFunc(simFth, halo.MinusSigmaCurve)

	cFunc0 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)
	cFunc5 = halo.ConcentrationFunc(halo.Bhattacharya2013, 0.5)

	colNames = []string {
		"m500c",
		"temp-z=0",
		"temp-z=0+",
		"temp-z=0-",
		"temp-z=0.5",
		"temp-z=0.5+",
		"temp-z=0.5-",
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

		h0 := halo.New(fTh, simPpt, cFunc0, halo.Biased, mass, 0.0)
		h0p := halo.New(fThP, simPpt, cFunc0, halo.Biased, mass, 0.0)
		h0m := halo.New(fThM, simPpt, cFunc0, halo.Biased, mass, 0.0)

		h5 := halo.New(fTh, simPpt, cFunc5, halo.Biased, mass, 0.5)
		h5p := halo.New(fThP, simPpt, cFunc5, halo.Biased, mass, 0.5)
		h5m := halo.New(fThM, simPpt, cFunc5, halo.Biased, mass, 0.5)

		outTable.AddRow(mass,
			h0.EWTemperature(tempBt, valPpt, h5p.C500.R),
			h0p.EWTemperature(tempBt, valPpt, h0p.C500.R),
			h0m.EWTemperature(tempBt, valPpt, h0m.C500.R),
			h5.EWTemperature(tempBt, valPpt, h5p.C500.R),
			h5p.EWTemperature(tempBt, valPpt, h5p.C500.R),
			h5m.EWTemperature(tempBt, valPpt, h5m.C500.R))
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "mass-temp.table"))
}
