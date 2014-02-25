package main

import (
	"os"
	"path"
	"math"

	"bitbucket.org/phil-mansfield/halo"
	"bitbucket.org/phil-mansfield/table"
)

const (
	z = 0
	steps = 300
)

var (
	fTh = halo.FThermalFunc(halo.Battaglia2013, halo.MeanCurve)
	cFunc = halo.ConcentrationFunc(halo.Bhattacharya2013, z)

	colNames = []string {
		"r/r500c",
		"P-Battaglia-m=1e14-z=0.0",
		"P-Battaglia-m=1e15-z=0.0",
		"P-Battaglia-m=1e14-z=0.2",
		"P-Battaglia-m=1e15-z=0.2",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]
	outFile := path.Join(outDir, "radial-pressure.table")
	t := table.NewOutTable(colNames...)

	h014 := halo.New(fTh, halo.BattagliaAGN2012,
		cFunc, halo.Corrected, 1e14, 0)
	h015 := halo.New(fTh, halo.BattagliaAGN2012,
		cFunc, halo.Corrected, 5e14, 0)
	h214 := halo.New(fTh, halo.BattagliaAGN2012,
		cFunc, halo.Corrected, 1e14, 0.2)
	h215 := halo.New(fTh, halo.BattagliaAGN2012,
		cFunc, halo.Corrected, 5e14, 0.2)

	minFracLog, maxFracLog := math.Log10(0.01), math.Log10(10)
	logWidth := (maxFracLog - minFracLog) / steps
	for fracLog := minFracLog; fracLog <= maxFracLog; fracLog += logWidth {
		frac := math.Pow(10, fracLog)

		r014 := h014.C500.R * frac
		r015 := h015.C500.R * frac
		r214 := h214.C500.R * frac
		r215 := h215.C500.R * frac

		t.AddRow(frac,
			h014.Pressure(halo.ThermalPressure, halo.BattagliaAGN2012,
			halo.ElectronPressure, r014),
			h015.Pressure(halo.ThermalPressure, halo.BattagliaAGN2012,
			halo.ElectronPressure, r015),
			h214.Pressure(halo.ThermalPressure, halo.BattagliaAGN2012,
			halo.ElectronPressure, r214),
			h215.Pressure(halo.ThermalPressure, halo.BattagliaAGN2012,
			halo.ElectronPressure, r215))
	}

	t.Write(table.KeepHeader, outFile)
}