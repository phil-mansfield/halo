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

	bt = halo.Corrected
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)

	colNames = []string {
		"r/r500c",

		"temp-ew-m=1e14-z=0",
		"temp-ew-m=1e15-z=0",
		"temp-ew-m=1e14-z=0.5",
		"temp-ew-m=1e15-z=0.5",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]

	outTable := table.NewOutTable(colNames...)

	cFunc0 := halo.ConcentrationFunc(halo.Bhattacharya2013, 0.0)
	cFunc5 := halo.ConcentrationFunc(halo.Bhattacharya2013, 0.5)

	h014 := halo.New(fTh, simPpt, cFunc0, halo.Corrected, 1e14, 0.0)
	h015 := halo.New(fTh, simPpt, cFunc0, halo.Corrected, 5e14, 0.0)
	h514 := halo.New(fTh, simPpt, cFunc5, halo.Corrected, 1e14, 0.5)
	h515 := halo.New(fTh, simPpt, cFunc5, halo.Corrected, 5e14, 0.5)

	minFracLog, maxFracLog := math.Log10(0.01), math.Log10(10)
	logWidth := (maxFracLog - minFracLog) / steps
	for fracLog := minFracLog; fracLog <= maxFracLog; fracLog += logWidth {
		frac := math.Pow(10, fracLog)

		println(
			h014.EWTemperature(halo.Corrected, valPpt, h015.C500.R * frac),
			h015.EWTemperature(halo.Corrected, valPpt, h015.C500.R * frac),
			h514.EWTemperature(halo.Corrected, valPpt, h514.C500.R * frac),
			h515.EWTemperature(halo.Corrected, valPpt, h515.C500.R * frac))

		outTable.AddRow(frac,
			h014.EWTemperature(halo.Corrected, valPpt, h015.C500.R * frac),
			h015.EWTemperature(halo.Corrected, valPpt, h015.C500.R * frac),
			h514.EWTemperature(halo.Corrected, valPpt, h514.C500.R * frac),
			h515.EWTemperature(halo.Corrected, valPpt, h515.C500.R * frac))
	}

	outTable.Write(table.KeepHeader, path.Join(outDir, "radial-temp.table"))
}