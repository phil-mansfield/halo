package main

import (
	"os"
	"path"
	"math"

	"bitbucket.org/phil-mansfield/halo"
	"bitbucket.org/phil-mansfield/table"
)

const (
	steps = 300

	simFth = halo.Battaglia2013
	simPpt = halo.BattagliaAGN2012
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)

	alphaNames = []string {
		"r/r500",

		"alpha-m=1e14-z=0",
		"alpha-m=5e14-z=0",
		"alpha-m=1e14-z=0.5",
		"alpha-m=5e14-z=0.5",
	}

	betaNames = []string {
		"r/r500c",

		"beta-m=1e14-z=0",
		"beta-m=5e14-z=0",
		"beta-m=1e14-z=0.5",
		"beta-m=5e14-z=0.5",
	}

	fthNames = []string {
		"r/r500c",

		"1/fTh-m=1e14-z=0",
		"1/fTh-m=5e14-z=0",
		"1/fTh-m=1e14-z=0.5",
		"1/fTh-m=5e14-z=0.5",
	}

	biasNames = []string {
		"r/r500c",

		"bias-m=1e14-z=0",
		"bias-m=5e14-z=0",
		"bias-m=1e14-z=0.5",
		"bias-m=5e14-z=0.5",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]

	alphaTable := table.NewOutTable(alphaNames...)
	betaTable := table.NewOutTable(betaNames...)
	fthTable := table.NewOutTable(fthNames...)
	biasTable := table.NewOutTable(biasNames...)

	cFunc0 := halo.ConcentrationFunc(halo.Bhattacharya2013, 0)
	cFunc5 := halo.ConcentrationFunc(halo.Bhattacharya2013, 0.5)

	h014 := halo.New(fTh, simPpt, cFunc0, halo.Corrected, 1e14, 0)
	h015 := halo.New(fTh, simPpt, cFunc0, halo.Corrected, 5e14, 0)
	h214 := halo.New(fTh, simPpt, cFunc5, halo.Corrected, 1e14, 0.2)
	h215 := halo.New(fTh, simPpt, cFunc5, halo.Corrected, 5e14, 0.2)

	minFracLog, maxFracLog := math.Log10(0.01), math.Log10(10)
	logWidth := (maxFracLog - minFracLog) / steps
	for fracLog := minFracLog; fracLog <= maxFracLog; fracLog += logWidth {
		frac := math.Pow(10, fracLog)
		r014 := h014.C500.R * frac
		r015 := h015.C500.R * frac
		r214 := h214.C500.R * frac
		r215 := h215.C500.R * frac

		alphaTable.AddRow(frac,
			h014.AlphaBias(r014),
			h015.AlphaBias(r015),
			h214.AlphaBias(r214),
			h215.AlphaBias(r215))

		betaTable.AddRow(frac,
			h014.BetaBias(r014),
			h015.BetaBias(r015),
			h214.BetaBias(r214),
			h215.BetaBias(r215))

		fthTable.AddRow(frac,
			h014.FThermal(r014),
			h015.FThermal(r015),
			h214.FThermal(r214),
			h215.FThermal(r215))

		biasTable.AddRow(frac,
			h014.BFrac(r014),
			h015.BFrac(r015),
			h214.BFrac(r214),
			h215.BFrac(r215))
	}

	alphaTable.Write(table.KeepHeader, path.Join(outDir, "radial-alpha.table"))
	betaTable.Write(table.KeepHeader, path.Join(outDir, "radial-beta.table"))
	fthTable.Write(table.KeepHeader, path.Join(outDir, "radial-fth.table"))
	biasTable.Write(table.KeepHeader, path.Join(outDir, "radial-bias.table"))
}