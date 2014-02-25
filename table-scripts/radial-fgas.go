package main


import (
	"os"
	"path"
	"math"

	"bitbucket.org/phil-mansfield/halo"
	"bitbucket.org/phil-mansfield/halo/cosmo"
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

	fgasNames = []string {
		"r/r500c",

		"fgas-m=1e14-z=0",
		"fgas-m=1e15-z=0",
		"fgas-m=1e14-z=0.5",
		"fgas-m=1e15-z=0.5",
	}

	normNames = []string {
		"r/r500c",

		"fgas-m=1e14-z=0/(OmegaB/OmegaM)",
		"fgas-m=1e15-z=0/(OmegaB/OmegaM)",
		"fgas-m=1e14-z=0.5/(OmegaB/OmegaM)",
		"fgas-m=1e15-z=0.5/(OmegaB/OmegaM)",
	}
)

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]

	println("WARNING: You still need to implement this quickly :3")

	fgasTable := table.NewOutTable(fgasNames...)
	normTable := table.NewOutTable(normNames...)

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
		r14 := h014.C500.R * frac
		r15 := h015.C500.R * frac

		frac014 := h014.GasEnclosed(bt, valPpt, r14) /
			h014.MassEnclosed(bt, r14)
		frac015 := h015.GasEnclosed(bt, valPpt, r15) /
			h015.MassEnclosed(bt, r15)
		frac514 := h514.GasEnclosed(bt, valPpt, r14) /
			h514.MassEnclosed(bt, r14)
		frac515 := h515.GasEnclosed(bt, valPpt, r15) /
			h515.MassEnclosed(bt, r15)

		fgasTable.AddRow(frac, frac014, frac015, frac514, frac515)

		normTable.AddRow(frac,
			frac014 / (cosmo.OmegaB / cosmo.OmegaM),
			frac015 / (cosmo.OmegaB / cosmo.OmegaM),
			frac514 / (cosmo.OmegaB / cosmo.OmegaM),
			frac515 / (cosmo.OmegaB / cosmo.OmegaM))
	}

	fgasTable.Write(table.KeepHeader,
		path.Join(outDir, "radial-density-frac.table"))
	normTable.Write(table.KeepHeader,
		path.Join(outDir, "radial-density-frac-norm.table"))
}