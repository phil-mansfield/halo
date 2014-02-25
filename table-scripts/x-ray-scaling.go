package main

import (
	"os"
	"path"
	"math"

	"bitbucket.org/phil-mansfield/halo"
	"bitbucket.org/phil-mansfield/table"
)

const (
	z = 0.0

	steps = 200
	simFth = halo.Battaglia2013
	ppt = halo.Planck2012

	bpbt = halo.ThermalPressure
	cpbt = halo.EffectivePressure
)

var (
	fTh = halo.FThermalFunc(simFth, halo.MeanCurve)
	cFunc = halo.ConcentrationFunc(halo.Bhattacharya2013, z)

	tempColNames = []string {
		"m500c-b",
		"m500c-c",
		"Tew500c-b",
		"Tew500c-c",
		"c/b",
	}

	yColNames = []string {
		"m500c-b",
		"m500c-c",
		"Y500c-b",
		"Y500c-c",
		"c/b",
	}

	fGasColNames = []string {
		"m500c-b",
		"m500c-c",
		"fg500c-b",
		"fg500c-c",
		"c/b",
	}
)

func fGas(h *halo.Halo, bt halo.BiasType) float64 {
	if bt == halo.Corrected {
		return h.GasEnclosed(bt, cpbt, ppt, h.C500.R) / 
			h.MassEnclosed(bt, h.C500.R)
	} else {
		return h.GasEnclosed(bt, bpbt, ppt, h.R500cBias) / 
			h.MassEnclosed(bt, h.R500cBias)
	}	
}

func main() {
	if len(os.Args) != 2 {
		panic("You must give exactly one argument.")
	}
	outDir := os.Args[1]
	tempTable := table.NewOutTable(tempColNames...)
	yTable := table.NewOutTable(yColNames...)
	fGasTable := table.NewOutTable(fGasColNames...)

	minMassLog, maxMassLog := math.Log10(1e13), math.Log10(1e15)
	logWidth := (maxMassLog - minMassLog) / steps

	for massLog := minMassLog; massLog <= maxMassLog; massLog += logWidth {
		bMass := math.Pow(10, massLog)

		h := halo.New(fTh, ppt, cFunc, halo.Biased, bMass, z)


		bTemp := h.EWTemperature(halo.Biased, bpbt, ppt, h.R500cBias)
		cTemp := h.EWTemperature(halo.Corrected, cpbt, ppt, h.C500.R)
		tempTable.AddRow(bMass, h.C500.M, bTemp, cTemp, cTemp / bTemp)

		bY := h.ThompsonY(bpbt, ppt, h.R500cBias)
		cY := h.ThompsonY(cpbt, ppt, h.C500.R)
		yTable.AddRow(bMass, h.C500.M, bY, cY, cY / bY)

		bFGas := fGas(h, halo.Biased)
		cFGas := fGas(h, halo.Corrected)
		fGasTable.AddRow(bMass, h.C500.M, bFGas, cFGas, cFGas / bFGas)
	}

	tempTable.Write(table.KeepHeader, path.Join(outDir, "mass-temp.table"))
	yTable.Write(table.KeepHeader, path.Join(outDir, "mass-y.table"))
	fGasTable.Write(table.KeepHeader, path.Join(outDir, "mass-fgas.table"))
}
