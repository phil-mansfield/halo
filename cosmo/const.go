package cosmo

const (
	OmegaM float64 = 0.27
	OmegaL         = 0.73
	OmegaB         = 0.0469
	OmegaR         = 0.0

	Sigma8 = 0.82

	H100 = 0.7
	H70  = H100 / 0.7

	EVMks  = 1.60217e-19
	KBMks  = 1.38066e-23
	MHyMks = 1.67223e-27
	MeMks  = 9.10939e-31

	H0Mks   = 3.24086e-18
	GMks    = 6.67259e-11
	MpcMks  = 3.08560e+22
	MSunMks = 1.98900e+30
	CMks    = 2.99792e+08

	SigmaTMks = 6.65246e-29

	XHy = 0.76
	YHe = 1.0 - XHy

	Mu         = 2.0*XHy + 3.0/4.0*YHe
	ElectronMu = XHy + 2.0/4.0*YHe
)
