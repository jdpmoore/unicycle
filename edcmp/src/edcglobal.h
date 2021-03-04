c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	GLOBAL CONSTANTS
c	================
c
c	NRECMAX = max. number of observation positions
c	NZMAX = max. number of the discrete source depths
c	NRMAX = max. number of the discrete radial diatances
c	NSMAX = max. number of the source rectangles
c	NPSMAX = max. number of discrete point sources each depth
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer NZMAX,NRMAX,NSMAX,NPSMAX,NRECMAX,NFIELDS
	parameter(NZMAX=300,NRMAX=4001)
	parameter(NSMAX=1000,NPSMAX=500000)
	parameter(NRECMAX=270000)
	parameter(NFIELDS=4)
