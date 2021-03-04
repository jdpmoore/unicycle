c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
c	global parameters:
c	  nzmax: max. interface index;
c	  lmax: max. no of total homogeneous layers (lmax <= nzmax-2);
c	  nrmax: max. no of traces;
c	  nzsmax: max. no of source depth samples.
c
        integer nzmax,lmax,nrmax,nzsmax 
        parameter(nzmax=1002,lmax=1000,nrmax=4001,nzsmax=300)
c
c	index parameters for bessel function table
c
	integer nbess,ndbess
	parameter(nbess=512,ndbess=256)
c
c	earth gravity in cm/s^2
c
	double precision g0
c	parameter(g0=0.982d+03)
	parameter(g0=0.d0)
