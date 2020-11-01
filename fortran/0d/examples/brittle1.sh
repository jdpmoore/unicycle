#!/bin/bash

# computes the slip evolution of a spring-slider system with radiation damping
# the spring stiffness is given by 
#
#    k = 7 pi / 16 G / R, 
#
# where R is the radius of the asperity

selfdir=$(dirname $0)

WDIR=$selfdir/brittle1

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time unicycle-0d-ratestate --epsilon 1e-6 --maximum-step 3.15e6 --maximum-iterations 200000 $* <<EOF
# output directory
$WDIR
# elastic moduli
30e3
# time interval
3.15e9
# Vpl radius
 1e-9   1e3
# mu0 sig    a     b    L   Vo  G/2Vs
  0.6 1e2 0.01 0.014 1e-3 1e-6      5
EOF

