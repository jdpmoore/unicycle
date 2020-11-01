#!/bin/bash

# earthquake cycle on a long vertical strike-slip fault 
# in antiplane strain with the radiation damping approximation

selfdir=$(dirname $0)

WDIR=$selfdir/brittle1

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

WIDTH=1e2
NW=500
GEOMETRY=`echo "" | awk -v w=$WIDTH -v nw=$NW  \
	'{for (j=1;j<=nw;j++){ 
		Vpl=1e-9; 
		x2=0;
		x3=(j-1)*w;
		printf "%3d %9.2e %9.2e %9.2e %9.2e %3d\n", 
		j,Vpl,x2,x3,w,90}}'`

FRICTION=`echo "$GEOMETRY" | awk -v nw=$NW 'function min(x,y){return (x<y)?x:y}{ 
	depth=$4; 
	V=$2; 
	L=0.01; 
	a=1e-2; 
	b=(5e3<= depth && depth<=15e3)?a+4.0e-3:a-4e-3; 
	mu0=0.6; 
	sig=1e2; 
	Vo=1e-6; 
	printf "%3d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d\n", 
	         NR,   mu0,   sig,     a,     b,     L,    Vo, 5}'`

cat <<EOF > $WDIR/in.param
# output directory
$WDIR
# elastic moduli
30e3
# time interval
3.15e10
# number of patches
`echo "$GEOMETRY" | wc -l`
# n  Vpl           x2           x3    width      dip
`echo "$GEOMETRY"`
# number of frictional patches
`echo "$GEOMETRY" | wc -l`
# n  mu0       sig         a         b         L        Vo     G/2Vs
`echo "$FRICTION"`
# number of strain volumes
0
# number of observation patches
4
# n   index (2, 5, 10, and 20 km depth)
  1      20
  2      50
  3     100
  4     200
# number of observation volumes
0
# number of observation points
5
# n name  x2 x3
  1 0001 5e2  0
  2 0002 5e3  0
  3 0003 1e4  0
  4 0004 2e4  0
  5 0005 4e4  0
EOF

time mpirun -n 2 unicycle-ap-viscouscycles \
	--epsilon 1e-6 \
	--export-netcdf \
	--maximum-step 3.15e6 \
	$* $WDIR/in.param





