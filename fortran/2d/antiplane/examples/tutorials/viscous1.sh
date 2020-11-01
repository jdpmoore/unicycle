#!/bin/bash

# earthquake cycle on a vertical strike-slip fault 
# in antiplane strain in a viscoelastic medium
# with the radiation damping approximation

selfdir=$(dirname $0)

WDIR=$selfdir/viscous1

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

WIDTH=1e2
NW=400
GEOMETRY=`echo "" | awk -v w=$WIDTH -v nw=$NW  \
	'{for (j=1;j<=nw;j++){ 
		Vpl=1e-9; 
		x2=0;
		z3=(j-1)*w;
		x3=(z3>=25e3)?z3+10e3:z3;
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

N2=50
N3=50
WIDTH=2e2
THICKNESS=2e2;
MANTLE=`echo "" | awk -v n2=$N2 -v n3=$N3 -v thickness=$THICKNESS -v width=$WIDTH \
	'{for (j3=1;j3<=n3;j3++){ 
		x3=25e3+(j3-1)*width;
		for (j2=1;j2<=n2;j2++){
			e12=1e-15;
			e13=0;
			x2=(j2-n2/2)*thickness;
			dip=90;
			printf "%4d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d\n", 
			j2+(j3-1)*n2, e12, e13, x2, x3, width, thickness, dip
		}
	}}'`

RHEOLOGY=`echo "$MANTLE" | awk -v nw=$NW 'function min(x,y){return (x<y)?x:y}{ 
	depth=$5; 
	Gm=1;
	gammadot0m=1e-12;
	n=1;
	Q=0;
	R=8.314;
	printf "%3d %10.2e %10.2e %10.2e %10.2e %10.2e\n", 
	         NR,   Gm,   gammadot0m,     n,     Q,     R}'`

TEMP=`echo "$MANTLE" | awk 'function min(x,y){return (x<y)?x:y}{ 
	depth=$5; 
	rhoc=1e-6;
	To=1200;
	printf "%3d %10.2e %10.2e\n", 
	         NR,   rhoc,   To}'`

OBS=`echo "" | awk -v n2=$N2 -v n3=$N3 \
	'{for (j3=1;j3<=n3;j3++){ 
		printf "%4d %d\n",j3,n2/2+(j3-1)*n2
	}}'`


cat <<EOF > $WDIR/in.param
# output directory
$WDIR
# elastic moduli
30e3 30e3
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
# number of cuboid strain volumes
`echo "$MANTLE" | wc -l`
# n       e12       e13        x2        x3    width thickness
`echo "$MANTLE"`
# number of nonlinear Maxwell strain volumes
`echo "$MANTLE" | wc -l`
# n        Gm gammadot0m         n         Q         R
`echo "$RHEOLOGY"`
# number of thermal strain volumes
`echo "$MANTLE" | wc -l`
# n     rhoc  temperature
`echo "$TEMP"`
# number of observation patches
4
# n   index (2, 5, 10, and 20 km depth)
  1      20
  2      50
  3     100
  4     200
# number of observation volumes
`echo "$OBS" | wc -l`
# n    index
`echo "$OBS"`
# number of observation points
5
# n name  x2 x3
  1 0001 5e2  0
  2 0002 5e3  0
  3 0003 1e4  0
  4 0004 2e4  0
  5 0005 4e4  0
EOF

#time mpirun -n 2 unicycle-ap-viscouscycles \
#	--epsilon 1e-6 \
#	--export-netcdf \
#	--maximum-step 3.15e6 \
#	$* $WDIR/in.param



echo "# exporting to $WDIR/e12.grd"
xRange=`ls $WDIR/volume-*.dat | wc -l | awk '{print $1}'`
for i in $WDIR/volume-*.dat; do
	if [ "" == "$yRange" ]; then
		yRange=`wc -l $i | awk '{print $1}'`
	fi
done
c=1
for i in $WDIR/volume-*.dat; do
	if [ "" == "$yRange" ]; then
		yRange=`wc -l $i | awk '{print $1}'`
	fi
	cat $i | awk -v i=$c '{print i,NR,$6}'
	#cat $i | awk -v i=$c '{print i,NR,log(sqrt($6**2+$7**2))/log(10)}'
	c=`echo $c | awk '{print $1+1}'`
done | xyz2grd -I1/1 -R1/$xRange/1/$yRange -G$WDIR/e12.grd



