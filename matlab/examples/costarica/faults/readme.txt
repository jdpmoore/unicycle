# reference point
echo -85.56 9.76 | proj +proj=utm +zone=16
657947.74	1079213.96

# Convert everything to .flt
awk '{for (i=0; i<19; i++){print $1}}' temp | paste sse07 -| proj +proj=utm +zone=16 | awk '{print ($1-657947.74)/1e3,($2-1079213.96)/1e3,$4,$3}' | awk 'function asin(x) { return atan2(x, sqrt(1-x*x)) } BEGIN{pi=atan2(1,0)*2}{for (i=1; i<=NF; i++){a[NR,i]= $i}}NF>p{p=NF} END {for (j=1; j<=NR ; j++) {k=int((j-1)/19); i=j-(k*19); if (i<19){len=sqrt((a[k*19+i+1,1]-a[k*19+i,1])**2+(a[k*19+i+1,2]-a[k*19+i,2])**2)} else {len=sqrt((a[k*19+i-1,1]-a[k*19+i,1])**2+(a[k*19+i-1,2]-a[k*19+i,2])**2)} if (k<18){width=sqrt((a[(k+1)*19+i,1]-a[k*19+i,1])**2+(a[(k+1)*19+i,2]-a[k*19+i,2])**2+(a[(k+1)*19+i,3]-a[k*19+i,3])**2)}else{width=sqrt((a[(k-1)*19+i,1]-a[k*19+i,1])**2+(a[(k-1)*19+i,2]-a[k*19+i,2])**2+(a[(k-1)*19+i,3]-a[k*19+i,3])**2)}; if (i<19){if (k<18){strike=(atan2(a[k*19+i+1,1]-a[k*19+i,1],a[k*19+i+1,2]-a[k*19+i,2])+atan2(a[(k+1)*19+i+1,1]-a[(k+1)*19+i,1],a[(k+1)*19+i+1,2]-a[(k+1)*19+i,2]))/2}else{strike=atan2(a[k*19+i+1,1]-a[k*19+i,1],a[k*19+i+1,2]-a[k*19+i,2])}} else {strike=atan2((a[k*19+i,1]-a[k*19+i-1,1]),(a[k*19+i,2]-a[k*19+i-1,2]))}; if (k<18){dip=asin((a[(k+1)*19+i,3]-a[k*19+i,3])/width)} else {dip=asin((a[k*19+i,3]-a[(k-1)*19+i,3])/width)} ;print j,a[j,4]/1e3,a[j,2]*1e3,a[j,1]*1e3,a[j,3]*1e3,len*1e3, width*1e3, strike*180/pi, dip*180/pi,43.8}}' > dixon+14_07.flt

# create .ned and .tri files
grep -v "#" Kyriakopoulos_etal_JGR_2015_slabInterface.xyz | proj +proj=utm +zone=16 | awk 'BEGIN{i=1}{if (100>$3 && 4<$3){print i,($2- 1079213.96),($1-657947.74),$3*1e3; i=i+1}}' > kyriakopoulos+15_receiver.ned
grep -v "#" kyriakopoulos+15_receiver.ned | awk '{print $2,$3}' | triangulate -V | awk '{print NR,$1+1,$2+1,$3+1,90}' > kyriakopoulos+15_receiver.tri
tri2vtp.sh -s 60e3 kyriakopoulos+15_receiver.tri > test.tri
mv test.tri kyriakopoulos+15_receiver.tri
tri2vtp.sh kyriakopoulos+15_receiver.tri

# create plate boundary
grep -v "#" plate_boundaries_ll.xy | proj +proj=utm +zone=16 -t">" | grep -v "*" | awk '{if (">"==substr($0,0,1)){print $0}else{printf "%12.2f %12.2f\n", $1-657947.74,$2-1079213.96,0}}' > plate_boundaries.xy
