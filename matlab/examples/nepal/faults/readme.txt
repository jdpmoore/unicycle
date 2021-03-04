
# reference point
echo 84.708 28.147 | proj +proj=utm +zone=45
274917.78       3115610.96

# convert the model of Wei (2015) to the .flt format
grep -v "#" static_out.dat | proj +proj=utm +zone=45 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip         x1         x2        x3   length  width strike   dip  rake"}{slip=$4/1e2;rake=$5;str0=$6;dip0=$7;L=$10;W=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-3115610.96)/1e3-s[1]*L/2;x2=($1-274917.78)/1e3-s[2]*L/2;x3=$3;printf "%3d %6.2f %10.2f %10.2f %9.1f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > wei+15_1.flt

# convert models to .kml files
makecpt -C/Users/sbarbot/Documents/src/relax/share/jet.cpt -T0/3.106/0.01 -Z -N > hayes15.cpt
flt2kml.sh -x 274917.78 -y 3115610.96 -z 45 -C hayes15.cpt hayes15.flt
makecpt -C/Users/sbarbot/Documents/src/relax/share/jet.cpt -T0/4.211/0.01 -Z -N > wei+15_1.cpt
flt2kml.sh -x 274917.78 -y 3115610.96 -z 45 -C wei+15_1.cpt wei+15_1.flt

# convert GOCAD MFT_ model to .ned/.tri
grep VRTX MFT_surface.ts | awk '{printf "%04.4d %f %f %f\n",NR,($4-3115610.96),($3-274917.78),$5}' > MFT_surface.ned
grep TRGL MFT_surface.ts | awk '{printf "%04.4d 0 %04.4d %04.4d %04.4d 90\n",NR,$2,$3,$4}' > MFT_surface.tri

grep TRGL MFT_v20.ts | awk '{printf "%04.4d %04.4d %04.4d %04.4d 90\n",NR,$2,$3,$4}' > MFT_receiver_v20.tri

# filter out distant mesh points
awk -v s=-70 'NR==FNR {i=0;for (j=2;j<=4;j++){a[j-1,0+$1]=$j}; next} {pi=atan2(1,0)*2;x1=a[1,$2];x2=a[2,$2];d=x1*cos(pi*s/180)+x2*sin(pi*s/180);if (20e3>d && -190e3<d){i=i+1;print i,$2,$3,$4,$5}}' MFT_receiver_v20.ned MFT_receiver_v20.tri > MFT_receiver_v20_short.tri
