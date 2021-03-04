
# reference point
echo 166.5 -46.0 | proj +proj=utm +zone=59
151566.87	-5103900.18

# convert the slip model of Beavan et al. (2010) to the Relax/Unicycle format
grep -v "#" Beavan2010CP3.dat | proj +proj=utm +zone=59 -r | awk 'BEGIN{Pi=3.14159265;print "# n    slip         x1         x2        x3   length  width strike   dip  rake"}{slip=$8;rake=$5;str0=$3;dip0=$4;L=$6;W=$7;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2+5103900.18)/1e3-s[1]*L/2+d[1]*W/2;x2=($1-151566.87)/1e3-s[2]*L/2+d[2]*W/2;x3=$9+d[3]*W/2;printf "%03.3d %7.3f %10.5f %10.5f %9.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1,x2,x3,L,W,str0,dip0,rake}' > beavan+10_km.flt
grep -v "#" Beavan2010CP3.dat | proj +proj=utm +zone=59 -r | awk 'BEGIN{Pi=3.14159265;print "# n    slip           x1           x2          x3   length  width strike   dip  rake"}{slip=$8;rake=$5;str0=$3;dip0=$4;L=$6;W=$7;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2+5103900.18)/1e3-s[1]*L/2+d[1]*W/2;x2=($1-151566.87)/1e3-s[2]*L/2+d[2]*W/2;x3=$9+d[3]*W/2;printf "%03.3d %7.3f %12.5f %12.5f %11.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > beavan+10.flt

grep -v "#" Beavan2010CP3.dat | proj +proj=utm +zone=59 -r | awk 'BEGIN{Pi=3.14159265;print "# n      x1           x2          x3   length  width strike   dip  rake"}{slip=$8;rake=$5;str0=$3;dip0=$4;L=$6;W=$7;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2+5103900.18)/1e3-s[1]*L/2+d[1]*W/2;x2=($1-151566.87)/1e3-s[2]*L/2+d[2]*W/2;x3=$9+d[3]*W/2;printf "%03.3d %12.5f %12.5f %11.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > beavan+10_receiver_patch.flt



grep -v "#" Beavan2010CP1.dat | proj +proj=utm +zone=59 -r | awk 'BEGIN{Pi=3.14159265;print "# n    slip           x1           x2          x3   length  width strike   dip  rake"}{slip=$8;rake=$5;str0=$3;dip0=$4;L=$6;W=$7;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2+5103900.18)/1e3-s[1]*L/2+d[1]*W/2;x2=($1-151566.87)/1e3-s[2]*L/2+d[2]*W/2;x3=$9+d[3]*W/2;printf "%03.3d %7.3f %12.5f %12.5f %11.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > beavan+10_CP1.flt


grep -v "#" Beavan2010CP1.dat | proj +proj=utm +zone=59 -r | awk 'BEGIN{Pi=3.14159265;print "# n      x1           x2          x3   length  width strike   dip  rake"}{slip=$8;rake=$5;str0=$3;dip0=$4;L=$6;W=$7;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2+5103900.18)/1e3-s[1]*L/2+d[1]*W/2;x2=($1-151566.87)/1e3-s[2]*L/2+d[2]*W/2;x3=$9+d[3]*W/2;printf "%03.3d %12.5f %12.5f %11.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > beavan+10_CP1_receiver_patch.flt                                              
