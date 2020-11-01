
# reference point
echo 130.726 32.782 | proj +proj=utm +zone=52
661641.27       3628438.35

# convert the Hayes model to the Unicycle format
grep -v "#" 20005iis.param | proj +proj=utm +zone=52 -r | awk '{if(NR>5){$1=($1-661641.27)/1e3;$2=($2-3628438.35)/1e3;print $0}}' | awk 'BEGIN{print "# nb slip(m)      x1(m)          x2         x3    length(m)   width strike  dip   rake";pi=atan2(1,0)*2}{str=$6*pi/180;dip=$7/180*pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;n[1]=sin(str)*cos(dip);n[2]=-cos(str)*cos(dip);d[3]=sin(dip);slip=$4/1e2;rake=$5;len=5;wid=2.9;x1=$2-s[1]*len/2;x2=$1-s[2]*len/2;x3=$3-s[3]*len;printf("%3.3d  %6.3f %11.2f %11.2f %10.4f      %5.2f %5.2f %6.1f %3.1f %6.2f\n",NR,slip,x1*1e3,x2*1e3,x3*1e3,len*1e3,wid*1e3,str*180/pi,dip*180/pi,rake)}' > ./faults/kumamoto-hayes16-prel.flt


# Convert Chi-Hsien's data to the Unicycle format
cat disp_pos3m.dat | awk '{print $2, $3, $1,$4, $5,$6,$7,$8,$9}' | proj +proj=utm +zone=52 | awk 'BEGIN{printf("%s %6s %6s %10s %11s %11s %11s %11s %11s %11s\n","ID","stationname", "x","y", "dE", "dN", "dU", "SE", "SN", "SU" )}{printf("%4.4d  %s %12.4f %12.4f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n", NR, $3,($2-3628438.35), ($1-661641.27), $4, $5, $6, $7, $8, $9)}' > disp_pos3m_gps.dat 

# Create the receiver _patch.flt file

awk 'NR<=2 {print} NR>2 {printf "%3s %10.2f %10.2f %10.4f %10.2f %10.2f %8.2f %8.2f %9.2f\n", $1,$3,$4,$5,$6,$7,$8,$9,$10}' swei16.flt > swei16_patch.flt