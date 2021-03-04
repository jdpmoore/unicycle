
# reference point
echo -155.273256184167 19.339057464198 | proj +proj=utm +zone=5
261179.55       2139913.61

# convert the slip distribution of Montgomery-Brown (2013) to .flt format
./transpose.sh PM.txt | paste - <(  awk '{print $1}' MHAT.dat  ) | awk 'BEGIN{print "# n  slip(m)        x1(km)            x2           x3   length(km)        width strike    dip   rake";pi=atan2(1,0)*2}{rake = atan2 ($9,$8)/pi*180; printf "%3.3d  %7.4f %13.6e %13.6e %12.5e %12.5e %12.5e %6.1f %4.2f %6.2f\n", NR, $11,$7, $6, $3, $1, $2, $5, $4, rake}' > montgomery-brown+13_980220.flt

