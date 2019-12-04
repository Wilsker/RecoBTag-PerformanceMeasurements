#!/bin/bash
numbins=4 #ptrange-1 entries
ini=0
fin=$numbins
ptrange=(30 50 70 100 140)
stat=("0 1 2 3 4 5 6 7 8 9 10 11")
#Go in the folder with nominal values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep central | awk '{printf "%.12f\n", $11}'
sfval=(
0.975910604000
0.973668038845
0.986115515232
0.978846848011
0.950889110565
0.959384560585
0.981989026070
0.961818277836
0.924486041069
0.951390326023
0.973088264465
0.980542778969
)
#Go in the folder with nominal values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep down/up | awk '{printf "%.12f\n", $11}'
sfdevtoterr=(
0.975910604000
0.973668038845
0.986115515232
0.978846848011
0.950889110565
0.959384560585
0.981989026070
0.961818277836
0.924486041069
0.951390326023
0.973088264465
0.980542778969
)
#Go in the folder with varied values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep central | awk '{printf "%.12f\n", $11}'
sfvalnewsys=(
0.958191871643
0.959637582302
0.988767743111
0.968010842800
0.929773688316
0.942144930363
0.962891995907
0.952567398548
0.889630854130
0.921991944313
0.945743024349
0.962389469147
)
for (( wp=0; wp<3; wp++)); do
 echo "wp $wp"
 #Statistics 
 p=0
 l=0
 r=1
 for p in $stat; do
  if [[ "$ini" -le "$p" ]] && [[ "$p" -lt "$fin" ]]; then
   #Get values of sys variation
   devnewsys=$(awk "BEGIN {printf \"%.12f\",(sqrt((${sfval[$p]}-${sfvalnewsys[$p]})*(${sfval[$p]}-${sfvalnewsys[$p]})))}") 
   relerror=$(awk "BEGIN {printf \"%.12f\",$devnewsys/${sfval[$p]}*100}")
   echo $relerror
  fi
  let p=p+1
 done
 let ini=ini+numbins
 let fin=fin+numbins
done
