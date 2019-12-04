#!/bin/bash
numbins=7 #ptrange-1 entries
ini=0
fin=$numbins
ptrange=(30 50 70 100 140 200 300 600)
stat=("0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20")
#Go in the folder with nominal values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep central | awk '{printf "%.12f\n", $11}'
sfval=(
0.979518592358
0.977694332600
0.988283693790
0.985508024693
0.974605739117
0.966486811638
0.981102526188
0.954221308231
0.972418606281
0.985319554806
0.976234972477
0.956586420536
0.961752116680
0.938968122005
0.932330667973
0.958764374256
0.978433191776
0.980090618134
0.979581475258
0.975249826908
1.046803474430
)
#Go in the folder with nominal values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep down/up | awk '{printf "%.12f\n", $11}'
sfdevtoterr=(
0.979518592358
0.977694332600
0.988283693790
0.985508024693
0.974605739117
0.966486811638
0.981102526188
0.954221308231
0.972418606281
0.985319554806
0.976234972477
0.956586420536
0.961752116680
0.938968122005
0.932330667973
0.958764374256
0.978433191776
0.980090618134
0.979581475258
0.975249826908
1.046803474430
)
#Go in the folder with varied values
#cp kin_calib.csv kin_calib_awk.csv and 0,$s/"//g in kin_calib_awk.csv
#cat kin_calib_awk.csv | grep central | awk '{printf "%.12f\n", $11}'
sfvalnewsys=(

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
   echo $devnewsys
   #Get updated values of total sys down/up
   originaluncert=$(awk "BEGIN {printf \"%.12f\",(sqrt((${sfval[$p]}-${sfdevtoterr[$p]})*(${sfval[$p]}-${sfdevtoterr[$p]})))}")
   newuncert=$(awk "BEGIN {printf \"%.12f\",(sqrt(($devnewsys*$devnewsys)+($originaluncert*$originaluncert)))}") 
   echo $originaluncert $devnewsys $newuncert
   valsys=$(awk "BEGIN {printf \"%.12f\",${sfval[$p]}+${variations[$p]}*${sfval[$p]}}") # }")
   echo "$wp, kin (emu), up_tWth, 0, -2.4, 2.4, ${ptrange[$l]}, ${ptrange[$r]}, 0, 1, \""$valsys"\""
   let l=l+1
   let r=r+1
  fi
  let p=p+1
 done
 let ini=ini+numbins
 let fin=fin+numbins
done
