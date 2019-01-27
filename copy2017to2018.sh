#! /bin/bash

SCRATCH2017="/scratch/ineuteli/analysis/LQ_2017"
SCRATCH2018="/scratch/ineuteli/analysis/LQ_2018"
CHANNEL="tautau"
SAMPLES="
TT/TTToSemiLeptonic
DY/DYJetsToLL_M-10to50
DY/DY4JetsToLL_M-50
WJ//WJetsToLNu
WJ//W1JetsToLNu
ST/ST_t-channel_top
ST/ST_t-channel_antitop
VV/WZ
"

for samplename in $SAMPLES; do
  
  filename1="$SCRATCH2017/${samplename}_${CHANNEL}.root"
  filename2="$SCRATCH2018/${samplename}_${CHANNEL}.root"  

  [ ! -e $filename ] && echo "Warning! $filename1 does not exist!" && continue 

  #echo "Copying $filename1 to $filename2!"
  cp -v $filename1 $filename2
  
done


