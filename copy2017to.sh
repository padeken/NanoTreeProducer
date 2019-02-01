#! /bin/bash

CHANNEL="mutau"
YEAR="2018"
while getopts "c:y:" option; do case "${option}" in
  c) CHANNELS="${OPTARG}";;
  y) YEAR="${OPTARG}";;
esac done

[ "$YEAR" != "2016" -o "$YEAR" != ] && echo ">>> ERROR! Year $YEAR not valid!" && exit 1

SCRATCHOLD="/scratch/ineuteli/analysis/LQ_2017"
SCRATCHNEW="/scratch/ineuteli/analysis/LQ_$YEAR"

if [ "$YEAR" == "2016" ]; then
  SAMPLES="
    DY/DY3JetsToLL_M-50
    DY/DY4JetsToLL_M-50
    WJ//W4JetsToLNu
    ST/ST_t-channel_top
  "
else
  SAMPLES="
    DY/DYJetsToLL_M-10to50
    DY/DY4JetsToLL_M-50
    WJ//WJetsToLNu
    ST/ST_t-channel_top
    ST/ST_t-channel_antitop
  "
fi

for samplename in $SAMPLES; do
  [[ $samplename = '#'* ]] && continue
  echo $samplename
  
  fileold="$SCRATCHOLD/${samplename}_${CHANNEL}.root"
  filenew="$SCRATCHNEW/${samplename}_${CHANNEL}.root"  
  [ ! -e $fileold ] && echo "Warning! $fileold does not exist!" && continue 
  
  #echo "Copying $fileold to $filenew!"
  cp -v $fileold $filenew
  
done
