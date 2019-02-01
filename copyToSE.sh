#! /bin/bash
# gfal-rm -r gsiftp://t3se01.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/ineuteli/analysis/LQ_2017/DY

CHANNEL="mutau"
YEAR="2018"
while getopts "c:y:" option; do case "${option}" in
  c) CHANNELS="${OPTARG}";;
  y) YEAR="${OPTARG}";;
esac done
function peval { echo ">>>   $@"; eval "$@"; }
#IGNORE="EES MES LTF JTF"
#SELECT="TES"
#GREP=""

# for pattern in $IGNORE; do
#   FILES0=`ls $FILES0 | grep -v $pattern`
# done
# 
# [[ ! $SELECT ]] && FILES=$FILES0
# for pattern in $SELECT; do
#   FILES+=" "`ls $FILES0 | grep $pattern`
# done
# echo ">>> files to be copied:"
# echo "$FILES"

for year in $YEAR; do
  [[ $year = '#'* ]] && continue
  
  OUTPUT="ineuteli/analysis/LQ_$year"
  PNFS_OUTPUT="/pnfs/psi.ch/cms/trivcat/store/user/$OUTPUT"
  XRD="root://t3dcachedb.psi.ch:1094"
    
  cd "/scratch/$OUTPUT"
  FILES=`ls */*.root`
  DIRS=`ls /scratch/$OUTPUT`
  
  # MAKE DIR
  TRY=0
  while [ ! -e $PNFS_OUTPUT -a $TRY -lt 4 ]; do
    # Warning: gfal tools are broken by initialization of CMSSW environment!
    if [ $TRY == 0 ]; then
      echo ">>> Warning! $PNFS_OUTPUT does not exist..."
      peval "gfal-mkdir -p gsiftp://t3se01.psi.ch/$PNFS_OUTPUT"
    else
      echo ">>> check again..."
    fi
    sleep 10; TRY=$((TRY+1))
  done
  [ ! -e $PNFS_OUTPUT ] && echo ">>>   ERROR! $PNFS_OUTPUT still does not exist after $TRY attempts..." && exit 1
  
  # MAKE DIR
  for d in $DIRS; do
    dir="$PNFS_OUTPUT/${d%%/}"
    if [ ! -e $dir ]; then
      # Warning: gfal tools are broken by initialization of CMSSW environment!
      echo ">>> Warning! $dir does not exist..."
      peval "gfal-mkdir -p gsiftp://t3se01.psi.ch/$dir"
      sleep 2;
    fi
  done
  
  # CHECK DIR FEW TIMES
  for d in $DIRS; do
    dir="$PNFS_OUTPUT/${d%%/}"
    TRY=0
    while [ ! -e $dir -a $TRY -lt 3 ]; do
      # Warning: gfal tools are broken by initialization of CMSSW environment!
      [ $TRY == 0 ] && echo ">>> Warning! $dir does not exist..." || echo ">>> check again..."
      sleep 5; TRY=$((TRY+1))
    done
    [ ! -e $dir ] && echo ">>>   ERROR! $dir does not exist..." && exit 1
  done
    
  i=0
  N=`echo $FILES | wc -w`
  for f in $FILES; do
    i=$((i+1))
    CMD="xrdcp -f $f ${XRD}/${PNFS_OUTPUT}/$f"
    echo
    echo ">>> ${i}/${N}: ${f}"
    peval "$CMD"
  done
  echo

done
