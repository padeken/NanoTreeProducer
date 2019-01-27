#! /bin/bash
# gfal-rm -r gsiftp://t3se01.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/ineuteli/analysis/LQ_2017/DY

YEAR="#2017 2018"
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
  if [ ! -e $PNFS_OUTPUT ]; then
    # Warning: gfal tools are broken by initialization of CMSSW environment!
    echo ">>> Warning! $PNFS_OUTPUT does not exist..."
    CMD="gfal-mkdir -p gsiftp://t3se01.psi.ch/$PNFS_OUTPUT"
    echo ">>>   $CMD"
    $CMD
    sleep 6
    [ ! -e $PNFS_OUTPUT ] && echo ">>>   ERROR! $PNFS_OUTPUT does not exist..." && exit 1
  fi
  
  # MAKE DIR
  for d in $DIRS; do
    dir="$PNFS_OUTPUT/${d%%/}"
    if [ ! -e $dir ]; then
      # Warning: gfal tools are broken by initialization of CMSSW environment!
      echo ">>> Warning! $dir does not exist..."
      CMD="gfal-mkdir -p gsiftp://t3se01.psi.ch/$dir"
      echo ">>>   $CMD"
      $CMD
      sleep 6
      [ ! -e $dir ] && echo ">>>   ERROR! $dir does not exist..." && exit 1
    fi
  done
  
  i=0
  N=`echo $FILES | wc -w`
  for f in $FILES; do
    i=$((i+1))
    CMD="xrdcp -f $f ${XRD}/${PNFS_OUTPUT}/$f"
    echo
    echo ">>> ${i}/${N}: ${f}"
    echo ">>> $CMD"
    $CMD
  done
  echo

done
