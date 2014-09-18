#!/bin/bash

#source /sandbox/cmss/cmsset_default.sh
source /afs/cern.ch/cms/cmsset_default.sh

startdir=T2_TMP_DIR/gain_calib_NUM
mkdir -p $startdir

mydir=CFGDIR
storedir=STOREDIR
indir=INDIR

cd $mydir

echo "Setting the environment ..."
eval `scramv1 runtime -sh`

cp Run_offline_DQM_NUM_cfg.py $startdir/
cd $startdir
echo -e "************************"
echo -e "  => ls: \n`ls`"
echo -e "************************\n\n"

echo -e "Copying file from storage to local ..."
file=GainCalibration_NUM_RUN.EXT
echo "(T2_CP $indir/$file $file)"
T2_CP $indir/$file $file
echo -e "************************"
echo -e "  => ls: \n`ls`"
echo -e "************************\n\n"

cat Run_offline_DQM_NUM_cfg.py |\
  sed "s#FILENAME#file:$file#" > torun.py
cat torun.py

echo -e "\n\n Running CMSSW job:"
cmsRun torun.py
cat *.log

echo -e "Copying output to pnfs:"
echo "(T2_CP "DQM*.root" ${storedir}/NUM.root)"
T2_CP DQM*.root ${storedir}/NUM.root

echo -e "end ... \n\n\n"

cp -f Summary.txt ${mydir}/TEXToutput/Summary_NUM.txt
cp -f SummaryPerDetID.txt ${mydir}/TEXToutput/SummaryPerDetID_NUM.txt
cp *.log ${mydir}/JOB_NUM/
cp *.txt ${mydir}/JOB_NUM/
cd $mydir
rm -fr $startdir

