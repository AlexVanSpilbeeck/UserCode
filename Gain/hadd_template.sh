#!/bin/bash

#source /sandbox/cmss/cmsset_default.sh
source /afs/cern.ch/cms/cmsset_default.sh

storedir=STOREDIR
storeoutput=T2_TMP_DIR/hadd_output
runningdir=RUNNINGDIR

#cp ${runningdir}/hadd.txt .
cd ${runningdir}
eval `scramv1 runtime -sh`

if [ -d ${storeoutput} ] ; then rm -r ${storeoutput} ; fi
mkdir -p ${storeoutput}

T2_RM ${storedir}/GainCalibration.root
hadd -f ${storeoutput}/GainCalibration.root `cat hadd.txt`
T2_CP ${storeoutput}/GainCalibration.root ${storedir}/GainCalibration.root 

