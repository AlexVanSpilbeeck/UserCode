#!/bin/bash

########################
##
##   Script to run the gain calibration on all 40 FEDs
##
########################

#################
# Sourcing castor/T2 utilities
source utils.sh

#################
# Set the extension for input root files

ext=dmp


usage(){

echo 'Usage :
./Run.sh -create   RUNNUMBER INPUTDIR STOREDIR  : will create the needed directories, python files to run the calib.
./Run.sh -submit   RUNNUMBER : will launch the 40 calibration jobs
./Run.sh -resubmit RUNNUMBER iJOB: will resubmit job iJOB , using the submit_iJOB.sh in the working directory.
./Run.sh -stage    RUNNUMBER : will stage all files needed that are stored on castor.
./Run.sh -hadd     RUNNUMBER : will hadd all 40 output files of the calib jobs into one file
./Run.sh -summary  RUNNUMBER : will launch the summary job.
./Run.sh -pdf      RUNNUMBER : will recompile the latex file to recreate the pdf summary.
./Run.sh -compare  RUNNUMBER1 FILE1 RUNNUMBER2 FILE2
     OR  -compare  RUNNUMBER1 RUNNUMBER2 : only if you have run -create/-submit/-hadd for both runs
./Run.sh -payload  RUNNUMBER : will produce the payloads.
'
exit


}


set_parameters(){
  storedir=$storedir/GainRun_$run
  echo "run : $run"
  echo "indir : $indir"
  echo "storedir : $storedir"
  if [ "$run" == "" ] || [ "$indir" == "" ] || [ "$storedir" == "" ];then usage ; fi
  #run=$1
  #indir=$2
  #storedir=$3
  runningdir=`pwd`/Run_$run
}

write_config(){
  echo "run = $run" 		>  $runningdir/config
  echo "indir = $indir" 	>> $runningdir/config
  echo "storedir = $storedir" 	>> $runningdir/config
}

read_config(){
  config=Run_$run/config
  if [ ! -f $config ];then echo "No config found for run $run. Make sure Run_$run exist in `pwd` ..." ; exit ;fi
  indir=`cat $config|grep -e "indir ="|awk '{printf $3}'`
  storedir=`cat $config|grep -e "storedir ="|awk '{printf $3}'`
  echo -e "Reading config file $config"
  echo -e "run : $run"
  echo -e "indir : $indir"
  echo -e "storedir : $storedir\n"
  runningdir=`pwd`/Run_$run
}

make_dir(){
  if [ ! -d $1 ] ; then mkdir $1
  else rm -r $1/*
  fi
}

lock(){
  if [ -f .lock_gaincalib ];then
    echo "Another instance of Run.sh is already running, so wait for it to finish !"
    echo "In case it is really not the case (like you previously killed it), remove the file \".lock_gaincalib\""
    exit
  else
    touch .lock_gaincalib
  fi
}


create(){
  
  #making running directory
  make_dir ${runningdir}

  #cleaning
  rm -f filelist.txt es.log
  
  #cleaning output dir
  set_specifics ${storedir}
  $T2_RM -r ${storedir}
  $T2_MKDIR -p ${storedir}
  $T2_CHMOD 0777 ${storedir}/ #TOFIX ?
  
  #chmod -R 0777 $indir

  if [ `is_on_castor $indir` -eq 1 ] ; then wait_for_staging ; fi
}


make_file_list(){

  #Copying template specific to gain calib to the general one used by Run_offline_DQM.csh
  #cp -f gaincalib_template_cfg.py client_template_calib_cfg.py

  #if [ `is_on_castor $indir` -eq 1 ] ; then wait_for_staging ; fi

  #making python files
  touch filelist.txt
  for i in `seq 0 39`;do
    file=GainCalibration_${i}_$run.$ext
    if [ `is_file_present $indir/$file` -eq 0 ];then echo "File $file is not present in $indir ...";continue;fi 

    echo "$indir/$file" > filelist.txt
    ./Run_offline_DQM.csh filelist.txt Calibration
    mv Run_offline_DQM_1_cfg.py $runningdir/Run_offline_DQM_${i}_cfg.py

    rm filelist.txt
  done
}

wait_for_staging(){
  echo -e "Files are on castor.\nStaging =====>"
  get_done=0
  need_to_wait=1
  while [ $need_to_wait -eq 1 ];do
    need_to_wait=0
    for i in `seq 0 39`;do
      file=GainCalibration_${i}_$run.$ext
      if [ `is_file_present $indir/$file` -eq 0 ];then echo "File $file is not present in $indir ...";continue;fi	 
      stager_qry -M $indir/$file
      if [ `is_staged $indir/$file` -eq 0 ];then
        need_to_wait=1
	if [ $get_done -eq 1 ] ; then break ; fi
      fi
    done
    get_done=1
    if [ $need_to_wait -eq 1 ];then
      echo "At least one file is not staged. Will sleep for 5min before trying again ..."
      sleep 300
    fi
  done
  echo -e "Staging is finished !"
}


do_seds(){
  #Doing sed stuff
echo
  #sed "s#STOREdir#${storedir}#" < makeone.c > ${runningdir}/TEXToutput/makeone.c

}

submit_calib(){

  set_specifics $indir
    
  cat submit_template.sh |\
    sed "s#INDIR#$indir#" |\
    sed "s/RUN/$run/" |\
    sed "s/.EXT/.$ext/" |\
    sed "s#STOREDIR#${storedir}#" |\
    sed "s#T2_CP#${T2_CP}#"   |\
    sed "s#T2_TMP_DIR#${T2_TMP_DIR}#"   |\
    sed "s#CFGDIR#${runningdir}#" > ${runningdir}/submit_template.sh

  cd $runningdir
  
  for i in `seq 0 39`;do
    cfg=${runningdir}/Run_offline_DQM_${i}_cfg.py

    make_dir ${runningdir}/JOB_${i}
    
    rm -f submit_${i}.sh
    cat submit_template.sh |\
      sed "s/NUM/${i}/"    > submit_${i}.sh

    #qsub -q localgrid@cream01 -j oe -N job_${i} -o ${runningdir}/JOB_${i}/stdout submit_${i}.sh
    submit_to_queue ${run}_${i} ${runningdir}/JOB_${i}/stdout submit_${i}.sh
     
  done
}

resubmit_job(){
  set_specifics $storedir 
  
  cd $runningdir
  if [ `is_file_present $storedir/$ijob.root` -eq 1 ];then
    echo "Output of job $ijob is already in $storedir."
    exit
  fi
  
  echo "Re-submitting job $ijob:"
  submit_to_queue ${run}_${ijob} ${runningdir}/JOB_${ijob}/stdout submit_${ijob}.sh
  
}

stage_all_files(){
  set_specifics $indir
  f_to_stage=''
  if [ `is_on_castor $indir` -eq 1 ];then
    for file in `nsls $indir`;do
      f_to_stage="$f_to_stage $indir/$file"
    done
  fi
  if [ `is_on_castor $storedir` -eq 1 ];then
    for file in `nsls $storedir`;do
      f_to_stage="$f_to_stage $storedir/$file"
    done
  fi
  
  if [ `is_on_castor $indir` -eq 1 ] || [ `is_on_castor $storedir` -eq 1 ];then
    stage_list_of_files $f_to_stage
  else
    echo "Nothing to stage ..."
  fi
}


submit_hadd(){
  set_specifics ${storedir}
  rm -f ${runningdir}/hadd.sh
  cat hadd_template.sh |\
    sed "s#STOREDIR#${storedir}#" |\
    sed "s#T2_TMP_DIR#$T2_TMP_DIR#" |\
    sed "s#T2_CP#$T2_CP#" |\
    sed "s#T2_RM#$T2_RM#" |\
    sed "s#RUNNINGDIR#${runningdir}#"  > ${runningdir}/hadd.sh
    
  cd ${runningdir}
  rm -f hadd.txt
  touch temp.txt
  $T2_RM ${storedir}/GainCalibration.root
  stage_dir ${storedir}
  for file in `$T2_LS ${storedir}`;do
    echo `file_loc ${storedir}/$file` >> temp.txt
  done
  echo `cat temp.txt` > hadd.txt
  cat temp.txt
  rm -f temp.txt

  echo "Submitting hadd job to batch ..."
  submit_to_queue "${run}_hadd" `pwd`/hadd_stdout hadd.sh
  #bsub -q cmscaf1nh -J job_1 < Hadd.csh

}

submit_summary(){

  set_specifics ${storedir}
  if [ `$T2_LS  $storedir/GainCalibration.root 2>&1|grep "No such"|wc -l` -eq 1 ]; then
    echo "File $storedir/GainCalibration.root is not present ..."; exit ; fi ;
  stage_list_of_files $storedir/GainCalibration.root

  #making directories
  make_dir ${runningdir}/TEXToutput
  make_dir $runningdir/Summary_Run$run
  
  cp -fr scripts/make_SummaryPlots.cc ${runningdir}/Summary_Run$run/make_SummaryPlots_template.cc
  cp -fr scripts/gain_summary.txt  ${runningdir}/Summary_Run$run/gain_summary_template.tex
  cp -fr scripts/TMean.* ${runningdir}/Summary_Run$run/.
  cp -fr scripts/PixelNameTranslator.* ${runningdir}/Summary_Run$run/.
  
  cd $runningdir/Summary_Run$run

  rm -fr make_SummaryPlots.cc gain_summary.tex
  sed "s#RUNNUMBER#$run#" < make_SummaryPlots_template.cc > make_SummaryPlots.cc
  sed "s#RUNNUMBER#$run#" < gain_summary_template.tex > gain_summary.tex
  #rm -fr gain_summary_template.tex make_SummaryPlots_template.cc

  #rm -fr $T2_TMP_DIR/*
  #$T2_CP `file_loc $storedir/GainCalibration.root` $T2_TMP_DIR/GainCalibration.root
  echo "(root -l -b -x make_SummaryPlots.cc+\"(\"`file_loc $storedir/GainCalibration.root`\")\" -q)"
  root -l -b -x make_SummaryPlots.cc+"(\"`file_loc $storedir/GainCalibration.root`\")" -q

  rm -f gain_summary_final_run_$run.tex
  sed '/TOREPLACE/,$ d' < gain_summary.tex > gain_summary_final_run_$run.tex
  cat texSummary_Run${run}.tex >> gain_summary_final_run_$run.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_run_$run.tex

  pdflatex gain_summary_final_run_$run.tex
  pdflatex gain_summary_final_run_$run.tex
}

submit_summary_new(){

  set_specifics ${storedir}
  if [ `$T2_LS  $storedir/GainCalibration.root 2>&1|grep "No such"|wc -l` -eq 1 ]; then
    echo "File $storedir/GainCalibration.root is not present ..."; exit ; fi ;
  stage_list_of_files $storedir/GainCalibration.root

  #making directories
  make_dir ${runningdir}/TEXToutput
  make_dir $runningdir/Summary_Run$run
  
  cp -fr scripts/make_ComparisonPlots.cc ${runningdir}/Summary_Run$run/make_ComparisonPlots.cc
  cp -fr scripts/gain_summary.txt  ${runningdir}/Summary_Run$run/gain_summary_template.tex
  cp -fr scripts/TMean.* ${runningdir}/Summary_Run$run/.
  cp -fr scripts/PixelNameTranslator.* ${runningdir}/Summary_Run$run/.
  
  cd $runningdir/Summary_Run$run

  rm -fr gain_summary.tex
  #sed "s#RUNNUMBER#$run#" < make_SummaryPlots_template.cc > make_SummaryPlots.cc
  sed "s#RUNNUMBER#$run#" < gain_summary_template.tex > gain_summary.tex
  #rm -fr gain_summary_template.tex make_SummaryPlots_template.cc

  #rm -fr $T2_TMP_DIR/*
  #$T2_CP `file_loc $storedir/GainCalibration.root` $T2_TMP_DIR/GainCalibration.root
  echo "(root -l -b -x make_ComparisonPlots.cc+\"(\"`file_loc $storedir/GainCalibration.root`\",\"$run\")\" -q)"
  root -l -b -x make_ComparisonPlots.cc+"(\"`file_loc $storedir/GainCalibration.root`\",\"$run\")" -q

  rm -f gain_summary_final_run_$run.tex
  sed '/TOREPLACE/,$ d' < gain_summary.tex > gain_summary_final_run_$run.tex
  cat texSummary_Run${run}.tex >> gain_summary_final_run_$run.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_run_$run.tex

  pdflatex gain_summary_final_run_$run.tex
  pdflatex gain_summary_final_run_$run.tex
  echo -e "\nPDF file:\n `pwd`/gain_summary_final_$run.pdf"
}


compile_pdf(){

  cp -fr scripts/gain_summary.txt  ${runningdir}/Summary_Run$run/gain_summary_template.tex
  cd $runningdir/Summary_Run$run
  
  rm -fr gain_summary.tex
  sed "s#RUNNUMBER#$run#" < gain_summary_template.tex > gain_summary.tex

  rm -f gain_summary_final_run_$run.tex
  sed '/TOREPLACE/,$ d' < gain_summary.tex > gain_summary_final_run_$run.tex
  cat texSummary_Run${run}.tex >> gain_summary_final_run_$run.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_run_$run.tex

  pdflatex gain_summary_final_run_$run.tex
  pdflatex gain_summary_final_run_$run.tex
  echo -e "\nPDF file:\n `pwd`/gain_summary_final_$run.pdf"
}

set_files_for_comparison(){
  if [ "$run2" == "" ] && [ "$file2" == "" ];then
    run2=$file1
  
    run=$run1
    read_config
    file1=$storedir/GainCalibration.root
  
    run=$run2
    read_config
    file2=$storedir/GainCalibration.root  
  fi
  
  echo -e "\n-----------------------------------------------------------"
  echo "Comparing run $run1 & run $run2"
  echo "File for run $run1: $file1"
  echo "File for run $run2: $file2"
  echo -e "-----------------------------------------------------------\n"
  
  
}


compare_runs(){

  #echo t $run1 t $file1 t $run2 t $file2

  if [ "$run1" == "0" ] || [ "$run2" == "0" ] || [ "$file1" == "" ] || [ "$file2" == "" ];then usage ; fi

  stage_list_of_files $file1 $file2
  
  dir=Comp_${run1}-${run2}
  make_dir $dir

  set_specifics $file1
  if [ `$T2_LS $file1 2>&1|grep "No such"|wc -l` -eq 1 ]; then
    echo "File $file1 is not present ..."; exit ; fi ;
  file1=${T2_FSYS}${file1}

  set_specifics $file2
  if [ `$T2_LS $file2 2>&1|grep "No such"|wc -l` -eq 1 ]; then
    echo "File $file2 is not present ..."; exit ; fi ;
  file2=${T2_FSYS}${file2}

  #cat scripts/make_ComparisonPlots.cc |\
  #  sed "s#RUNNUMBER1#${run1}#" |\
  #  sed "s#RUNNUMBER2#${run2}#" > $dir/make_ComparisonPlots.cc
  
  cp -fr scripts/make_ComparisonPlots.cc $dir/.
  cp -fr scripts/TMean.* $dir/.
  cp -fr scripts/PixelNameTranslator.* $dir/.
  
  cd $dir
  
  echo "( root -l -b -q make_ComparisonPlots.cc+\"(\"$file1\",\"$run1\",\"$file2\",\"$run2\")\" )"
  root -l -b -q make_ComparisonPlots.cc+"(\"$file1\",\"$run1\",\"$file2\",\"$run2\")"
  
  
  cp -fr ../scripts/gain_summary.txt  gain_summary_template.tex
  rm -fr gain_summary.tex
  sed "s#RUNNUMBER#$run1-$run2#" < gain_summary_template.tex > gain_summary.tex

  rm -f gain_summary_final_$run1-$run2.tex
  sed '/TOREPLACE/,$ d' < gain_summary.tex > gain_summary_final_$run1-$run2.tex
  cat texSummary_Run${run1}-$run2.tex >> gain_summary_final_$run1-$run2.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_$run1-$run2.tex

  pdflatex gain_summary_final_$run1-$run2.tex &>  latex.log
  pdflatex gain_summary_final_$run1-$run2.tex >> latex.log 2>&1
  echo -e "\nPDF file:\n `pwd`/gain_summary_final_$run1-$run2.pdf"
}


make_payload(){

  #check if CondTools is present
  if [ ! -d ../../../CondTools/SiPixel/test ];then
    echo "You need to check-out the CondTools/SiPixel package:"
    echo "cvs co -r MY_CMSSW_VERSION CondTools/SiPixel"
    exit
  fi


  cd ../../../CondTools/SiPixel/test
  if [ ! -f SiPixelGainCalibrationReadDQMFile_cfg.py ];then
    echo "You are missing SiPixelGainCalibrationReadDQMFile_cfg.py. Please fix this !!" ; exit
  fi
  cp SiPixelGainCalibrationReadDQMFile_cfg.py original_SiPixelGainCalibrationReadDQMFile_cfg.py
  
  
  if [ ! -f prova.db ];then
    echo "prova.db not found ... That is not normal. Re-checkout CondTools/SiPixel ..." ; exit
  fi
  
  set_specifics $storedir
  stage_list_of_files $storedir/GainCalibration.root
  
  file=$T2_TMP_DIR/GainCalibration_$run.root
  $T2_CP $storedir/GainCalibration.root $file
  #set_specifics $file
  
  ###########################################   OFFLINE PAYLOAD
  
  
  payload=prova_GainRun${run}_31X.db
  cp prova.db $payload
  payload_root=Summary_payload_Run${run}.root
  
  echo "RM: `$T2_RM $storedir/$payload`"
  echo "RM: `$T2_RM $storedir/$payload_root`"
  
  #Changing some parameters in the python file:
  cat original_SiPixelGainCalibrationReadDQMFile_cfg.py |\
    sed "s#file:///tmp/rougny/test.root#`file_loc $file`#"  |\
    sed 's#useMeanWhenEmpty = cms.untracked.bool(False)#useMeanWhenEmpty = cms.untracked.bool(True)#'|\
    sed "s#sqlite_file:prova.db#sqlite_file:$T2_TMP_DIR/${payload}#" |\
    sed "s#/tmp/rougny/histos.root#$T2_TMP_DIR/$payload_root#" > SiPixelGainCalibrationReadDQMFile_cfg.py
    
    
  #cat SiPixelGainCalibrationReadDQMFile_cfg.py
  
  echo -e "\n--------------------------------------"
  echo "Making the payload for offline:"
  echo "  $storedir/$payload"
  echo "  ==> Summary root file: $payload_root"
  echo -e "--------------------------------------\n"
  
  echo "  (\" cmsRun SiPixelGainCalibrationReadDQMFile_cfg.py \" )"
  cmsRun SiPixelGainCalibrationReadDQMFile_cfg.py
  
  $T2_CP $T2_TMP_DIR/$payload $storedir/$payload
  $T2_CP $T2_TMP_DIR/$payload_root $storedir/$payload_root
   
  rm -f $T2_TMP_DIR/${payload}
  rm -f $T2_TMP_DIR/${payload_root}
  
  
  
  ###########################################   HLT PAYLOAD
 
  payload=prova_GainRun${run}_31X_HLT.db
  cp prova.db $payload
  payload_root=Summary_payload_Run${run}_HLT.root
  
  echo "RM: `$T2_RM $storedir/$payload`"
  echo "RM: `$T2_RM $storedir/$payload_root`"
  
  
  #Changing some parameters in the python file:
  cat original_SiPixelGainCalibrationReadDQMFile_cfg.py |\
    sed "s#file:///tmp/rougny/test.root#`file_loc $file`#"  |\
    sed 's#useMeanWhenEmpty = cms.untracked.bool(False)#useMeanWhenEmpty = cms.untracked.bool(True)#'|\
    sed "s#sqlite_file:prova.db#sqlite_file:$T2_TMP_DIR/${payload}#" |\
    sed "s#/tmp/rougny/histos.root#$T2_TMP_DIR/$payload_root#"|\
    sed "s#cms.Path(process.readfileOffline)#cms.Path(process.readfileHLT)#"|\
    sed "s#record = cms.string('SiPixelGainCalibrationOfflineRcd')#record = cms.string('SiPixelGainCalibrationForHLTRcd')#"|\
    sed "s#GainCalib_TEST_offline#GainCalib_TEST_hlt#" > SiPixelGainCalibrationReadDQMFile_cfg.py
    
    
  #cat SiPixelGainCalibrationReadDQMFile_cfg.py
  
  echo -e "\n--------------------------------------"
  echo "Making the payload for offline:"
  echo "  $storedir/$payload"
  echo "  ==> Summary root file: $payload_root"
  echo -e "--------------------------------------\n"
  
  echo "  (\" cmsRun SiPixelGainCalibrationReadDQMFile_cfg.py \" )"
  cmsRun SiPixelGainCalibrationReadDQMFile_cfg.py
  
  $T2_CP $T2_TMP_DIR/$payload $storedir/$payload
  $T2_CP $T2_TMP_DIR/$payload_root $storedir/$payload_root
  
  
  
  ####################################
   
  rm -f $T2_TMP_DIR/${payload}
  rm -f $T2_TMP_DIR/${payload_root}
  rm -f $file  
  
  cp original_SiPixelGainCalibrationReadDQMFile_cfg.py SiPixelGainCalibrationReadDQMFile_cfg.py
}




create=0
submit=0
resubmit=0
stage=0
hadd=0
summary=0
pdf=0
compare=0
verbose=0
ijob=-1
prova=0

run1=0
file1=""
run2=0
file2=""

#lock

if [ $# -eq 0 ] ; then usage ; fi
for arg in $* ; do
  case $arg in
    -create)       create=1   ; run=$2 ; indir=$3 ; storedir=$4 ; shift ; shift ; shift ; shift ;;
    -submit)       submit=1   ; run=$2 ; shift ;;
    -resubmit)     resubmit=1 ; run=$2 ; ijob=$3  ; shift ;;
    -stage)       stage=1    ; run=$2 ; shift ;;
    -hadd)         hadd=1     ; run=$2 ; shift ;;
    -summary)      summary=1  ; run=$2 ; shift ;;
    -pdf)          pdf=1      ; run=$2 ; shift ;;
    -compare)      compare=1  ; run=0  ; run1=$2 ; file1=$3 ; run2=$4 ; file2=$5 ; shift ; shift ; shift ; shift ; shift ;;
    -payload)      prova=1  ; run=$2 ; shift ;;
    -v)            verbose=1  ; shift ;;
    -help)         usage ;;
    *)             ;;
  esac
done

if [ "$run" == "" ] ; then usage ; fi

if [ $create -eq 1 ];then
  set_parameters
  create
  write_config
  make_file_list
fi

if [ $submit -eq 1 ];then
  read_config
  submit_calib
fi

if [ $resubmit -eq 1 ];then
  read_config
  resubmit_job
fi

if [ $stage -eq 1 ];then
  read_config
  stage_all_files
fi

if [ $hadd -eq 1 ];then
  read_config
  submit_hadd
fi

if [ $summary -eq 1 ];then
  read_config
  submit_summary_new
fi


if [ $pdf -eq 1 ];then
  read_config
  compile_pdf
fi


if [ $compare -eq 1 ];then
  set_files_for_comparison
  compare_runs
fi


if [ $prova -eq 1 ];then
  read_config
  make_payload
fi

rm -f .lock_gaincalib

