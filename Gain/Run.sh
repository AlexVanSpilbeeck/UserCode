#!/bin/bash

########################
##
##   Script to run the gain calibration on all 40 FEDs
##
##
##   HowTo : source Run_GainCalibration.sh RUNNUMBER INPUTDIR STOREDIR
##
##
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
./Run.sh -create  RUNNUMBER INPUTDIR STOREDIR  : will create the needed directories, python files to run the calib.
./Run.sh -submit  RUNNUMBER : will launch the 40 calibration jobs
./Run.sh -hadd    RUNNUMBER : will hadd all 40 output files of the calib jobs into one file
./Run.sh -summary RUNNUMBER : will launch the summary job.
./Run.sh -pdf     RUNNUMBER : will recompile the latex file to recreate the pdf summary.
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
  echo "run : $run"
  echo "indir : $indir"
  echo "storedir : $storedir"
  runningdir=`pwd`/Run_$run
}

make_dir(){
  if [ ! -d $1 ] ; then mkdir $1
  else rm -r $1/*
  fi
}


create(){
  
  #making running directory
  make_dir ${runningdir}

  #cleaning
  rm -f filelist.txt client_template_calib_cfg.py es.log
  
  #cleaning output dir
  set_specifics ${storedir}
  $T2_RM -r ${storedir}
  $T2_MKDIR ${storedir}
  $T2_CHMOD 0777 ${storedir}/ #TOFIX ?
  
  #chmod -R 0777 $indir

  if [ `is_on_castor $indir` -eq 1 ] ; then wait_for_staging ; fi
}


make_file_list(){

  #Copying template specific to gain calib to the general one used by Run_offline_DQM.csh
  cp -f gaincalib_template_cfg.py client_template_calib_cfg.py

  #if [ `is_on_castor $indir` -eq 1 ] ; then wait_for_staging ; fi

  #making python files
  touch filelist.txt
  for i in `seq 0 39`;do
    file=GainCalibration_${i}_$run.$ext
    if [ `$T2_LS $indir/$file 2>&1|grep "No such"|wc -l` -eq 1 ];then echo "File $file is not present in $indir ...";continue;fi 

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
      file=$indir/GainCalibration_${i}_$run.$ext
      if [ `$T2_LS $indir/$file 2>&1|grep "No such"|wc -l` -eq 1 ];then echo "File $file is not present in $indir ...";continue;fi	 
      stager_qry -M $indir/GainCalibration_${i}_$run.$ext
      if [ `is_staged $file` -eq 0 ];then
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
    submit_to_queue job_${i} ${runningdir}/JOB_${i}/stdout submit_${i}.sh
     
  done
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
  for file in `$T2_LS ${storedir}`;do
    echo `file_loc ${storedir}/$file` >> temp.txt
  done
  echo `cat temp.txt` > hadd.txt
  cat temp.txt
  rm -f temp.txt

  echo "Submitting hadd job to batch ..."
  submit_to_queue "hadd_job" `pwd`/hadd_stdout hadd.sh
  #bsub -q cmscaf1nh -J job_1 < Hadd.csh

}

submit_summary(){

  set_specifics ${storedir}

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
  cat Summary_Run${run}number.txt >> gain_summary_final_run_$run.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_run_$run.tex

  pdflatex gain_summary_final_run_$run.tex
  pdflatex gain_summary_final_run_$run.tex
}


compile_pdf(){

  cp -fr scripts/gain_summary.txt  ${runningdir}/Summary_Run$run/gain_summary_template.tex
  cd $runningdir/Summary_Run$run
  
  rm -fr gain_summary.tex
  sed "s#RUNNUMBER#$run#" < gain_summary_template.tex > gain_summary.tex

  rm -f gain_summary_final_run_$run.tex
  sed '/TOREPLACE/,$ d' < gain_summary.tex > gain_summary_final_run_$run.tex
  cat Summary_Run${run}number.txt >> gain_summary_final_run_$run.tex
  sed '1,/TOREPLACE/d'< gain_summary.tex >> gain_summary_final_run_$run.tex

  pdflatex gain_summary_final_run_$run.tex
  pdflatex gain_summary_final_run_$run.tex
}



create=0
submit=0
hadd=0
summary=0
pdf=0
verbose=0

if [ $# -eq 0 ] ; then usage ; fi
for arg in $* ; do
  case $arg in
    -create)       create=1   ; run=$2 ; indir=$3 ; storedir=$4 ; shift ; shift ; shift ; shift ;;
    -submit)       submit=1   ; run=$2 ; shift ;;
    -hadd)         hadd=1     ; run=$2 ; shift ;;
    -summary)      summary=1  ; run=$2 ; shift ;;
    -pdf)          pdf=1      ; run=$2 ; shift ;;
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

if [ $hadd -eq 1 ];then
  read_config
  submit_hadd
fi

if [ $summary -eq 1 ];then
  read_config
  submit_summary
fi


if [ $pdf -eq 1 ];then
  read_config
  compile_pdf
fi


