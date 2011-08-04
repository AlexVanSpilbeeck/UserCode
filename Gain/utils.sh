#################
#Queues
q_lxplus=cmscaf1nd
q_T2=localgrid@cream01


#order: name stdout submit.sh
submit_to_queue(){
 if [ `working_on_lxplus` -eq 1 ]; then
   bsub -q $q_lxplus -J $1 -eo $2 < $3
 else
   qsub -q $q_T2 -j oe -N $1 -o $2 $3
 fi
}


##############################################
# Your T2 Specifics (if not on castor/lxplus)
set_T2_specifics(){
  T2_LS='ls'
  T2_RM='rm'
  T2_CP='dccp -d 2'
  T2_MKDIR='mkdir'
  T2_CHMOD='chmod'
  T2_FSYS='dcap://'
  T2_CHECK='/pnfs/iihe/'
  T2_TMP_DIR='/scratch/'
  if [ $verbose -eq 1 ] ; then echo "Setting T2 specifics ..." ; fi
}

set_castor_specifics(){
  T2_LS='nsls'
  T2_RM='nsrm'
  T2_CP='rfcp'
  T2_MKDIR='rfmkdir'
  T2_CHMOD='rfchmod'
  T2_FSYS='rfio://'
  T2_TMP_DIR='/tmp/$USER/'
  if [ $verbose -eq 1 ] ; then echo "Setting castor specifics ..." ; fi
}

set_local_specifics(){
  T2_LS='ls'
  T2_RM='rm'
  T2_CP='cp'
  T2_MKDIR='mkdir'
  T2_CHMOD='chmod'
  T2_FSYS=''
  T2_TMP_DIR='/scratch/'
  if [ $verbose -eq 1 ] ; then echo "Setting local specifics ..." ; fi
}
set_specifics(){
  if [ $verbose -eq 1 ] ; then echo "Specifics for $1" ; fi
  if      [ `is_on_castor $1` -eq 1 ]; then set_castor_specifics ;
  else if [ `is_on_T2 $1` -eq 1 ]; then set_T2_specifics;
  else  set_local_specifics ;
  fi ; fi ;
}


is_on_castor(){
  if [ `echo $1 | grep -e "/castor/cern.ch/"|wc -l` -eq 1 ];then echo 1
  else echo 0
  fi
}

is_on_T2(){
  if [ `echo $1 | grep -e "$T2_CHECK"|wc -l` -eq 1 ];then echo 1
  else echo 0
  fi
}

file_loc(){
  file=$1
  if [ `is_on_castor $file` -eq 1 ];then
    echo "rfio://$file"
  else
    if [ `is_on_T2 $file` -eq 1 ];then
      echo "${T2_FSYS}$file"
    else
      echo $file
    fi
  fi
}

working_on_lxplus(){
  if [ `uname -a | grep lxplus | wc -l` -gt 0 ];then echo 1
  else echo 0
  fi
}

is_staged(){
  status=`stager_qry -M $1`
  if [ `echo $status|grep -e "STAGED"|wc -l` -eq 1 ] || [ `echo $status|grep -e "CANBEMIGR"|wc -l` -eq 1 ];then
    echo 1
  else
    if [ `echo $status|grep -e "STAGEIN"|wc -l` -eq 1 ];then
      echo 0
    else
      stager_get -M $1 &> /dev/null
      echo 0
    fi
  fi
}

