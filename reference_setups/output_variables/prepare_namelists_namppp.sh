#!/bin/bash

# This script is called as
# prepare_namelist_namppp.sh input_file
#
# The input_file must be organised with at least five columns as in the example below.
# Any extra columns are ignored.
# Any line starting with ! is ignored (can be used for comments).
# It is checked that the grib code and the short variable name exist in RCA_var_file (e.g. rca35name.txt).
# Output frequency is specified with one of dd pp qq hh ss cl ml,
# where ml represents model-level files (fc-files without any extension).
# The output file, namelists_namppp.dat, should be placed in your WORK_DIR or in
# the directory where you submit your run_script.
#
# 250 105 1	t2m_i	qq
# 15 105 2	t2max	dd
# 16 105 2	t2min	dd
# 33 109 40	u40	ml
#

RCA_var_file=''
if [[ -f $RCA_var_file ]]; then
  nvarrca=`wc -l $RCA_var_file | awk '{print $1}'`
  RCA_var_file_exist=true
else
  echo "  "
  echo "********************************************************************"
  echo " Warning: No valid RCA_var_file is specified!!"
  echo " No validation of correctly specified GRIB codes will be performed."
  echo "********************************************************************"
  echo "  "
  RCA_var_file_exist=false
fi

noutvar=`wc -l $1 | awk '{print $1}'`

let "nndd=0"
let "nnpp=0"
let "nnqq=0"
let "nnhh=0"
let "nnss=0"
let "nncl=0"
let "nnml=0"

# Loop over variables in input list $1
for ((i=1;i<=noutvar;i++)); do
   text=(`sed -n "$i p" $1`)
   par=${text[0]}

   if [[ ${par:0:1} != '!' ]]; then

      type=${text[1]}
      lev=${text[2]}
      shortname=${text[3]}
      outfreq=${text[4]}

# Check each line in $1 against RCA_var_file
      if [[ $RCA_var_file_exist == "true" ]];then
         varfit=false
         for ((jj=1;jj<=nvarrca;jj++)); do
            textrca=(`sed -n "$jj p" $RCA_var_file`)
            parrca=${textrca[0]}
            if [[ $par == $parrca ]]; then
               levrca=${textrca[2]}
               if [[ $lev == $levrca ]]; then
                  typerca=${textrca[1]}
                  shortnamerca=${textrca[3]}
                  if [[ $type == $typerca && $shortname == $shortnamerca ]]; then
                     varfit=true
                  fi
               fi
            fi
         done
      else
         varfit=true
      fi
#
      if [[ $varfit == "false" || $outfreq != dd && $outfreq != pp && $outfreq != qq && $outfreq != hh && $outfreq != ss && $outfreq != cl && $outfreq != ml ]]; then
         echo "****************************"
         echo "Not a valid input line:"
         echo ${text[*]}
         echo "****************************"
      else
         if [ $type == 100 ]; then
            level=`expr $lev \* 100`
            echo $par" "      >> namppp_tmp_100_iwmomlp_$outfreq
#            echo $level". "   >> namppp_tmp_100_alevmlp_$outfreq
            echo $level" "   >> namppp_tmp_100_alevmlp_$outfreq
   
         else
            if [ $outfreq == dd ]; then
               nnnn=$nndd
            elif [ $outfreq == pp ]; then
               nnnn=$nnpp
            elif [ $outfreq == qq ]; then
               nnnn=$nnqq
            elif [ $outfreq == hh ]; then
               nnnn=$nnhh
            elif [ $outfreq == ss ]; then
               nnnn=$nnss
            elif [ $outfreq == cl ]; then
               nnnn=$nncl
            elif [ $outfreq == ml ]; then
               nnnn=$nnml
            fi
            nnnn=`expr $nnnn + 1`
            echo "iwmoslp("$nnnn")="$par     >> namppp_tmp_$outfreq
            echo "ltypslp("$nnnn")="$type    >> namppp_tmp_$outfreq
#            echo "alevslp("$nnnn")="$lev"."  >> namppp_tmp_$outfreq
            echo "alevslp("$nnnn")="$lev  >> namppp_tmp_$outfreq
            if [ $outfreq == dd ]; then
              nndd=$nnnn
            elif [ $outfreq == pp ]; then
              nnpp=$nnnn
            elif [ $outfreq == qq ]; then
              nnqq=$nnnn
            elif [ $outfreq == hh ]; then
              nnhh=$nnnn
            elif [ $outfreq == ss ]; then
              nnss=$nnnn
            elif [ $outfreq == cl ]; then
              nncl=$nnnn
            elif [ $outfreq == ml ]; then
              nnml=$nnnn
            fi
         fi
      fi
   fi
done

# Arrange the pressure level output parameters
for outfreq in dd pp qq hh ss; do
   if [ -f namppp_tmp_100_iwmomlp_$outfreq ]; then
      sort namppp_tmp_100_iwmomlp_$outfreq | uniq > test1
      nwmomlp=`wc -l test1 | awk '{print $1}'`
      echo -n "iwmomlp=" > namppp_tmp_100_iwmomlp_$outfreq
      for i in `cat test1`; do
         echo -n $i", " >> namppp_tmp_100_iwmomlp_$outfreq
      done
      echo " " >> namppp_tmp_100_iwmomlp_$outfreq
      sort namppp_tmp_100_alevmlp_$outfreq | uniq > test1
      nlevmlp=`wc -l test1 | awk '{print $1}'`
      echo -n "alevmlp=" > namppp_tmp_100_alevmlp_$outfreq
      for i in `cat test1`; do
         echo -n $i", " >> namppp_tmp_100_alevmlp_$outfreq
      done
      echo "nlevmlp="$nlevmlp >> namppp_tmp_100_$outfreq
      echo "ltypmlp=100" >> namppp_tmp_100_$outfreq
      cat namppp_tmp_100_alevmlp_$outfreq >> namppp_tmp_100_$outfreq
      echo " " >> namppp_tmp_100_$outfreq
      echo "nwmomlp="$nwmomlp >> namppp_tmp_100_$outfreq
      cat namppp_tmp_100_iwmomlp_$outfreq >> namppp_tmp_100_$outfreq
      rm namppp_tmp_100_iwmomlp_$outfreq namppp_tmp_100_alevmlp_$outfreq
      rm test1
   fi
done


echo "&nampp" > namelists_namppp.dat
nppstr=0
if [[ -f namppp_tmp_100_dd || -f namppp_tmp_dd ]]; then
  echo "&namppp" >> namppp_dd
  echo "lunppfp=802" >> namppp_dd
  echo "sufixp='dd'" >> namppp_dd
  if [ -f namppp_tmp_100_dd ]; then
    cat namppp_tmp_100_dd >> namppp_dd
  else
    echo "nlevmlp=0" >> namppp_dd
    echo "ltypmlp=100" >> namppp_dd
    echo "nwmomlp=0" >> namppp_dd
  fi
  if [ -f namppp_tmp_dd ]; then
    echo "nslp="$nndd >> namppp_dd
    cat namppp_tmp_dd >> namppp_dd
  else
    echo "nslp=0" >> namppp_dd
  fi
  echo "/" >> namppp_dd
  (( nppstr += 1 ))
  echo "suff("$nppstr")='dd'" >> namelists_namppp.dat
fi

if [[ -f namppp_tmp_100_pp || -f namppp_tmp_pp ]]; then
  echo "&namppp" >> namppp_pp
  echo "lunppfp=801" >> namppp_pp
  echo "sufixp='pp'" >> namppp_pp
  if [ -f namppp_tmp_100_pp ]; then
    cat namppp_tmp_100_pp >> namppp_pp
  else
    echo "nlevmlp=0" >> namppp_pp
    echo "ltypmlp=100" >> namppp_pp
    echo "nwmomlp=0" >> namppp_pp
  fi
  if [ -f namppp_tmp_pp ]; then
    echo "nslp="$nnpp >> namppp_pp
    cat namppp_tmp_pp >> namppp_pp
  else
    echo "nslp=0" >> namppp_pp
  fi
  echo "/" >> namppp_pp
  (( nppstr += 1 ))
  echo "suff("$nppstr")='pp'" >> namelists_namppp.dat
fi

if [[ -f namppp_tmp_100_qq || -f namppp_tmp_qq ]]; then
  echo "&namppp" >> namppp_qq
  echo "lunppfp=803" >> namppp_qq
  echo "sufixp='qq'" >> namppp_qq
  if [ -f namppp_tmp_100_qq ]; then
    cat namppp_tmp_100_qq >> namppp_qq
  else
    echo "nlevmlp=0" >> namppp_qq
    echo "ltypmlp=100" >> namppp_qq
    echo "nwmomlp=0" >> namppp_qq
  fi
  if [ -f namppp_tmp_qq ]; then
    echo "nslp="$nnqq >> namppp_qq
    cat namppp_tmp_qq >> namppp_qq
  else
    echo "nslp=0" >> namppp_qq
  fi
  echo "/" >> namppp_qq
  (( nppstr += 1 ))
  echo "suff("$nppstr")='qq'" >> namelists_namppp.dat
fi

if [[ -f namppp_tmp_100_hh || -f namppp_tmp_hh ]]; then
  echo "&namppp" >> namppp_hh
  echo "lunppfp=807" >> namppp_hh
  echo "sufixp='hh'" >> namppp_hh
  if [ -f namppp_tmp_100_hh ]; then
    cat namppp_tmp_100_hh >> namppp_hh
  else
    echo "nlevmlp=0" >> namppp_hh
    echo "ltypmlp=100" >> namppp_hh
    echo "nwmomlp=0" >> namppp_hh
  fi
  if [ -f namppp_tmp_hh ]; then
    echo "nslp="$nnhh >> namppp_hh
    cat namppp_tmp_hh >> namppp_hh
  else
    echo "nslp=0" >> namppp_hh
  fi
  echo "/" >> namppp_hh
  (( nppstr += 1 ))
  echo "suff("$nppstr")='hh'" >> namelists_namppp.dat
fi

if [[ -f namppp_tmp_100_ss || -f namppp_tmp_ss ]]; then
  echo "&namppp" >> namppp_ss
  echo "lunppfp=804" >> namppp_ss
  echo "sufixp='ss'" >> namppp_ss
  if [ -f namppp_tmp_100_ss ]; then
    cat namppp_tmp_100_ss >> namppp_ss
  else
    echo "nlevmlp=0" >> namppp_ss
    echo "ltypmlp=100" >> namppp_ss
    echo "nwmomlp=0" >> namppp_ss
  fi
  if [ -f namppp_tmp_ss ]; then
    echo "nslp="$nnss >> namppp_ss
    cat namppp_tmp_ss >> namppp_ss
  else
    echo "nslp=0" >> namppp_ss
  fi
  echo "/" >> namppp_ss
  (( nppstr += 1 ))
  echo "suff("$nppstr")='ss'" >> namelists_namppp.dat
fi

if [ -f namppp_tmp_cl ]; then
  echo "&namppp" >> namppp_cl
  echo "lunppfp=805" >> namppp_cl
  echo "prefixp='cl'" >> namppp_cl
  echo "sufixp='gb'" >> namppp_cl
  echo "ltypmlp=-1" >> namppp_cl
  echo "nlevmlp=0" >> namppp_cl
  echo "nslp="$nncl >> namppp_cl
  cat namppp_tmp_cl >> namppp_cl
  echo "/" >> namppp_cl
  (( nppstr += 1 ))
  echo "suff("$nppstr")='gb'" >> namelists_namppp.dat
fi

if [ -f namppp_tmp_ml ]; then
  echo "&namppp" >> namppp_ml
  echo "lunppfp=806" >> namppp_ml
  echo "prefixp='fc'" >> namppp_ml
  echo "sufixp='  '" >> namppp_ml
  echo "ltypmlp=100" >> namppp_ml
  echo "nlevmlp=0" >> namppp_ml
  echo "nwmomlp=0" >> namppp_ml
  echo "nslp="$nnml >> namppp_ml
  cat namppp_tmp_ml >> namppp_ml
  echo "/" >> namppp_ml
  (( nppstr += 1 ))
  echo "suff("$nppstr")='  '" >> namelists_namppp.dat
fi
echo "nppstr="$nppstr >> namelists_namppp.dat
echo "/" >> namelists_namppp.dat

rm *_tmp_*
cat namppp_* >> namelists_namppp.dat
rm namppp_*

