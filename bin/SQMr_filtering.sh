#!/bin/bash
#SBATCH -t 00:50:00
#SBATCH -J SQMr_filter
#SBATCH --mem=0 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/SQMr_filt_report.out

## INITIAL VARIABLES
export BINPATH=~/bin/MinION_SqueezeMeta
export WDIR=/mnt/lustre/scratch/home/otras/ini/alg/Results/Metagenome/METALGEN_SQMReads

duplicated=false

export paramfile=$BINPATH/Rparams.R
export phenotmod=$BINPATH/Rphenotmods.R

## MODULE LOADING AND SCRIPT RUNNING
module load cesga/2018 gcc/6.4.0 R/3.5.1 parallel jdk/8u181
if [ ! -d $WDIR/Filtered ]; then mkdir $WDIR/Filtered; fi

## Replace name of duplicated column names adding '_old', and remove '_new' from new duplicates.
if [ "$duplicated" = true ]; then
  for i in $WDIR/SQMreads/*combined*; do
    echo "Processing duplicated columns in $i"
    perl -i~ -pale'
    next if $. != 1;
    my %map = map { /^(.*)_new\z/s ? ( $1 => "$1_old", $_ => $1 ) : () } @F;
    @F = map { $map{$_} // $_ } @F;
    $_ = "@F";
    s/\s/\t/g; s/^/\t/;
    ' $i
    sed -r -i 's/SH-19/SH-19_old/g' $i
    sed -r -i 's/SH-27/SH-27_old/g' $i
    sed -r -i 's/SH-32/SH-32_old/g' $i
  done
  if [ ! -d $WDIR/SQMreads/Backups ]; then mkdir $WDIR/SQMreads/Backups; fi
  mv $WDIR/SQMreads/*~ $WDIR/SQMreads/Backups/
elif [ "$duplicated" = false ]; then
  echo "Correction of duplicated column names not selected"
else
  echo "ERROR: duplicated value not valid. Please choose booleans: true or false"
  exit 1
fi

## 1.1) SAMPLES FILTER:
Rscript $BINPATH/Filter1_Samples.R $paramfile

## 1.2) TAXONOMY FILTER:
if [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = taxonomy ]; then
  Rscript $BINPATH/Filter2_Taxa.R $paramfile
elif [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = keggs ]; then
  echo "datatype = keggs -> No taxonomy filter done"
else
  echo "ERROR: Rparams datatype value not valid. Please choose: 'taxonomy' or 'keggs'"
  exit 1 
fi

## 2.1) PREVALENCE FILTER: Making data subsets
Rscript $BINPATH/Filter3_Subsets.R $paramfile $phenotmod

## 2.2) PREVALENCE FILTER: RandomForest
jid1=$(sbatch $BINPATH/SQMr_randomforest.sh)
