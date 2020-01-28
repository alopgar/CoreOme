#!/bin/bash
#SBATCH -t 00:20:00
#SBATCH -J 2.2_RndmFrst_%a
#SBATCH --mem=0 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --ntasks=20
#SBATCH --array=1-20
#SBATCH --output=2.2_RndmFrst_%a.out
#SBATCH --error=2.2_RndmFrst_%a.out

## MODULE LOADING AND SCRIPT RUNNING
module load cesga/2018 gcc/6.4.0 R/3.5.1

## 2.2) PREVALENCE FILTER: RandomForest
SUBDIR=$WDIR/Filtered/RFSubsets
if [ ! -d $SUBDIR ]; then mkdir $SUBDIR; fi

subset=$(ls $SUBDIR | egrep "RanFor_thr[0-9]{2}\.txt" | sed "s|.txt||" | sed -n ${SLURM_ARRAY_TASK_ID}p)
Rscript $BINPATH/Filter3_RandomForest.R $paramfile $SUBDIR/$subset.txt
