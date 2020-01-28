#!/bin/bash
#SBATCH -t 00:01:00
#SBATCH -J SQMr_filter 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/SQMr_exe.out
#SBATCH --error=/mnt/netapp2/Home_FT2/home/otras/ini/alg/SQMr_exe.out

# For no sbatch output, put '--ouptut=/dev/null'.
export BINPATH=~/bin/MinION_SqueezeMeta
export WDIR=/mnt/lustre/scratch/home/otras/ini/alg/Results/Metagenome/METALGEN_SQMReads
export SUBDIR=$WDIR/Filtered/RFSubsets

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
job11="$BINPATH.job"
echo "#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -J 1.1_Samples
#SBATCH --mem=10GB 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --ntasks=1
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/1.1_Samples.out
#SBATCH --error=/mnt/netapp2/Home_FT2/home/otras/ini/alg/1.1_Samples.out
Rscript $BINPATH/Filter1_Samples.R $paramfile" > $job11

jid11=$(sbatch $job11)
jid11=$(grep -Eo "[0-9]*"<<< $jid11) 

## 1.2) TAXONOMY FILTER:
job12="$BINPATH.job"
echo "#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -J 1.2_Taxonomy
#SBATCH --mem=2GB 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --ntasks=1
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/1.2_Taxonomy.out
#SBATCH --error=/mnt/netapp2/Home_FT2/home/otras/ini/alg/1.2_Taxonomy.out
Rscript $BINPATH/Filter2_Taxa.R $paramfile" > $job12

if [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = taxonomy ]; then
    jid12=$(sbatch --dependency=afterok:$jid11 $job12)
    jid12=$(grep -Eo "[0-9]*"<<< $jid12)
elif [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = keggs ]; then
    echo "datatype = keggs -> No taxonomy filter done"
else
    echo "ERROR: Rparams datatype value not valid. Please choose: 'taxonomy' or 'keggs'"
    exit 1 
fi

## 2.1) PREVALENCE FILTER: Making data subsets
job21="$BINPATH.job"
echo "#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -J 2.1_Subsets
#SBATCH --mem=0 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/2.1_Subsets.out
#SBATCH --error=/mnt/netapp2/Home_FT2/home/otras/ini/alg/2.1_Subsets.out
Rscript $BINPATH/Filter3_Subsets.R $paramfile $phenotmod" > $job21

if [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = taxonomy ]; then
    jid21=$(sbatch --dependency=afterok:$jid12 $job21)
elif [ $(grep "datatype <-" $paramfile | awk -F'"' '{print $2}') = keggs ]; then
    jid21=$(sbatch --dependency=afterok:$jid11 $job21)
fi
jid21=$(grep -Eo "[0-9]*"<<< $jid21)

## 2.2) PREVALENCE FILTER: RandomForest
jid22=$(sbatch --dependency=afterok:$jid21 $BINPATH/SQMr_randomforest.sh)
jid22=$(grep -Eo "[0-9]*"<<< $jid22)

### 2.3) PREVALENCE FILTER: OOB table & VIP from filtered features
jid23=$(sbatch --dependency=afterok:$jid22 $BINPATH/SQMr_OOB_VIP.sh)
jid23=$(grep -Eo "[0-9]*"<<< $jid23)

cd $WDIR/Filtered
tar -czvf SQMr_filters_taxonomy.tar.gz SQMr_filtered_genera.csv SQMr_important_Lprev_genera.csv R_report.txt OOB_table.tsv summary_stats.tsv
tar -czvf SQMr_filters_keggs.tar.gz SQMr_filtered_kegg.csv SQMr_important_Lprev_kegg.csv R_report.txt OOB_table.tsv summary_stats.tsv
