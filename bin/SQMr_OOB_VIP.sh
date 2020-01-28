#!/bin/bash
#SBATCH -t 00:20:00
#SBATCH -J 2.3_OOB_VIP
#SBATCH --mem=0 
#SBATCH --partition=thinnodes,cola-corta
#SBATCH --ntasks=1
#SBATCH --output=/mnt/netapp2/Home_FT2/home/otras/ini/alg/2.3_OOB.out
#SBATCH --error=/mnt/netapp2/Home_FT2/home/otras/ini/alg/2.3_OOB.out

## 2.3) GET OOB ERRORS: and join to previous OOB_table.
subset=$(ls $SUBDIR | egrep "RanFor_thr[0-9]{2}_out\.txt")

for f in $subset; do
  awk 'BEGIN { FS = "\t"; OFS = "\t"; ORS = "\n" } 
  function basename(file) {sub(".*/", "", file); return file} /Mean of squared residuals/ {print $2, basename(FILENAME)}' $SUBDIR/$f | sed 's/_out//'
done > $WDIR/Filtered/OOBerr.tmp
sed -i '1i OOB_error\tRF_file' $WDIR/Filtered/OOBerr.tmp

paste -d '\t' $WDIR/Filtered/OOB_table.tmp $WDIR/Filtered/OOBerr.tmp > $WDIR/Filtered/OOB_table.tsv
rm $WDIR/Filtered/OOB_table.tmp $WDIR/Filtered/OOBerr.tmp

# Extract subsets not failed in RanFog (if used), sort by minimum OOB_error (6th column) and take the first row:
# a) best: File name (7th column); b) thr: Threshold level (2nd column)
best=$(sort -k6 -n $WDIR/Filtered/OOB_table.tsv | sed -n 2p | awk -F '\t' '{print $7}')
thr=$(sort -k6 -n $WDIR/Filtered/OOB_table.tsv | sed -n 2p | awk -F '\t' '{print $2}')

Rscript $BINPATH/Filter3_VIP.R $paramfile $phenotmod $best
