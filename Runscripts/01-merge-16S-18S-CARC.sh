#!/bin/bash

#These are for running this script on USC's CARC. Not needed/would need to be altered for users HPC.
#SBATCH --account=fuhrman_1138
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

#get variables
source 515FY-926R.cfg
timestamp=`date +"%y%m%d-%H%M"`

#these `ls` expressions will generate more than one file name if you export your data more than once!
cp 02-PROKs/10-exports/*$studyName.16S.all-16S-seqs.with-tax.tsv .
cp 02-EUKs/15-exports/*$studyName.18S.all-18S-seqs.with-PR2-tax.tsv . 
cp 02-PROKs/03-DADA2d/*stats.tsv .
cp 02-EUKs/08-DADA2d/*stats.tsv .

#This script requires the following R package:
#Tidyverse (https://www.tidyverse.org/)
#Please open 

#run script
module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r

Rscript --verbose ./scripts/01-1-correct_16S_18S_ASV.R

rm *.16S.all-16S-seqs.with-tax.tsv
rm *-18S-seqs.with-PR2-tax.tsv
rm *stats.tsv

mkdir 03-Merged

mv corrected_18S_16S_counts.tsv 03-Merged/

cd 03-Merged/

original_file="corrected_18S_16S_counts.tsv"
new_file="${studyName}_${timestamp}_corrected_18S_16S_counts.tsv"

mv $original_file $new_file



