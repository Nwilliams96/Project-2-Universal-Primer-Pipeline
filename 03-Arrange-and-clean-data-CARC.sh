#!/bin/bash

#SBATCH --account=fuhrman_1138
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00

#Automated script for adding Proportal classifications to an ASV table
source 515FY-926R.cfg
timestamp=`date +"%y%m%d-%H%M"`

#Move all the required files for "03-1-arrange-and-clean-sample-data.R" into this working directory
cp 03-Merged/*corrected_18S_16S_counts_ProPortal.tsv . #Merged data from last script
cp ../Longhurst.output.csv . #Longhurst information - this is gained from Bror Johnsons script.
cp ../grump_sat_matched_20240714.csv . #Bror Euphotic zone calculation.

#Move all the required files for "03-2-long-data-preparation.R"
cp 02-PROKs/03-DADA2d/*16S.dna-sequences.fasta .
cp 02-EUKs/08-DADA2d/*18S.dna-sequences.fasta .

module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r

Rscript --verbose ./scripts/03-1-arrange-and-clean-sample-data.R

original_file_sample_data="sample-metadata.csv"
new_file_sample_data="${studyName}-sample-metadata.csv"
mv $original_file_sample_data $new_file_sample_data

Rscript --verbose ./scripts/03-2-long-data-preparation.R

original_file_asv="asv_long.csv"
new_file_asv="${studyName}_${timestamp}_asv_long.csv"
mv $original_file_asv $new_file_asv

original_file_asv="asv_long_CMAP.csv"
new_file_asv="${studyName}_${timestamp}_asv_long_CMAP.csv"
mv $original_file_asv $new_file_asv

mkdir 04-Formatted
mv *asv_long* 04-Formatted

mkdir 05-Sequences
mv *sequences* 05-Sequences

rm *corrected_18S_16S_counts_ProPortal.tsv
