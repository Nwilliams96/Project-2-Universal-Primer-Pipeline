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

module purge
module load gcc/12.3.0
module load git/2.42.0

#clone into external repo
git clone https://github.com/jcmcnch/ProPortal-ASV-Annotation.git
mkdir -p ProPortal-intermediates

#get ASV table names
counts=`ls 03-Merged/*corrected_18S_16S_counts.tsv`

#get DNA sequence file name for denoised 16S
dnaseqs=`ls 02-PROKs/03-DADA2d/*dna-sequences.fasta`

#Need to make seqkit environment
eval "$(conda shell.bash hook)"
conda activate seqkit

#get ASV ids from denoised data
grep "Synechococcales" $counts | cut -f1 > ProPortal-intermediates/Synechococcales.ids
seqkit grep -f ProPortal-intermediates/Synechococcales.ids $dnaseqs > ProPortal-intermediates/Synechococcales.fasta

conda deactivate

module load python

#enter repo
cd ProPortal-ASV-Annotation
cp ../ProPortal-intermediates/Synechococcales.fasta ASVs-2-classify
sbatch /scripts/07-blast-all-datasets.sh
./scripts/08-classify-ASVs-with-ProPortal.py blast-results/Synechococcales.blastout.tsv > ../ProPortal-intermediates/Synechococcales.Proportal-classified.tsv

cd ..

cp ProPortal-intermediates/Synechococcales.Proportal-classified.tsv .
cp 03-Merged/*corrected_18S_16S_counts.tsv .

module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r

Rscript --verbose ./scripts/02-9-attach-proportal-assignment.R

original_file="corrected_18S_16S_counts_ProPortal.tsv"
new_file="${studyName}_${timestamp}_corrected_18S_16S_counts_ProPortal.tsv"

mv $original_file $new_file

mv *corrected_18S_16S_counts_ProPortal.tsv 03-Merged

rm Synechococcales.Proportal-classified.tsv

