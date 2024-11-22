#!/usr/bin/env Rscript

#Purpose: normalize raw counts of 16S and 18S ASVs by dividing the 16S reads and 18S reads by the % passing DADA2, 
#then converting to proportions and multiplying counts of 18S sequences by the bias against them and adjusting the 16S sequences to match that bias (as user specifies).
#Required packages: Rscript, args, tidyverse
#Output: This returns a file with ASV relative abundances out of (16S + 18S).
#Note: we recommend assuming a 2-fold bias against 18S sequences, which has been found with Illumina HiSeq or MiSeq data (Yeh et al. 2018)
#This script must be run from the base directory (the folder that contains 02-PROKs/ and 02-EUKs/)
#Author: Nathan Lloyd Robert Williams
#Final version: 11.16.2022

#Note that there will be Eukaryotes in your 16S. These come from the 16S chloroplast, and are phytoplankton.

#Set directory for libraries.
myPaths <- .libPaths()
myPaths <- c(myPaths, '/home1/nathanwi/R/x86_64-pc-linux-gnu-library/4.3')
.libPaths(myPaths)

#Load dependencies
library(tidyverse)

#Load in data
# List and read files, combine into a single data.table
asv_table <- list.files(pattern = '*corrected_18S_16S_counts.tsv') %>%
  # Use map to read and clean each file
  map_dfr(~ read_delim(.x, delim = "\t", col_names = TRUE)) %>%
  # Convert to tibble for easier row handling
  as_tibble()

ProPortal_Assigned_Cyanobacteria <- read_tsv("Synechococcales.Proportal-classified.tsv", col_names = FALSE)
colnames(ProPortal_Assigned_Cyanobacteria) <- c("ASV_hash", "ProPortal_ASV_Ecotype")

asv_table <- asv_table %>% left_join(ProPortal_Assigned_Cyanobacteria) %>%
  select(1:2, ProPortal_ASV_Ecotype, everything())

write_tsv(asv_table, "corrected_18S_16S_counts_ProPortal.tsv")

q()
