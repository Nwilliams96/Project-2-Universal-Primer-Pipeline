#!/usr/bin/env Rscript

#Set directory for libraries.
myPaths <- .libPaths()
myPaths <- c(myPaths, '/home1/nathanwi/R/x86_64-pc-linux-gnu-library/4.3')
.libPaths(myPaths)

#Load dependencies
library(tidyverse)

#Load in data
# List and read files, combine into a single data.table
asv_table <- list.files(pattern = '*corrected_18S_16S_counts_ProPortal.tsv') %>%
  # Use map to read and clean each file
  map_dfr(~ read_delim(.x, delim = "\t", col_names = TRUE)) %>%
  # Convert to tibble for easier row handling
  as_tibble()

#Import the sample data
sample_data_list <- list.files(pattern = '*metadata.csv')
sample_data_df <- lapply(sample_data_list, readr::read_csv)
sample_data <- data.table::rbindlist(sample_data_df, use.names = TRUE, fill = TRUE)

#Merge Date and Time to be in the format that CMAP needs
sample_data <- sample_data %>% separate(Date, into = c("Day","Month", "Year"), sep = "/")
sample_data$Year <- as.numeric(sample_data$Year)
sample_data <- sample_data%>% mutate(Year_CMAP = Year + 2000)
sample_data$Year <- as.character(sample_data$Year)
sample_data <- sample_data %>% unite(Date_CMAP, c("Year_CMAP","Month","Day"), sep='-', remove=F) #unite year-month-day and then a "T" and then the time
sample_data <- sample_data %>% unite(Date_Time, c("Date_CMAP","Time"), sep='T', remove=F)
sample_data <- sample_data %>% mutate(Date_Time = paste0(Date_Time, "Z"))
sample_data <- sample_data %>% select(-c("Date_CMAP","Year_CMAP"))
sample_data <- sample_data %>% filter(!Cruise_ID %in% (NA))
sample_data <- sample_data %>% unite(Date, c("Day","Month","Year"), sep='/', remove=F)

#Call the different seasons
Northern.Hemisphere.Spring <- sample_data %>% filter(Latitude >= 0, Month >= 3, Month <= 5) %>% mutate(Season = c("Spring"))
Northern.Hemisphere.Summer <- sample_data %>% filter(Latitude >= 0, Month >= 6, Month <= 8) %>% mutate(Season = c("Summer"))
Northern.Hemisphere.Autumn <- sample_data %>% filter(Latitude >= 0, Month >= 9, Month <= 11) %>% mutate(Season = c("Autumn"))
Northern.Hemisphere.Winter.1 <- sample_data %>% filter(Latitude >= 0, Month >= 12) %>% mutate(Season = c("Winter"))
Northern.Hemisphere.Winter.2 <- sample_data %>% filter(Latitude >= 0, Month <= 2) %>% mutate(Season = c("Winter"))
Southern.Hemisphere.Spring <- sample_data %>% filter(Latitude < 0, Month >= 9, Month <= 11) %>% mutate(Season = c("Spring"))
Southern.Hemisphere.Summer.1 <- sample_data %>% filter(Latitude < 0, Month >= 12) %>% mutate(Season = c("Summer"))
Southern.Hemisphere.Summer.2 <- sample_data %>% filter(Latitude < 0, Month <= 2) %>% mutate(Season = c("Summer"))
Southern.Hemisphere.Autumn <- sample_data %>% filter(Latitude < 0, Month >= 3, Month <= 5) %>% mutate(Season = c("Autumn"))
Southern.Hemisphere.Winter <- sample_data %>% filter(Latitude < 0, Month >= 6, Month <= 8) %>% mutate(Season = c("Winter"))

Seasons <-bind_rows(Northern.Hemisphere.Spring, Northern.Hemisphere.Summer, Northern.Hemisphere.Autumn, Northern.Hemisphere.Winter.1, Northern.Hemisphere.Winter.2, Southern.Hemisphere.Spring, Southern.Hemisphere.Summer.1, Southern.Hemisphere.Summer.2, Southern.Hemisphere.Autumn, Southern.Hemisphere.Winter)

sample_data <- sample_data %>% left_join(Seasons) %>% distinct(SampleID, .keep_all = TRUE)

#Longhurst_province_csv - this also should be added to your folder.
longhurst <- read_csv("Longhurst.output.csv")

#Join in the sample data
sample_data <- sample_data %>% left_join(longhurst)

#Load in Bror's Euphotic zone estimate
Euphotic <- read_csv("grump_sat_matched_20240714.csv")

#Join this into the sample data
Euphotic <- Euphotic %>%
  select(-c("chlor_a","kd_490","sst","Time","Date"))
sample_data <- sample_data %>% 
  left_join(Euphotic, by = c("Latitude", "Longitude", "Depth","Cruise_ID"))

write_csv(sample_data, "sample-metadata.csv")

q()
