
# Title: Acute Heat Damage Biomarkers Experiment 1 (Manuscript 1 ~ Baseline Biomarkers & Acute Heat Mortality)

# Purpose: This script was used to merge data from qPCR (mitochondria and telomeres), hormone data (CORT), with the experiment master spreadsheet 

# Clear memory
rm(list=ls(all = TRUE))

# Set working directory (assumes data file is available in same location as code)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

library(dplyr)
library(readr)
library(stringr)

# Read CSV Files

MasterSpreadsheet <- read.csv("AHDB_Exp1_BloodData.csv", header=T, sep = ",", as.is=T)
str(MasterSpreadsheet)

MitoTelo <- read.csv("AHDB1_Manuscript1_ESEBqPCR_FinalData.csv", header=T, sep = ",", as.is=T)
str(MitoTelo)

CORT <- read.csv("AHDB_CORT_FileForMasterDataset.csv", header=T, sep = ",", as.is=T)
str(CORT)

# Clean MitoTelo Dataset 
# 1. Set Sample Column Name to ID_Band,
# 2. Remove B in Sample Column (B = Baseline) 
# There exists a seperate timepoint column hence B in Sample column (now ID_Band) can be removed)

MitoTelo_clean <- MitoTelo %>%
  mutate(
    ID_Band = as.numeric(str_extract(Sample, "\\d+"))  # keeps numeric part corresponding to bird_ID 
  ) %>%
  select(-Sample)  # Removes old Sample column
head(MitoTelo_clean)

# Subset Mito Columns Needed
mito_subset <- MitoTelo_clean[, c("ID_Band", "TimePoint", "SQ.Mean_SCNAG", "SQ.Mean_mtDNA", "mtDNA.Mean",
                        "Cq.Mean_SCNAG", "SQ.Mean_Telomeres", "Cq.Mean_Telomeres", 
                        "Telomeres.per.cell")]

# Merge by ID_Band and TimePoint
MitoBlood_Merge <- merge(MasterSpreadsheet, mito_subset, 
                      by = c("ID_Band", "TimePoint"), all.x = TRUE)
head(MitoBlood_Merge)

# Clean CORT Dataset
# 1. Align Column Names to match MitoBlood Merge
# 2. Change ID Band to integer 

# Fix Column Names
colnames(CORT)[colnames(CORT) == "BirdID"] <- "ID_Band"
colnames(CORT)[colnames(CORT) == "Timepoint"] <- "TimePoint"

# Ensure matching types
CORT$ID_Band <- as.integer(CORT$ID_Band)

# Subset CORT Columns Needed

CORT_subset <- CORT[, c("PlateNumber", "TubeID", "ID_Band", "TimePoint", "Extractor", "CORT.ng.mL.", "Notes")]


# Merge CORT into MasterSpreadsheet by ID_Band and TimePoint
AHDB_Exp1_Final <- merge(MitoBlood_Merge, CORT_subset, 
                      by = c("ID_Band", "TimePoint"), all.x = TRUE)

#Export to csv
write.csv(AHDB_Exp1_Final, 
          file = "AHDB_Exp1_MasterSpreadSheet_Final.csv", 
          row.names = FALSE)

#Subset Data for MS1
AHDB_Exp1_MS1 <- AHDB_Exp1_Final %>%
  filter(TimePoint == "Baseline",
         Treatment     != "A")

#Export to csv
write.csv(AHDB_Exp1_MS1, 
          file = "AHDB_Exp1_MS1.csv", 
          row.names = FALSE)
