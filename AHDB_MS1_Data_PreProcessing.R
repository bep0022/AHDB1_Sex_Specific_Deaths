
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

colnames(mito_subset)[3:9] <- c(
  "Blood_SQ.Mean_SCNAG",
  "Blood_SQ.Mean_mtDNA",
  "Blood_mtDNA.Mean",
  "Blood_Cq.Mean_SCNAG",
  "Blood_SQ.Mean_Telomeres",
  "Blood_Cq.Mean_Telomeres",
  "Blood_Telomeres.per.cell"
)

# Merge by ID_Band and TimePoint
MitoBlood_Merge <- merge(MasterSpreadsheet, mito_subset, 
                      by = c("ID_Band", "TimePoint"), all.x = TRUE)
head(MitoBlood_Merge)

# Clean CORT Dataset
# 1. Align Column Names to match MitoBlood Merge
# 2. Change ID Band to integer 

# Fix Column Names
colnames(CORT)[colnames(CORT) == "BirdID"] <- "ID_Band"

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

### Data Processing & Filtering - qPCR Measures of mtDNA and Telomeres
# Import raw qPCR output for each plate 1) the Single Copy Autosomal Gene (EEF2) and mtDNA gene multiplex and 2) the Telomere reaction. 

SCNAG <- read.csv("AHDB1_Manuscript1_Deaths_Multiplex_2025-08-01 13-38-10_795BR20744 -  Quantification Cq Results.csv")
dim(SCNAG)

Telo <- read.csv("AHDB1_Manuscript1_Deaths_Telomere_Heidinger_2025-08-01 11-52-51_795BR20744 -  Quantification Cq Results.csv")
dim(Telo)
```

# Edit 

# Concatenate data across runs for the same samples
Plate1 <- rbind(SCNAG, Telo)
dim(Plate1)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate1$PlateID <- "Plate1"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate1$Target[Plate1$Fluor == "VIC"] <- "scnag"
Plate1$Target[Plate1$Fluor == "FAM"] <- "mtdna"
Plate1$Target[Plate1$Fluor == "SYBR"] <- "telomeres"
```

# Identify and remove outliers (across the three replicates)

# Add a column called "Diff_AVG-Cq"
Plate1$Diff_AVG_Cq <- abs(Plate1$Cq - Plate1$Cq.Mean)

## Add a column called "Flag_outlier"
Plate1$Flag_outlier <- ifelse(Plate1$Diff_AVG_Cq > 0.4001, "yes", "no")
HighCq <- Plate1[Plate1$Diff_AVG_Cq > "0.3001", ]
VeryHighCq <- Plate1[Plate1$Diff_AVG_Cq > "0.4001", ]

# Report number of rows to be removed
print(paste("Number of rows with High Cq:", nrow(VeryHighCq)))
# Make a table called  of the high cq values
HighCq_Samples<- HighCq[c("Well", "Sample", "Fluor", "Diff_AVG_Cq" )]
print(HighCq_Samples)


# Filter and remove rows where Diff_AVG_Cq > 0.4001 for each SampleID
# USed to be max(AVG) but error said that col did not exist
dim(Plate1)
Plate1 <- Plate1 %>%
  group_by(Sample, Target) %>%
  filter(!(Diff_AVG_Cq > 0.4001 & Diff_AVG_Cq == max(Diff_AVG_Cq)))
dim(Plate1)

```

# Additional identification and removal of outliers

# Recalculate Cq Mean and Add a column called "Cq Mean2"
Plate1 <- Plate1 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    SQ.Mean2 = mean(Starting.Quantity..SQ.)) %>%
  ungroup()
dim(Plate1)

# Add a column called "Diff_AVG-Cq_2" that contains the difference between new Cq and new mean
Plate1$Diff_AVG_Cq_2 <- abs(Plate1$Cq - Plate1$Cq.Mean2)

# Add a column called "Flag_outlier_2"
Plate1$Flag_outlier_2 <- ifelse(Plate1$Diff_AVG_Cq_2 > 0.401, "yes", "no")
VeryHighCq_2 <- Plate1[Plate1$Diff_AVG_Cq > "0.401", ]
# Report number of rows to be removed
print(paste("Number of rows with High Cq:", nrow(VeryHighCq_2)))
# Make a table called  of the high cq values
HighCq_2_Samples<- HighCq[c("Well", "Sample", "Fluor", "Diff_AVG_Cq" )]
print(HighCq_2_Samples)

# Remove the rows with a high Cq Outliers
# USed to be max(AVG) but error said that col did not exist

Plate1 <- Plate1 %>%
  group_by(Sample, Target) %>%
  filter(!(Diff_AVG_Cq_2 > 0.4001 & Diff_AVG_Cq_2 == max(Diff_AVG_Cq_2)))

dim(Plate1)

###############  If any "Samples" do not have more than two rows (replicates left), remove
Plate1 <- Plate1 %>%
  group_by(Sample, Target) %>%
  filter(n() >= 2) %>%
  ungroup()
dim(Plate1)


# Remove negative controls and standards

# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate1 <- Plate1 %>%
  filter(!str_detect(Content, "Std-*")) %>%
  filter(!str_detect(Content, "NTC")) %>%
  filter(!str_detect(Sample, "Neg")) %>%
  filter(!str_detect(Sample, "STD*"))
dim(Plate1)


# Subset dataset based on the value in "Target" column

unique_targets <- unique(Plate1$Target)
# Create a list to store the subset dataframes
subset_dfs <- list()
# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate1, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}


# Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value

Plate1_SCNAG<-print(subset_dfs[["scnag"]])
Plate1_SCNAG<-Plate1_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate1_SCNAG <- Plate1_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate1_mtDNA<-print(subset_dfs[["mtdna"]])
Plate1_mtDNA<-Plate1_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate1_mtDNA <- Plate1_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate1_Telomeres<-print(subset_dfs[["telomeres"]])
Plate1_Telomeres<-Plate1_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate1_Telomeres <- Plate1_Telomeres %>% 
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)


#  Make a final MPX dataset for by merging the Target datasets horizontally, in rows.

Plate1_FinalMPX <- merge(Plate1_SCNAG, Plate1_mtDNA, by = c("PlateID", "Well", "Sample"))

# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate1_FinalMPX$mtDNA <- (Plate1_FinalMPX$SQ.Mean_mtDNA / Plate1_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate1_FinalMPX <- Plate1_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()


# Merge final MPX with Telomeres

## Reduce datasets to single row per individual containing only the columns we want.
Plate1_FinalMPX <- distinct(Plate1_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate1_FinalTelo <- distinct(Plate1_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate1_FinalData <- merge(Plate1_FinalMPX, Plate1_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate1_FinalData <- Plate1_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)


# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate1_FinalData <- Plate1_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()


# Merge Final Data with Trait MetaData for your individuals
# Load in Data
Trait <- read.csv("Trait_MetaData.csv")
dim (Trait)

# Merge both datasets
FinalData <- merge(Plate1_FinalData, Trait, by = c("Sample"))



# Write the final data file for this plate
write.csv(file = "AHDB1_Manuscript1_ESEBqPCR_FinalData.csv", FinalData, row.names = FALSE)

