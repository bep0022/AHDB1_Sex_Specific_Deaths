#### AHDB1_Sex_Specific_Deaths

### This repository contains data and R Markdowns for assessing potential biomarkers as predictors of acute heat stress induced death in Zebra Finches. A preliminary analysis of the significance of the deaths by age category and sex was performed. After identified as significant, these results informed our sample selection for analysis of the possible predictors of heat stress induced deaths. The data includes quantitative real-time PCR (qPCR) measurements of telomere length and mitochondrial DNA copy number. We found that neither are statistically significant predictors (p-value = > 0.05).

## Overview
```
├── AHDB1_Death_TS.Rmd
├── AHDB1_Death_TS.md
├── AHDB1_Death_TS.pdf
├── AHDB1_Death_TS_files
│   └── figure-gfm
│       ├── unnamed-chunk-11-1.png
│       ├── unnamed-chunk-15-1.png
│       ├── unnamed-chunk-5-1.png
│       └── unnamed-chunk-8-1.png
├── AHDB1_Sex_Specific_Deaths.Rproj
├── AHDB1_mtDNA_Telo.Rmd
├── AHDB1_mtDNA_Telo.md
├── AHDB1_mtDNA_Telo.pdf
├── AHDB1_mtDNA_Telo_files
│   └── figure-gfm
│       ├── unnamed-chunk-23-1.png
│       ├── unnamed-chunk-24-1.png
│       ├── unnamed-chunk-25-1.png
│       ├── unnamed-chunk-26-1.png
│       ├── unnamed-chunk-27-1.png
│       ├── unnamed-chunk-28-1.png
│       ├── unnamed-chunk-30-1.png
│       └── unnamed-chunk-31-1.png
├── AHDB_MasterDataSheet.csv
├── Death_TS.nb.html
├── Manuscript1_DEATH_SCNAG_2025-04-20 19-13-14_795BR20744 -  Quantification Cq Results.csv
├── Manuscript1_DEATH_Telo_2025-04-20 20-35-37_795BR20744 -  Quantification Cq Results.csv
├── Plate1_FinalData.csv
├── README.md
├── Trait_MetaData.csv
└── Younger.csv
```

### Repository Outline and Summary
## R Markdowns
In this respository, we have opted for the use of R Markdowns instead of R Scripts.
- AHDB1_Death_TS.Rmd : This R Markdown is the preliminary analysis of the deaths observed, and whether there is a significant effect of sex and age category. 
- AHDB1_mtDNA_Telo.Rmd : This R Markdown is for the analysis of telomere length and mtDNA copy number qPCR data. This data is from baseline blood samples taken 14 days prior to an acute heat stress event (43°C for 5 hours). This dataset contains samples from individuals that survived this treatment, as well as individuals that died. We are testing whether telomere length and/or mtDNA copy number were predicitve of this future fate.

## Input Data
For the AHDB1_Death_TS.Rmd :
- AHDB_MasterDataSheet.csv : This is the specifying the samples of interest by their ID_Band. Indicated by the Died_InTrt column we have whether they did "1" or did not died "0" due to heat stress. Other important metrics such as AgeCategory, Sex, and Treatment are included.  
- Younger.csv : This is a subset of the master data sheet that only includes younger individuals.

For the AHDB1_mtDNA_Telo.Rmd : 
- Manuscript1_DEATH_SCNAG_2025-04-20 19-13-14_795BR20744 -  Quantification Cq Results.csv : This is the file from the qPCR software that is produced for the Cq quantification following the multiplex reaction of an autosomal gene (EEF2) and the mtDNA gene (ND5). 
- Manuscript1_DEATH_Telo_2025-04-20 20-35-37_795BR20744 -  Quantification Cq Results.csv : This is the file from the qPCR software that is produced for the Cq quantification following the telomere reaction.
- Trait_MetaData.csv : This is a file that contains phenotypic data of the Samples. It is merged with the qPCR data post filtering of outliers, in preparation for input for the telomere and mtDNA statistical analysis. 

## Output Data
This repository contains the output data files. 
- Plate1_FinalData.csv : This is the output from the AHDB1_mtDNA_Telo.Rmd file. This file contains the merged contents of the filtered telomere and mtDNA qPCR data and the Trait_MetaData.csv file. 
