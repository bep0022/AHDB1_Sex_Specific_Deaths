### AHDB1_Sex_Specific_Deaths

### This repository contains data and R Markdowns for assessing potential biomarkers as predictors of acute heat stress induced mortality in Zebra Finches. A secondary analysis of the significance of the deaths by age category and sex was performed using Firth's logistic regression that accounts for the small sample size and complete separation of Sex groups in the observed data. After identified as significant, these results informed our sample selection for analysis of the possible predictors of heat stress induced deaths. For our possible predictors our samples included baseline blood samples (taken 2 weeks prior to heat stress) of only younger females. Upon blood collection glucose, ketone, and hematocrit levels were measured. We found in this analysis that none were statistically significant predictors of acute heat stress induced mortality (p-value = > 0.05). The analysis also includes quantitative real-time PCR (qPCR) measurements of telomere length and mitochondrial DNA copy number. We found that neither are statistically significant predictors of acute heat stress induced mortality (p-value = > 0.05).

### Repository Outline and Summary
## Raw Data
- [AHDB_Exp1_BloodData.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_Exp1_BloodData.csv) : This is the raw data file for the blood measures of glucose, ketone, and hematocrit that were recorded during the experiment.
- [AHDB_CORT_FileForMasterDataset.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_CORT_FileForMasterDataset.csv) : This is the raw CORT data recorded during quantification.
- [AHDB1_Manuscript1_Deaths_Multiplex_2025-08-01 13-38-10_795BR20744 -  Quantification Cq Results.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_Deaths_Multiplex_2025-08-01%2013-38-10_795BR20744%20-%20%20Quantification%20Cq%20Results.csv) : This is the file from the qPCR software that is produced for the Cq quantification following the multiplex reaction of an autosomal gene (EEF2) and the mtDNA gene (ND5). 
- [AHDB1_Manuscript1_Deaths_Telomere_Heidinger_2025-08-01 11-52-51_795BR20744 -  Quantification Cq Results.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_Deaths_Telomere_Heidinger_2025-08-01%2011-52-51_795BR20744%20-%20%20Quantification%20Cq%20Results.csv) : This is the file from the qPCR software that is produced for the Cq quantification following the telomere reaction.
- [AHDB1_Manuscript1_ESEBqPCR_FinalData.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_ESEBqPCR_FinalData.csv) : This is the processed and combined data for the potential biomarkers measured using qPCR.


## R Code
- [AHDB_MS1_Data_PreProcessing.R](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_MS1_Data_PreProcessing.R) : This R code merges the raw data into a meta data file for usage in the analysis.
- [AHDB_Exp1_MS1_Analysis.Rmd](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_Exp1_MS1_Analysis.Rmd) : This R Markdown is the analysis. First, the analysis of the deaths observed, predicted proportion of deaths, and whether there is a significant effect of sex and age category. Then, the analysis of the potential biomarkers of glucose, ketone, hematocrit, CORT, telomere length, and mtDNA copy number. The potential biomarker data is from baseline blood samples taken 14 days prior to an acute heat stress event (43°C for 5 hours). This dataset contains samples from younger adult females that survived this treatment, as well as individuals that died to test whether any biomarkers were predictive of this future fate.


## Input Data
For the AHDB_Exp1_MS1_Analysis.Rmd : 
- [AHDB_Exp1_MS1.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_Exp1_MS1.csv) : This is the meta data file containing all potential biomarker metrics and specifying the samples of interest by their ID_Band. Indicated by the Died_InTrt column we have whether they did "yes" or did not died "no" due to heat stress. Other important metrics such as AgeCategory, Sex, and Treatment are included.  


## Software
All data processing and analyses were conducted on RStudio.
