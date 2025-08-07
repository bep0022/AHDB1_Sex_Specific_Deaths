### AHDB1_Sex_Specific_Deaths

### This repository contains data and R Markdowns for assessing potential biomarkers as predictors of acute heat stress induced mortality in Zebra Finches. A secondary analysis of the significance of the deaths by age category and sex was performed using Firth's logistic regression that accounts for the small sample size and complete separation of Sex groups in the observed data. After identified as significant, these results informed our sample selection for analysis of the possible predictors of heat stress induced deaths. For our possible predictors our samples included baseline blood samples (taken 2 weeks prior to heat stress) of only younger females. Upon blood collection glucose, ketone, and hematocrit levels were measured. We found in this analysis that none were statistically significant predictors of acute heat stress induced mortality (p-value = > 0.05). The analysis also includes quantitative real-time PCR (qPCR) measurements of telomere length and mitochondrial DNA copy number. We found that neither are statistically significant predictors of acute heat stress induced mortality (p-value = > 0.05).

## Overview
```
add file tree?
```

### Repository Outline and Summary
## R Markdowns
In this respository, we have opted for the use of R Markdowns instead of R Scripts.
- [AHDB_Manuscript1_Processing&Analysis.Rmd](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_MS1_Processing/AHDB_Manuscript1_Processing%20%26%20Analysis.Rmd) : This R Markdown is the analysis of the deaths observed, predicted proportion of deaths, and whether there is a significant effect of sex and age category. Also this includes the analysis of the blood measures of glucose, ketone, and hematocrit.
- [Final_AHDB1_mtDNA_Telomeres.Rmd](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/Final_AHDB1_mtDNA_Telomeres.Rmd) : This R Markdown is for the analysis of telomere length and mtDNA copy number qPCR data. This data is from baseline blood samples taken 14 days prior to an acute heat stress event (43°C for 5 hours). This dataset contains samples from younger adult females that survived this treatment, as well as individuals that died. We are testing whether telomere length and/or mtDNA copy number were predictive of this future fate.

## Input Data
For the For the AHDB_Manuscript1_Processing&Analysis.Rmd : 
- [AHDB_Exp1_BloodData.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_MS1_Processing/AHDB_Exp1_BloodData.csv) : This is the file specifying the samples of interest by their ID_Band. Indicated by the Died_InTrt column we have whether they did "yes" or did not died "no" due to heat stress. Other important metrics such as AgeCategory, Sex, and Treatment are included.  


For the AHDB1_mtDNA_Telo.Rmd : 
- [AHDB1_Manuscript1_Deaths_Multiplex_2025-08-01 13-38-10_795BR20744 -  Quantification Cq Results.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_Deaths_Multiplex_2025-08-01%2013-38-10_795BR20744%20-%20%20Quantification%20Cq%20Results.csv) : This is the file from the qPCR software that is produced for the Cq quantification following the multiplex reaction of an autosomal gene (EEF2) and the mtDNA gene (ND5). 
- [AHDB1_Manuscript1_Deaths_Telomere_Heidinger_2025-08-01 11-52-51_795BR20744 -  Quantification Cq Results.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_Deaths_Telomere_Heidinger_2025-08-01%2011-52-51_795BR20744%20-%20%20Quantification%20Cq%20Results.csv) : This is the file from the qPCR software that is produced for the Cq quantification following the telomere reaction.
- [Trait_MetaData.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/Trait_MetaData.csv) : This is a file that contains phenotypic data of the Samples. It is merged with the qPCR data post filtering of outliers, in preparation for input for the telomere and mtDNA statistical analysis. 

## Output Data
For the AHDB_Manuscript1_Processing&Analysis.Rmd : 
- [Combined_Blood_Plot.png](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB_MS1_Processing/Combined_Blood_Plot.png) : This is the final plots from the models ran for Glucose, Ketone, and Hematocrit values and their lack of significance (p-value = > 0.05) between Died_InTrt.

For the AHDB1_mtDNA_Telo.Rmd : 
- [AHDB1_Manuscript1_ESEBqPCR_FinalData.csv](https://github.com/bep0022/AHDB1_Sex_Specific_Deaths/blob/main/AHDB1_Manuscript1_ESEBqPCR_FinalData.csv) : This is the output from the AHDB1_mtDNA_Telo.Rmd file. This file contains the merged contents of the filtered telomere and mtDNA qPCR data and the Trait_MetaData.csv file. 


## Software
All data processing and analyses were conducted on RStudio.
