##EJOR 2010

This repo stores some routines to assess the validity of the R code in reproducing the MATLAB code from [Optimal bandwidth selection for conditional efficiency measures: A data-driven approach](http://www.sciencedirect.com/science/article/pii/S0377221709002148) by Luiza Badin, Cinzia Daraio , Léopold Simar.
The simulation assessment is performed using the DGPs illustrated in the abovementioned paper.
###Description of files in the Repo


* README.md : this file
* seeds.csv : seeds used in validation to generate data. If used, the results presented in "validation_assessment" can be reproduced 
* validation.m : script that reads the datasets produced by validation_ban.R and computes LSCV for every observation in the sample using Ker_LSCV_OUT.m
* validation.R : script that produces the dataset using the DGPs considered in the abovementioned article, saves datasets in .mat format and computes LSCV for every observation in the sample using Ker_LSCV_OUT.R
* validation_assessment.Rmd : R markdown document that illustrates the validation procedures and results 