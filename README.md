# Introduction

This analysis produces results included in our forthcoming paper 'Generalized Mantel--Haenszel Estimators for
Simultaneous Differential Item Functioning Tests' by Liu et al. 
The main results are the application of the derived simultaneous test of DIF to data from the 2012 
Programme for International Student Assessment (PISA).

# Installation

This uses R 3.6.3. Install these versions of [R](https://cran.r-project.org/bin/windows/base/old/3.6.3/) and [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html). Then run the setup file to install dependencies.

`
source('setup.R')
`

The raw data used in this analysis is available via [Kaggle](https://www.kaggle.com/samuelh/pisa-2012).
The file download-data.sh provides the kaggle API command to download this.

# Output

To reproduce the final output run after installing and downloading the dataset run `main.R`.
This will generate several files in 'data/output'.
These are

* pisa_results.RData - the main results including our test results
* table csvs - the output tables of the PISA section are provided in CSV format