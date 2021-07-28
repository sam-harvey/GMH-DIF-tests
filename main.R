source('renv/activate.R')

#Processes and subsets PISA data for use in analysis
source('analysis/create_pisa_analysis_data.R')
#Runs the results
source('analysis/pisa.R')
#Creates output tables
source('analysis/tables.R')