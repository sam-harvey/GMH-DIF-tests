dir.create('data')
dir.create('output')
dir.create('output/tables')

source('renv/activate.R')
renv::restore()