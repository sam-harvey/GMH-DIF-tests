dir.create('data')
dir.create('output')
dir.create('output/tables')
dir.create('output/plots')

source('renv/activate.R')
renv::restore()