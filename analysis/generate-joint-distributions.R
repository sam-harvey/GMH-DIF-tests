# Apply MLE estimator / Global invariance test to demonstrate value of GMH test.

source('libraries.R')
source("analysis/Odds.R") # that's where the UTI data comes from
source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
source("analysis/Lfct.r")
source("analysis/DIF.fcts.r")
source("analysis/pisa_analysis.R")
source("analysis/simulation_study.R")

K = c(20, 100)
m = c(2, 3, 4, 5, 7, 10
      # ,50
)
gamma = c(1, 10)

df_dist_variables = crossing(K, m, gamma)

test = pmap(list(df_dist_variables$K, df_dist_variables$m, df_dist_variables$gamma),
            function(x, y, z){
              generate_joint_distribution(K=x, m=y, Gamma=z)
            })