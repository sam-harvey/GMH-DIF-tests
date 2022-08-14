rm(list = ls())

source('libraries.R')
source("analysis/Odds.R") # that's where the UTI data comes from
source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
source("analysis/Lfct.r")
source("analysis/DIF.fcts.r")
source("analysis/pisa_analysis.R")
source('analysis/simulation_study.R')
source('analysis/mle-estimate-fns.R')

#Generate simulation scenarios
source('analysis/create-simulation-scenarios.R')

df_dist_variables = create_simulation_scenarios(K = c(20),
                            m = c(50),
                            m0 = c(10),
                            gamma = c(1, 10),
                            FH = c(0, 1, 5, 10, 20),
                            N_ref = c(500, 60),
                            OR = 1.5,
                            sim = 1e3,
                            mu_delta = 0)

walk(
  1:nrow(df_dist_variables),
     possibly(function(x){
       experiment_params = df_dist_variables[x,]
       
       file.name.gen <- paste("Sim_GenDif_m",experiment_params$m,"_K",experiment_params$K,"_Gamma",experiment_params$gamma,"_OR",experiment_params$OR,".RData",sep="")
       file.name.gen = glue("data/joint-distributions/{file.name.gen}")
       
       if(!file.exists(file.name.gen)){
         generate_joint_distribution(K=experiment_params$K, 
                                     m=experiment_params$m0,
                                     Gamma=experiment_params$gamma,
                                     OR=experiment_params$OR)
       }
       
       run_simulations(N_ref = experiment_params$N_ref,
                       N_foc = experiment_params$N_foc,
                       K=experiment_params$K, 
                       m=experiment_params$m,
                       m0 = experiment_params$m0,
                       FH=experiment_params$FH,
                       OR=experiment_params$OR,
                       sim=experiment_params$sim,
                       gamma = experiment_params$gamma,
                       mu_delta = experiment_params$mu_delta,
                       simulation_file = experiment_params$simulation_file,
                       split = T,
                       parallel = T)
       
       cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
     }, otherwise = 'Failed'))
