
source('libraries.R')
source("analysis/Odds.R") # that's where the UTI data comes from
source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
source("analysis/Lfct.r")
source("analysis/DIF.fcts.r")
source("analysis/pisa_analysis.R")
source('analysis/simulation_study.R')
source('analysis/mle-estimate-fns.R')

#Generate simulation data
# debugonce(dif_simulation)
# debugonce(dif_simulate)

# debugonce(dif_simulation)

K = c(20, 100)
m = c(50)
m0 = c(10)
gamma = c(1, 10)
FH = c(0, 1, 5
       ,10, 20
       )
N_ref = c(1000, 500)
OR = 1.5
sim = 10
# N_foc=c(200, 500)

df_dist_variables = crossing(K, m, m0, gamma, FH, N_ref, OR, sim) %>% 
  mutate(N_foc = case_when(N_ref == 1000 ~ 200,
                           N_ref == 500 ~ 500)) %>% 
  mutate(simulation_file = glue("data/simulations/Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sim}.RData") %>% 
           as.character())

#Generate simulated test performance
walk(1:nrow(df_dist_variables),
     function(x){
       experiment_params = df_dist_variables[x,]
       
       dif_sims = dif_simulation(
         Gamma=experiment_params$gamma,
         m0=experiment_params$m0,
         m=experiment_params$m,
         mu.delta=0,
         FH=experiment_params$FH,
         K=experiment_params$K,
         N_ref=experiment_params$N_ref,
         N_foc=experiment_params$N_foc,
         OR=1.5,
         BT=1,
         sim= 1e1,
         zeros=FALSE,
         within.group=F,
         seed_val=1,
         use_bt=F
       )
       
       cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
       
       gc()
     })

#Find MLE estimates
walk(1:nrow(df_dist_variables),
     function(x){
       experiment_params = df_dist_variables[x,]
       
       mle_estimation(sim_results_path=experiment_params$simulation_file,
                      m=experiment_params$m,
                      parallel=T)
       
       cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
       
       gc()
     })
