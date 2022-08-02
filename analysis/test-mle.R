#Test that MLE performs well with small number of strata
#and lack of highly correlated items

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

experiment_params = create_simulation_scenarios(N_ref = 1000,
                            K=10, 
                            m=50,
                            m0=10,
                            FH=1,
                            OR=1.5,
                            sim=1e1,
                            gamma = 1,
                            mu_delta = 0) %>% 
  mutate(N_foc = N_ref) %>% 
  mutate(version = 2) %>% 
  mutate(sim_name = glue('mle_example_scenario_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sim}_version{version}.RData'),
         sim_results_path = glue('data/simulations/mle_example_scenario_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sim}_version{version}.RData'),
         mle_name = glue('data/mle/mle_results_mle_example_scenario_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sim}_version{version}.RData'),
         output_file_name = glue('data/mle-summary/gamma/mle_results_mle_example_scenario_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sim}_version{version}.csv'))

# experiment_params = data.frame(
#   Gamma=1,
#   K=5, 
#   m=50,
#   m0 = v11,
#   mu.delta=0,
#   FH=5,
#   K=5,
#   # N_ref=500,
#   # N_foc=500,
#   OR=1.5,
#   BT=1,
#   sim= 1e1
# )

file.name.gen <- paste("Sim_GenDif_m",experiment_params$m,"_K",experiment_params$K,"_Gamma",experiment_params$gamma,"_OR",experiment_params$OR,".RData",sep="")
file.name.gen = glue("data/joint-distributions/{file.name.gen}")

if(!file.exists(file.name.gen)){
  
  generate_joint_distribution(K=experiment_params$K, 
                              m=experiment_params$m0,
                              Gamma=experiment_params$gamma,
                              OR=experiment_params$OR)
}

dif_sims = dif_simulation(
  Gamma=experiment_params$gamma,
  m0=experiment_params$m0,
  m=experiment_params$m,
  mu.delta=experiment_params$mu_delta,
  FH=experiment_params$FH,
  K=experiment_params$K,
  N_ref=experiment_params$N_ref,
  N_foc=experiment_params$N_foc,
  OR=experiment_params$OR,
  sim= experiment_params$sim,
  zeros=FALSE,
  within.group=F,
  seed_val=1,
  use_bt=F,
  sim_name = experiment_params$sim_name
)

# debugonce(simulation_results_to_glm_input)

mle_estimation(sim_results_path=experiment_params$sim_results_path,
               parallel=T,
               split=T,
               m=experiment_params$m,
               method = 'fast'
               # method = 'base'
               # , full_model_formula = response ~ question:k + question:C(group, contr = contr.SAS(2)) - 1,
               # reduced_model_formula = response ~  question:k - 1,
               # , full_model_formula = response ~ question + question:C(group, contr = contr.SAS(2)) - 1,
               # reduced_model_formula = response ~  question - 1,
               )

calculate_mle_summary_stats(gamma=experiment_params$gamma,
                            m=experiment_params$m,
                            FH=experiment_params$FH,
                            K=experiment_params$K,
                            N_ref=experiment_params$N_ref,
                            N_foc=experiment_params$N_foc,
                            mu_delta=experiment_params$mu_delta,
                            OR=experiment_params$OR,
                            sims = experiment_params$sim,
                            sim_name = experiment_params$mle_name,
                            output_file_name = experiment_params$output_file_name)

load(experiment_params$mle_name)
experiment_gamma_results = read_csv(experiment_params$output_file_name)

experiment_gamma_results %>% 
  ggplot + 
  geom_density(aes(x=gamma_estimates)) 
  # coord_cartesian(xlim = c(-5,5))

data.frame(
  coefficient = model_results[[1]]$model_matrix %>% colnames(),
  value = model_results[[1]]$full_model_fit$coefficients
)

data.frame(
  coefficient = model_results[[1]]$model_matrix %>% colnames(),
  value = model_results[[1]]$full_model_fit$coefficients
)


map(model_results,
    function(x){
      data.frame(
        coefficient = x$model_matrix %>% colnames(),
        value =x$full_model_fit$coefficients
      ) %>% 
        as.tbl()
    }) %>% 
  bind_rows() %>% 
  group_by(coefficient) %>% 
  summarise(mean(value)) %>% 
  filter(str_detect(coefficient, 'gr'))

model_results[[1]]$chisq_test

map(model_results, ~.$chisq_test) %>% 
  reduce(c) %>% 
  pchisq(., 10) 
