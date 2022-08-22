#Demostrate that MLE performs well with small number of strata
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

# Normal sample per strata
df_dist_variables = bind_rows(
  create_simulation_scenarios(N_ref = 1000,
                              K = 5, 
                              m = 10,
                              m0 = 10,
                              FH = 5,
                              OR = 1.5,
                              sim = 1e2,
                              gamma = 1,
                              mu_delta = 0),
  # Normal sample per strata - high correlation
  create_simulation_scenarios(N_ref = 1000,
                              K = 5, 
                              m = 10,
                              m0 = 10,
                              FH = 5,
                              OR = 1.5,
                              sim = 1e2,
                              gamma = 10,
                              mu_delta = 0),
  #Sparse sample per strata
  create_simulation_scenarios(N_ref = 100,
                              K = 10, 
                              m = 50,
                              m0 = 10,
                              FH = 5,
                              OR = 1.5,
                              sim = 1e2,
                              gamma = 1,
                              mu_delta = 0), 
    create_simulation_scenarios(N_ref = 10,
                                K = 10, 
                                m = 50,
                                m0 = 10,
                                FH = 5,
                                OR = 1.5,
                                sim = 1e2,
                                gamma = 1,
                                mu_delta = 0)
  ) %>%
  mutate(N_foc = N_ref)



walk(
  # 4,
  1:4,
  # 1,
  # 2:3,
  # 4,
  # 1:nrow(df_dist_variables),
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
      sim_name = glue('mle_example_scenario_{x}.RData')
    )
    
    # debugonce(simulation_results_to_glm_input)
    
    mle_estimation(sim_results_path= glue('data/simulations/mle_example_scenario_{x}.RData'),
                   parallel=T,
                   split=T,
                   m=experiment_params$m,
                   method = 'fast'
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
                                sim_name = glue('data/mle/mle_results_mle_example_scenario_{x}.RData'),
                                output_file_name = glue('data/mle-summary/gamma/mle_results_mle_example_scenario_{x}.csv'))
    
    cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
  }, otherwise = 'Failed'))



df_dist_variables= df_dist_variables %>% 
  mutate(n_rows = row_number()) %>% 
  mutate(power_results = glue('data/mle-summary/power/mle_summary_power_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.RData'),
         gamma_results = glue('data/mle-summary/gamma/mle_results_mle_example_scenario_{n_rows}.csv'))


#Calculate P value for each simulation
power_statistics = map(df_dist_variables$power_results,
                       function(x){
                         load(x)
                         sum((1- pchisq(chisq_stats, 50)) < 0.05)/length(chisq_stats)
                       }) %>% 
  reduce(c)

df_power_stats = data.frame(
  simulation_file = df_dist_variables$simulation_file,
  power_statistics= power_statistics
)

df_sim_gamma_results = map(df_dist_variables$gamma_results,
                           function(x){
                             read_csv(x) %>% 
                               mutate(gamma_results = x)
                           }) %>% 
  bind_rows()

# df_sim_gamma_results %>% 
#   group_by(scenario, covers_true_value) %>% 
#   summarise(n = n(), mean(gamma_estimates)) %>% 
#   mutate(n_perc = n/sum(n)) %>% 
#   ungroup

df_coverage_stats = df_sim_gamma_results %>% 
  count(gamma_results, covers_true_value) %>% 
  group_by(gamma_results) %>% 
  mutate(coverage_perc = n/sum(n)) %>% 
  filter(covers_true_value) %>% 
  ungroup()

df_coverage_true_value = df_sim_gamma_results %>% 
  group_by(gamma_results) %>% 
  summarise(mean_true_value = mean(true_value))

df_gamma_estimates_conditional = df_sim_gamma_results %>% 
  group_by(gamma_results, true_value) %>% 
  summarise(gamma_hat = mean(gamma_estimates)) %>% 
  pivot_wider(names_from = true_value, values_from = gamma_hat)

df_gamma_estimates = df_sim_gamma_results %>% 
  group_by(gamma_results) %>% 
  summarise(gamma_hat = mean(gamma_estimates)) 

df_results = df_dist_variables %>% 
  left_join(df_power_stats) %>% 
  left_join(df_coverage_stats) %>% 
  left_join(df_gamma_estimates) %>% 
  left_join(df_coverage_true_value)

df_results %>% 
  select(m, K, N_ref, N_foc, gamma, FH, power_statistics, coverage_perc, gamma_hat, mean_true_value)

df_results %>% 
  select(K, N_ref, N_foc, gamma, FH, power_statistics, coverage_perc, gamma_hat, mean_true_value) %>% 
  arrange(m, K, gamma, desc(N_ref), FH) %>% 
  mutate(across(where(is.numeric),~ round(., 4))) %>% 
  write_csv('output/tables/example-simulation-results-formatted.csv')
