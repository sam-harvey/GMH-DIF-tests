# Summarise simulation results and write publication tables

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

df_dist_variables = create_simulation_scenarios(K = c(20, 5),
                                                m = c(50),
                                                m0 = c(10),
                                                gamma = c(1, 10),
                                                FH = c(0, 1, 5, 10, 20),
                                                N_ref = c(500, 60),
                                                OR = 1.5,
                                                sim = c(500, 250),
                                                mu_delta = 0)

df_dist_variables = df_dist_variables %>% 
  filter((K == 20 & sim == 500) | (K == 5 & sim == 250 & N_ref == 500)) %>%
  mutate(mle_results = glue('data/mle/mle_results_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.Rdata'),
         gamma_results = glue('data/mle-summary/gamma/mle_summary_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.csv'),
         power_results = glue('data/mle-summary/power/mle_summary_power_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.RData'))


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
  # select(sim, K, N_ref, N_foc, gamma, FH, power_statistics, coverage_perc, gamma_hat, mean_true_value) %>% 
  select(K, N_ref, N_foc, gamma, FH, power_statistics, coverage_perc) %>% 
  arrange(K, gamma, desc(N_ref), FH) %>% 
  mutate(across(where(is.numeric),~ round(., 4))) %>% 
  write_csv('output/tables/simulation-results-formatted.csv')
