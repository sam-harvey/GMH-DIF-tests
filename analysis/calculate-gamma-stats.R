#Find MLE stats for summary table
source('libraries.R')
source('analysis/mle-estimate-fns.R')

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
mu_delta = c(0,1)

df_dist_variables = crossing(K, m, m0, mu_delta, gamma, FH, N_ref, OR, sim) %>% 
  mutate(N_foc = case_when(N_ref == 1000 ~ 200,
                           N_ref == 500 ~ 500)) %>% 
  mutate(simulation_file = glue("data/mle/mle_results_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.RData") %>% 
           as.character())


# debugonce(calculate_mle_summary_stats)

walk(
  1:nrow(df_dist_variables),
  # 61,
  # 21,
     function(x){
       experiment_params = df_dist_variables[x,]
       
       calculate_mle_summary_stats(gamma=experiment_params$gamma,
                                   m=experiment_params$m,
                                   FH=experiment_params$FH,
                                   K=experiment_params$K,
                                   N_ref=experiment_params$N_ref,
                                   N_foc=experiment_params$N_foc,
                                   mu_delta=experiment_params$mu_delta,
                                   OR=experiment_params$OR,
                                   sims = experiment_params$sim)
       
       cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
       
       gc()
     })

gamma_estimate_files = list.files('data/mle-summary/gamma', full.names = T, include.dirs = FALSE)
gamma_estimate_files = gamma_estimate_files[!dir.exists(gamma_estimate_files)]

df_mle_gamma_estimates = map(gamma_estimate_files,
    function(x){
      read_csv(x) %>%
        mutate(file_name = x)
      }) %>% 
  bind_rows()

df_mle_gamma_estimates %>% 
  filter(true_value<0) %>% 
  count(file_name, covers_true_value) %>% 
  group_by(file_name) %>% 
  mutate(perc_n = n/sum(n)) %>% 
  ungroup %>% 
  filter(covers_true_value)

df_mle_gamma_estimates %>% 
  filter(true_value<0) %>% 
  count(covers_true_value) %>% 
  ungroup %>% 
  mutate(perc_n = n/sum(n)) %>% 
  filter(covers_true_value)

df_mle_gamma_estimates %>% 
  group_by(file_name, question) %>% 
  summarise(mean(gamma_estimates), mean(true_value))

df_mle_gamma_estimates %>% 
  filter(true_value < 0) %>% 
  group_by(file_name, true_value) %>% 
  summarise(mean(gamma_estimates))

df_mle_gamma_estimates %>% 
  # filter(true_value < 0) %>% 
  group_by(true_value) %>%
  summarise(mean(gamma_estimates))

# Join to dist variables to summarise in tables

df_dist_variables %>% 
  mutate(simulation = basename(simulation_file)) %>% 
  mutate(simulation = str_extract(simulation, '(?<=Sim).*(?=\\.)')) %>% 
  inner_join(df_mle_gamma_estimates %>% 
               mutate(simulation = basename(file_name)) %>% 
               mutate(simulation = str_extract(simulation, '(?<=Sim).*(?=\\.)'))) %>% 
  select(K, N_ref, N_foc, gamma, FH, mu_delta, simulation, coefficient, true_value, gamma_estimates, gamma_estimate_lower, gamma_estimate_upper, covers_true_value) %>% 
  group_by(mu_delta, K, N_ref, N_foc, gamma, FH) %>% 
  summarise(covers_true_value_perc = sum(covers_true_value)/n(),
            covers_true_value_th = sum(covers_true_value * (true_value == 0)) / sum(true_value == 0),
            covers_true_value_fh = sum(covers_true_value * (true_value != 0)) / sum(true_value != 0),
            average_fh_coefficient = weighted.mean(gamma_estimates, true_value != 0),
            average_th_coefficient = weighted.mean(gamma_estimates, true_value == 0)) %>% 
  ungroup() %>% 
  arrange(mu_delta, gamma, K, desc(N_ref), FH) %>% 
  write_csv('output/tables/simulation-gamma-stats-formatted.csv')


# Do similar thing for deviance stats
power_estimate_files = list.files('data/mle-summary/power', full.names = T, include.dirs = FALSE)
power_estimate_files = power_estimate_files[!dir.exists(power_estimate_files)]

df_mle_gamma_estimates = map(power_estimate_files,
                             function(x){
                               load(x)
                               chisq_stats
                             }) %>% 
  reduce(c)
