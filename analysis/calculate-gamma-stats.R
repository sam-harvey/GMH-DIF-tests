#Find MLE stats for summary table
source('libraries.R')


calculate_mle_summary_stats = function(
  gamma=1,
  m=50,
  FH=1,
  K=20,
  N_ref=1000,
  N_foc=200,
  OR=1.5,
  sims = 100){
  
  # load('data/mle/mle_results_Sim_GenDif_m50_K100_Gamma1_OR1.5_FH0_Nref500_Nfoc500_sims10.RData')
  load(as.character(glue('data/mle/mle_results_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sims}.RData')))
  
  chisq_stats = lapply(model_results,
                       function(x){
                         x$chisq_test
                       }) %>% 
    reduce(c)
  
  # Power of test
  pchisq(chisq_stats, 50)
  sum(pchisq(chisq_stats, 50) < 0.05) / length(chisq_stats)
  
  save(chisq_stats,
       file = glue('data/mle-summary/power/mle_summary_power_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sims}.RData'))
  
  #Find estimated \gamma_j = \beta_1j (we have constraint \beta_2j = 0)
  df_simulation_summary = map(1:length(model_results),
                              function(x){
                                
                                #Find the coefficients corresponding to gamma_j estimates
                                model_matrix = model_results[[x]]$model_matrix
                                model_colnames = model_matrix %>% colnames 
                                gamma_estimate_cols = model_colnames %>% str_detect('C\\(group') %>% which
                                gamma_estimates = model_results[[x]]$full_model_fit$coefficients[gamma_estimate_cols]
                                
                                #Find COV not provided directly by fastLR()
                                fitted_probabilities = model_results[[1]]$full_model_fit$fitted.values
                                fitted_cov = t(model_matrix) %*% (diag(length(fitted_probabilities)) * fitted_probabilities) %*% model_matrix
                                gamma_estimate_col_var = diag(fitted_cov[gamma_estimate_cols,gamma_estimate_cols])
                                
                                fh_values = case_when(FH > 0 ~ c(rep(1.5, FH), rep(0, m-FH)),
                                                      FH == 0 ~ rep(0, m))
                                
                                data.frame(simulation = x,
                                           coefficient = model_colnames[gamma_estimate_cols],
                                           gamma_estimates = gamma_estimates) %>% 
                                  mutate(question = str_extract(coefficient, '(?<=Q)[0-9]+') %>% as.numeric()) %>% 
                                  arrange(question) %>% 
                                  mutate(true_value = fh_values) %>% 
                                  mutate(gamma_estimate_lower = gamma_estimates - 1.96 * sqrt(gamma_estimate_col_var),
                                         gamma_estimate_upper = gamma_estimates + 1.96 * sqrt(gamma_estimate_col_var)) %>% 
                                  mutate(covers_true_value = gamma_estimate_lower <= true_value & true_value <= gamma_estimate_upper) %>% 
                                  as.tbl()
                                
                              }) %>% 
    bind_rows()
  
  file_name = glue('data/mle-summary/gamma/mle_summary_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_sims{sims}.csv')
  write_csv(x=df_simulation_summary,
            file=file_name)
}

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

experiment_params = df_dist_variables[1,]

debugonce(calculate_mle_summary_stats)

walk(1:nrow(df_dist_variables),
     function(x){
       experiment_params = df_dist_variables[x,]
       
       calculate_mle_summary_stats(gamma=experiment_params$gamma,
                                   m=experiment_params$m,
                                   FH=experiment_params$FH,
                                   K=experiment_params$K,
                                   N_ref=experiment_params$N_ref,
                                   N_foc=experiment_params$N_foc,
                                   OR=experiment_params$OR,
                                   sims = experiment_params$sim)
       
       cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
       
       gc()
     })

df_mle_gamma_estimates = map(list.files('data/mle-summary/gamma', full.names = T),
    function(x){
      read_csv(x) %>%
        mutate(file_name = x)
      }) %>% 
  bind_rows()


df_mle_gamma_estimates %>% 
  filter(true_value>0) %>% 
  count(file_name, covers_true_value) %>% 
  group_by(file_name) %>% 
  mutate(perc_n = n/sum(n)) %>% 
  ungroup %>% 
  filter(!covers_true_value)


df_mle_gamma_estimates %>% 
  group_by(file_name, question) %>% 
  summarise(mean(gamma_estimates), mean(true_value))

df_mle_gamma_estimates %>% 
  filter(true_value > 0) %>% 
  group_by(file_name, true_value) %>% 
  summarise(mean(gamma_estimates))


df_mle_gamma_estimates %>% 
  filter(true_value > 0) %>% 
  group_by(true_value) %>% 
  summarise(mean(gamma_estimates))
