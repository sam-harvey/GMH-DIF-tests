#' Title
#'
#' @param simulation_results 
#'
#' @return
#' @export
#'
#' @examples
simulation_results_to_glm_input = function(simulation_results,
                                           split=T){
  
  y_sims = map(simulation_results, function(x){x[['y.sim']]}) %>% 
    reduce(rbind)
  
  colnames(y_sims) = c('k', 'group', paste0('Q', 1:m))
  
  df_input = y_sims %>% 
    as.data.frame() %>% 
    as_tibble() 
  
  df_input = df_input %>% 
    pivot_longer(cols = starts_with('Q'), 
                 names_to = 'question',
                 values_to = 'response')
  
  df_input = df_input  %>% 
    mutate(across(c(group, k), as.factor),
           response = as.logical(response))
  
  df_glm_input = df_input %>% 
    group_by(k, group, question) %>% 
    #Create a label of the N_k samples
    mutate(sample_label = row_number()) %>% 
    ungroup()
  
  if(split){
    df_model_return = map(unique(df_glm_input$sample_label),
                            function(x){
                              df_model_data = df_glm_input %>% 
                                filter(sample_label == x)
                              
                              return(df_model_data)
                            })
  } else{
    df_model_return = df_glm_input
  }
  
  return(df_model_return)
}

#' Title
#'
#' @param df_model_data 
#'
#' @return
#' @export
#'
#' @examples
mle_dif_estimate = function(df_model_data){
  #Set SAS-style constraint in GLM to \beta_{ref_group=2}j = 0
  #GLM coefficients returned are then \gamma_{focal}j
  full_model_frame = model.frame(response ~ question:k + question:C(group, contr = contr.SAS(2)) - 1,
                                 # response ~ question:k + question:group - 1,
                                 df_model_data)
  
  full_model_matrix = model.matrix(object = full_model_frame,
                                   data = df_model_data)
  
  fit_full = fastLR(x=full_model_matrix,
                    y=df_model_data$response)
  
  reduced_model_frame = model.frame(response ~ question:k + question - 1,
                                    df_model_data)
  
  reduced_model_matrix = model.matrix(object = reduced_model_frame,
                                      data = df_model_data)
  
  fit_reduced = fastLR(x=reduced_model_matrix,
                       y=df_model_data$response)
  
  full_deviance = fit_full$loglikelihood * 2
  reduced_deviance = fit_reduced$loglikelihood * 2
  
  chisq_test = full_deviance - reduced_deviance
  
  results = list(model_matrix = as(full_model_matrix, 'dgCMatrix'),
                 full_model_fit = fit_full,
                 chisq_test = chisq_test)
  
  return(results)
}

#' Title
#'
#' @param sim_results_path 
#' @param output_path 
#'
#' @return
#' @export
#'
#' @examples
mle_estimation = function(sim_results_path = 'data/simulations/Sim_GenDif_m50_K20_Gamma1_OR1.5_FH0_Nref500_Nfoc500_sims10000.RData',
                          output_path = glue('data/mle/mle_results_{basename(sim_results_path)}'),
                          parallel=F,
                          split=T){
  load(sim_results_path)
  
  df_model_data_split = simulation_results_to_glm_input(simulation_results,
                                                        split = split)
  
  if(!split){
    df_model_data_split = list(df_model_data_split)
  }
  
  if(parallel){
    sim_cluster = makeCluster(detectCores()-1)
    clusterEvalQ(sim_cluster, {
      source("libraries.R", local = TRUE)
      source("analysis/mle-estimate-fns.R", local = TRUE)
    })
    
    model_results = parLapply(
      sim_cluster,
      df_model_data_split,
      function(df_model_data){
        mle_dif_estimate(df_model_data)
      })
    
  } else{
    
    model_results = lapply(
      df_model_data_split,
      function(df_model_data){
        mle_dif_estimate(df_model_data)
      })
    
  }
  
  
  save(model_results,
       file = output_path)
}

calculate_mle_summary_stats = function(
  gamma=1,
  m=50,
  FH=1,
  K=20,
  N_ref=1000,
  N_foc=200,
  mu_delta=0,
  OR=1.5,
  sims = 100){
  
  # load('data/mle/mle_results_Sim_GenDif_m50_K100_Gamma1_OR1.5_FH0_Nref500_Nfoc500_sims10.RData')
  load(as.character(glue('data/mle/mle_results_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.RData')))
  
  chisq_stats = lapply(model_results,
                       function(x){
                         x$chisq_test
                       }) %>% 
    reduce(c)
  
  # Power of test
  pchisq(chisq_stats, 50)
  sum(pchisq(chisq_stats, 50) < 0.05) / length(chisq_stats)
  
  save(chisq_stats,
       file = glue('data/mle-summary/power/mle_summary_power_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.RData'))
  
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
                                sparse_id_mat = sparse_identity(length(fitted_probabilities))
                                V = sparse_id_mat * fitted_probabilities * (1 - fitted_probabilities)
                                fitted_cov = t(model_matrix) %*% (V) %*% model_matrix
                                fitted_cov = solve(fitted_cov)
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
  
  file_name = glue('data/mle-summary/gamma/mle_summary_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.csv')
  write_csv(x=df_simulation_summary,
            file=file_name)
}

sparse_identity = function(n){
  bandSparse(n, n, 0, list(rep(1, n+1))) 
}


