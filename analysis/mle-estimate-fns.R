#' Title
#'
#' @param simulation_results 
#'
#' @return
#' @export
#'
#' @examples
simulation_results_to_glm_input = function(simulation_results,
                                           split=T,
                                           m=50){
  
  y_sims = map(1:length(simulation_results),
               function(x){
                 df_sim_results = simulation_results[[x]]$'y.sim' %>% 
                   as.data.frame()
                   
                   colnames(df_sim_results) = c('k', 'group', paste0('Q', 1:m))
                   
                   df_sim_results = df_sim_results %>% 
                     mutate(sample_label = x) %>% 
                     dplyr::select(sample_label, everything())
                   
                   return(df_sim_results)
                 }) %>% 
    bind_rows()
  
  df_input = y_sims %>% 
    as_tibble() 
  
  df_input = df_input %>% 
    pivot_longer(cols = starts_with('Q'), 
                 names_to = 'question',
                 values_to = 'response')
  
  df_input$group = factor(df_input$group, levels = unique(df_input$group))
  df_input$k = factor(df_input$k, levels = unique(df_input$k))
  df_input$response = as.logical(df_input$response)
  
  if(split){
    df_model_return = split(df_input, df_input$sample_label)
  } else{
    df_model_return = df_input
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
mle_dif_estimate = function(df_model_data, 
                            method = c('fast', 'base'),
                            full_model_formula = response ~ question:k + question:C(group, contr = contr.SAS(2)) - 1,
                            reduced_model_formula = response ~ question:k + question - 1){
  
  if(method == 'fast'){
  #Set SAS-style constraint in GLM to \beta_{ref_group=2}j = 0
  #GLM coefficients returned are then \gamma_{focal}j
  full_model_frame = model.frame(full_model_formula,
                                 df_model_data)
  
  full_model_matrix = model.matrix(object = full_model_frame,
                                   data = df_model_data)
  
  fit_full = fastLR(x=full_model_matrix,
                    y=df_model_data$response)
  
  reduced_model_frame = model.frame(reduced_model_formula,
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
  } else if(method == 'base'){
    
    fit_full = glm(data=df_model_data,
                   formula = full_model_formula,
                   family = 'binomial')
    
    fit_reduced = glm(data=df_model_data,
                      formula = reduced_model_formula,
                      family = 'binomial')
    
    full_deviance = fit_full$deviance
    reduced_deviance = fit_reduced$deviance
    
    chisq_test = full_deviance - reduced_deviance
    
    results = list(model_matrix = as(model.matrix(fit_full), 'dgCMatrix'),
                   full_model_fit = fit_full,
                   chisq_test = chisq_test)
  }
  
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
                          split=T,
                          method = c('fast', 'base'),
                          full_model_formula = response ~ question:k + question:C(group, contr = contr.SAS(2)) - 1,
                          reduced_model_formula = response ~ question:k + question - 1,
                          m=50){
  load(sim_results_path)
  
  df_model_data_split = simulation_results_to_glm_input(simulation_results,
                                                        split = split,
                                                        m=m)
  
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
        mle_dif_estimate(df_model_data, 
                         method = method,
                         full_model_formula = full_model_formula,
                         reduced_model_formula = reduced_model_formula)
      })
    
  } else{
    
    model_results = lapply(
      df_model_data_split,
      function(df_model_data){
        mle_dif_estimate(df_model_data, 
                         method = method,
                         full_model_formula = full_model_formula,
                         reduced_model_formula = reduced_model_formula)
      })
    
  }
  
  save(model_results,
       file = output_path)
}

#https://github.com/cran/LRQMM/blob/master/R/spginv.R
spginv<-function (x)
{
  Xsvd <- sparsesvd::sparsesvd(x)
  Positive <- Xsvd$d > max(sqrt(.Machine$double.eps) * Xsvd$d[1L], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(x)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
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
  sims = 100,
  sim_name = glue('data/mle/mle_results_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.RData'),
  output_file_name = glue('data/mle-summary/gamma/mle_summary_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.csv')){
  
  # load('data/mle/mle_results_Sim_GenDif_m50_K100_Gamma1_OR1.5_FH0_Nref500_Nfoc500_sims10.RData')
  load(as.character(sim_name))
  
  chisq_stats = lapply(model_results,
                       function(x){
                         x$chisq_test
                       }) %>% 
    reduce(c)
  
  # Power of test
  # pchisq(chisq_stats, 50)
  # sum(pchisq(chisq_stats, 50) < 0.05) / length(chisq_stats)
  
  save(chisq_stats,
       file = glue('data/mle-summary/power/mle_summary_power_Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sims}.RData'))
  
  #Find estimated \gamma_j = \beta_1j (we have constraint \beta_2j = 0)
  df_simulation_summary = map(1:length(model_results),
                              possibly(function(x){
                                
                                #Find the coefficients corresponding to gamma_j estimates
                                model_matrix = model_results[[x]]$model_matrix
                                model_colnames = model_matrix %>% colnames 
                                gamma_estimate_cols = model_colnames %>% str_detect('C\\(group') %>% which
                                gamma_estimates = model_results[[x]]$full_model_fit$coefficients[gamma_estimate_cols]
                                
                                #Find COV not provided directly by fastLR()
                                fitted_probabilities = model_results[[x]]$full_model_fit$fitted.values
                                sparse_id_mat = sparse_identity(length(fitted_probabilities))
                                V = sparse_id_mat * fitted_probabilities * (1 - fitted_probabilities)
                                fitted_cov = t(model_matrix) %*% (V) %*% model_matrix
                                #TODO try catch etc
                                fitted_cov = solve(fitted_cov)
                                # fitted_cov = spginv(fitted_cov)
                                gamma_estimate_col_var = diag(fitted_cov[gamma_estimate_cols,gamma_estimate_cols])
                                
                                # fh_values = case_when(FH > 0 ~ c(rep(1.5, FH), rep(0, m-FH)),
                                #                       FH == 0 ~ rep(0, m))
                                
                                # These are on \gamma scale see start of section 3
                                fh_values = case_when(FH > 0 ~ c(rep(-log(OR), FH), rep(0, m-FH)),
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
                                
                              }, otherwise = NULL)) %>% 
    bind_rows()
  
  write_csv(x=df_simulation_summary,
            file=output_file_name)
}

sparse_identity = function(n){
  bandSparse(n, n, 0, list(rep(1, n+1))) 
}

setup_data_dirs = function(){
  dirs = c("data",
    "data/animal",
    "data/booklet",
    "data/data-dictionaries",
    "data/eap",
    "data/joint-distributions",
    "data/mle",
    "data/mle-summary",
    "data/mle-summary/gamma",
    "data/mle-summary/gamma/archive",
    "data/mle-summary/power",
    "data/mle-summary/power/archive",
    "data/mle/archive",
    "data/New folder",
    "data/raw-responses",
    "data/responses",
    "data/simulations")
  
  walk(dirs, ~dir.create(.))
}

run_simulations = function(N_ref = 1000,
                           N_foc = 200,
                           K=5, 
                           m=50,
                           m0 = 10,
                           FH=5,
                           OR=1.5,
                           sim= 1e1,
                           gamma = 1,
                           mu_delta = 0,
                           simulation_file,
                           skip_simulation_stage=F,
                           split = F,
                           parallel=T){
  
  experiment_params = create_simulation_scenarios(N_ref = N_ref,
                                                  K=K, 
                                                  m=m,
                                                  m0 = m0,
                                                  FH=FH,
                                                  OR=OR,
                                                  sim=sim,
                                                  gamma = gamma,
                                                  mu_delta = mu_delta) %>% 
    mutate(N_foc = N_ref)
  
  if(!skip_simulation_stage){
    dif_sims = dif_simulation(
      Gamma=gamma,
      m0=m0,
      m=m,
      mu.delta=mu_delta,
      FH=FH,
      K=K,
      N_ref=N_ref,
      N_foc=N_foc,
      OR=OR,
      sim=sim,
      zeros=FALSE,
      within.group=F,
      seed_val=1,
      use_bt=F
    )
  }
  
  mle_estimation(sim_results_path=simulation_file,
                 parallel=parallel,
                 split=split,
                 m=m,
                 method = 'fast'
                 # method = 'base'
  )
  
  calculate_mle_summary_stats(gamma=gamma,
                              m=m,
                              FH=FH,
                              K=K,
                              N_ref=N_ref,
                              N_foc=N_foc,
                              mu_delta=mu_delta,
                              OR=OR,
                              sims=sim)
  
  
  cat(glue('completed simulation: {x}/{nrow(df_dist_variables)}\n'))
  
}
