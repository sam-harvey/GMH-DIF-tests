#' Title
#'
#' @param simulation_results 
#'
#' @return
#' @export
#'
#' @examples
simulation_results_to_glm_input = function(simulation_results){
  
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
  
  df_model_data_split = map(unique(df_glm_input$sample_label),
                            function(x){
                              df_model_data = df_glm_input %>% 
                                filter(sample_label == x)
                              
                              return(df_model_data)
                            })
  
  return(df_model_data_split)
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
  
  full_deviance = fit_full$loglikelihood * -2
  reduced_deviance = fit_reduced$loglikelihood * -2
  
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
                          m = 10,
                          parallel=F){
  load(sim_results_path)
  
  df_model_data_split = simulation_results_to_glm_input(simulation_results)
  
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