# Apply MLE estimator / Global invariance test to demonstrate value of GMH test.

source('libraries.R')
source("analysis/Odds.R") # that's where the UTI data comes from
source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
source("analysis/Lfct.r")
source("analysis/DIF.fcts.r")
source("analysis/pisa_analysis.R")
source('analysis/simulation_study.R')

#Generate simulation data
# debugonce(dif_simulation)
# debugonce(dif_simulate)


m = 10

dif_sims = dif_simulation(Gamma=1,
               m0=10,
               m=50,
               mu.delta=0,
               FH=1,
               K=20,
               N_ref=1000,
               N_foc=200,
               OR=1.5,
               BT=1,
               sim= 100,
               zeros=FALSE,
               within.group=F,
               seed_val=1)

load('data/simulations/Sim_GenDif_m10_K100_Gamma1_OR1.5_sim_results.RData')

dif_sims = simulation_results

# Fit logisitic GLM log(pi_{j|gk} / (1 - pi_{j|gk})) = \beta_{gj} + \mu_{jk}
# Under H_0 : \gamma_1 = ... = \gamma_m = 0
# Where \gamma_j = \beta_{1j} - \beta_{2j} => \beta_{1j} = \beta{2j}
# For g = 1,2
# j = 1, ..., m
# K = 1, ..., K

m = 10

y_sims = map(dif_sims, function(x){x[['y.sim']]}) %>% 
  reduce(rbind)

df_input = y_sims %>% 
  as.data.frame() %>% 
  as_tibble() 

colnames(df_input) = c('k', 'group', paste0('Q', 1:m))

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
      
      df_model_data
    })

sim_cluster = makeCluster(detectCores()-1)
clusterEvalQ(sim_cluster, {
  source("libraries.R")
})

model_results = parLapply(
  sim_cluster,
  df_model_data_split,
  function(df_model_data){
    
    full_model_frame = model.frame(response ~ question:k + question:group - 1,
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
    
    return(chisq_test)
  }) %>% 
    reduce(c)

save(model_results,
     file='test_results.Rdata')

#Under H_0: D_0 - D_1 ~ chisq_m

pchisq(model_results, 2)

sum(pchisq(model_results, 2) < 0.05) / length(model_results)
