source('libraries.R')

df_eap_psg_a = read_excel('data/eap/Reading passage A Animal B Conflict with gender and background info.xlsx')
# df_eap_animal = read_excel('data/eap/Reading passage on Animal with gender info.xlsx')

df_eap_psg_a$student_total = df_eap_psg_a %>% 
  select(contains('item')) %>% 
  apply(1, sum)

#Calculate strata
df_eap_psg_a = df_eap_psg_a %>% 
  mutate(strata = ceiling((student_total - 10) /2) + 1)

# transform to GLM format
eap_to_glm = function(df_in){
  df_in = df_in %>% 
    select(Gender, strata, starts_with('A')) %>% 
    pivot_longer(cols = starts_with('A'), 
                 names_to = 'question',
                 values_to = 'response') %>% 
    rename(k = strata,
           group = Gender)
  
  df_in = df_in  %>% 
    mutate(across(c(group, k), as.factor),
           response = as.logical(response))
  
  return(df_in)
}

eap_deviance_test = function(df_in){
  #Set SAS-style constraint in GLM to \beta_{ref_group=2}j = 0
  #GLM coefficients returned are then \gamma_{focal}j
  full_model_frame = model.frame(response ~ question:k + question:C(group, contr = contr.SAS(2)) - 1,
                                 # response ~ question:k + question:group - 1,
                                 df_in)
  
  full_model_matrix = model.matrix(object = full_model_frame,
                                   data = df_in)
  
  fit_full = fastLR(x=full_model_matrix,
                    y=df_in$response)
  
  reduced_model_frame = model.frame(response ~ question:k + question - 1,
                                    df_in)
  
  reduced_model_matrix = model.matrix(object = reduced_model_frame,
                                      data = df_in)
  
  fit_reduced = fastLR(x=reduced_model_matrix,
                       y=df_in$response)
  
  full_deviance = fit_full$loglikelihood * 2
  reduced_deviance = fit_reduced$loglikelihood * 2
  
  chisq_test = full_deviance - reduced_deviance
  return(chisq_test)
}

df_input = eap_to_glm(df_eap_psg_a)
chisq_stat = eap_deviance_test(df_input)
1 - pchisq(chisq_stat, 15)

# Get LR test value
lrtest(reduced_model, full_model)

# Get boot values
bt = 1e3

bt_test = map(1:bt,
    function(x){
      bt_sample = sample(1:nrow(df_eap_psg_a),
             nrow(df_eap_psg_a),
             replace = T)
      df_sample = df_eap_psg_a[bt_sample,]
      df_sample = eap_to_glm(df_sample)
      
      eap_deviance_test(df_sample)
    })

sum(bt_test > chisq_stat)/length(bt_test)
