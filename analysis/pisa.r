# This runs the analysis of PISA data using proposed estimators

rm(list = ls())

source("libraries.R")
source("analysis/Odds.R") # that's where the UTI data comes from
source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
source("analysis/Lfct.r")
source("analysis/DIF.fcts.r")
source("analysis/pisa_analysis.R")
source('analysis/mle-estimate-fns.R')

# analysis_files = list.files('analysis')
# analysis_files = analysis_files[str_detect(analysis_files, '.R')]

df_pisa_aus_science <- read_csv("data/pisa_aus_science_book_1_responses.csv")

set.seed(1)

## DATA SETS
alpha <- 0.05
within.group <- F
zeros <- F
m <- 35

dat <- as.data.frame(df_pisa_aus_science)

dat[, 3:(m + 2)] <- as.matrix(dat[, 3:(m + 2)])

hilf <- dat[, 3:(m + 2)]
hilf[hilf == 2 | hilf == 1] <- 1
hilf[hilf == 8] <- 0
ind <- apply(hilf == 7, 1, mean, na.rm = TRUE) > 0

dat <- dat[!ind, ]
dat[, 3:(m + 2)] <- hilf[!ind, ]

to_purify = dat

dat[, "Total"] <- apply(dat[, 3:(m + 2)], 1, sum)
dat[, "group"] <- bin_data(dat[, "Total"], bins = c(-Inf, seq(5, 30, by = 5), Inf), binType = "explicit")


# Item purification -------------

# Approach 1
#Find DIH items via MH statistic
dif_detection = dichoDif(Data = to_purify[,-1],
                         group = 1,
                         thrTID = 1,
                         # MHstat = 'logOR',
                         focal.name = 'Female',
                         purify = T,
                         method = 'MH')

dif_items = dif_detection$DIFitems

#Calculate on same scale with 1 for all DIF items?
# purified_dat$Total = purified_dat$Total + length(dif_detection$DIFitems)
# purified_dat[, "group"] <- bin_data(purified_dat[, "Total"], bins = c(-Inf, seq(5, 30, by = 5), Inf), binType = "explicit")

# Approach 2

# Get from calculation in normal analysis without purification
pisa_analysis(dat,
              file_path_out = 'output/pisa_normal_analysis.RData',
              m = 35,
              K = 7)

load('output/pisa_results.RData')

dif_items = which(tests$pvalues<0.05)


#Remove these while generating the partition of the ability distribution
purified_dat = to_purify[, -1*(dif_items+2)]

#Calculate the partition of the ability distribution
purified_dat[, "Total"] <- apply(purified_dat[, 3:ncol(purified_dat)], 1, sum)

#Or change bins with 9 DIF items we have one uneven bin
purified_dat[, "group"] <- bin_data(purified_dat[, "Total"], bins = c(-Inf, seq(5, 20, by = 5), Inf), binType = "explicit")

#--------------

# Draw comparison of purification ------

df_purification_comparison = data.frame(purified = purified_dat$Total,
                                        raw = dat$Total)

df_purification_comparison %>% 
  count(purified, raw) %>% 
  group_by(raw) %>% 
  mutate(perc_raw = n/sum(n)) %>% 
  ungroup %>% 
  ggplot() +
  geom_point(aes(x = raw, y = purified, size = perc_raw))

#-----

pisa_analysis(dat,
              file_path_out = 'output/pisa_normal_analysis.RData',
              m = 35,
              K = 7)

pisa_analysis(dat,
              file_path_out = 'output/pisa_reduced_groups.RData',
              m = 35,
              K = 5)

pisa_analysis(purified_dat,
              file_path_out = 'output/pisa_purified_results.RData',
              m = 35 - length(dif_items),
              K = 5)

normal_dat_purified_group = dat
normal_dat_purified_group$group = purified_dat$group

pisa_analysis(normal_dat_purified_group,
              file_path_out = 'output/pisa_purified_grouping.RData',
              m = 35,
              K = 5)

# Perform same LR test for DIF as with simulation study in mle-estimation
# Do the same for sections 5.1/5.2

pisa_mat_to_model_data = function(mat_in){
  df_purified_pisa = mat_in %>% 
    as.data.frame() %>% 
    pivot_longer(cols = starts_with('PS'), 
                 names_to = 'question',
                 values_to = 'response') %>% 
    rename(group = Gender,
           k = group) %>% 
    mutate(group = factor(group, levels = unique(.$group)),
           question = factor(question, levels = unique(.$question)),
           response = as.logical(response))
  
  return(df_purified_pisa)
}

df_purified_pisa = pisa_mat_to_model_data(purified_dat)

model_results = mle_dif_estimate(df_purified_pisa)

full_model = glm(formula = response ~ question:k + question:C(group, contr = contr.SAS(2)) ,
    data = df_purified_pisa, 
    family = 'binomial')

reduced_model = glm(formula = response ~ question:k + question ,
    data = df_purified_pisa, 
    family = 'binomial')

lrtest(full_model, reduced_model)

save(model_results,
     file = 'data/mle/mle_results_pisa.rda')

pisa_deviance = function(df_model_data){
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
  
  return(chisq_test)
}

bt = 1e3

pisa_bootstrap_deviance = map(1:bt,
    function(x){
      sampled_stds = sample(1:nrow(purified_dat),
                            size = nrow(purified_dat),
                            replace = T)
      
      df_model_data = pisa_mat_to_model_data(purified_dat[sampled_stds,])
      
      model_results = pisa_deviance(df_model_data)
      cat(glue('finished {x}/{bt}'))
      
      return(model_results)
    }
) %>% 
  reduce(c)

sum(pisa_bootstrap_deviance > model_results$chisq_test)/length(pisa_bootstrap_deviance)

