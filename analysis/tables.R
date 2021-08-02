rm(list = ls())
source('libraries.R')

load('output/pisa_results.RData')
df_pisa_aus = read_csv('data/pisa_aus_science_book_1_responses.csv')

#Create individual question DIF test values
questions = df_pisa_aus %>% colnames
questions = questions[3:length(questions)]

df_output_table = data.frame(estimator = c(L),
                             se_estimator = c(sqrt(VarL)),
                             pvalues = c(tests$pvalues))
                             
df_output_table = df_output_table %>% 
  mutate(question = questions)

df_output_table %>% 
  write_csv('output/tables/pisa-estimate-table.csv')

df_output_table %>% 
  mutate_if(is.numeric, ~round(., 4)) %>% 
  write_csv('output/tables/pisa-estimate-table-rounded.csv')

#Create OR pair comparison
questions = df_pisa_aus %>% colnames
questions = questions[3:length(questions)]

colnames(ORpair.matrix) = questions
rownames(ORpair.matrix) = questions

colnames(pval.ORpair.matrix)  = questions
rownames(pval.ORpair.matrix)  = questions

#Comparing fraction of questions pairs with 5% significant ORs 
#for items from the same parent question to all items
locations = which(pval.ORpair.matrix < 0.05 & upper.tri(pval.ORpair.matrix))

cols = floor(locations /35) + 1
rows = locations %% 35

data.frame(q1 = questions[rows],
           q2 = questions[cols]) %>% 
  mutate(q1_parent_q = str_extract(q1, "PS[0-9]+"),
         q2_parent_q = str_extract(q2, "PS[0-9]+")) %>% 
  mutate(test = q1_parent_q == q2_parent_q) %>% 
  count(test)

combn(questions, 2) %>%
  t %>% 
  data.frame %>% 
  mutate(q1_parent_q = str_extract(X1, "PS[0-9]+"),
         q2_parent_q = str_extract(X2, "PS[0-9]+"))%>% 
  mutate(test = q1_parent_q == q2_parent_q) %>% 
  count(test)

#Pairwise test for items with DIF detected by individual test
question_index = which(tests$pvalues < 0.05)

lower_cv_tbl = tests$lowerdiffL[question_index, question_index] %>% t %>% round(4)
lower_cv_tbl[is.na(lower_cv_tbl)] = 0

upper_cv_tbl = tests$upperdiffL[question_index, question_index] %>% round(4)
upper_cv_tbl[is.na(upper_cv_tbl)] = 0

cv_tbl = upper_cv_tbl + lower_cv_tbl

question_labels = questions[question_index]

rownames(cv_tbl) = question_labels
colnames(cv_tbl) = question_labels
diag(cv_tbl) = '-'

cv_tbl = cv_tbl %>% 
  as.data.frame() %>% 
  mutate(question = row.names(.)) %>% 
  dplyr::select(question, everything()) 

cv_tbl %>% 
  write_csv('output/tables/dif-pairwise-dif-example.csv')

ci_both_ge_zero = (t(lower_cv_tbl > 0) & (upper_cv_tbl > 0))[upper.tri(upper_cv_tbl)]
ci_both_le_zero = (t(lower_cv_tbl < 0) & (upper_cv_tbl < 0))[upper.tri(upper_cv_tbl)]

#All questions with significant detected
ci_both_ge_zero | ci_both_le_zero

#Table of OR pair comparison
n_items = 7

question_labels = df_output_table$question[1:n_items]

or_tbl = ORpair.matrix[1:n_items, 1:n_items] %>% round(4)
colnames(or_tbl) = question_labels

or_tbl[lower.tri(or_tbl)] = '-'
diag(or_tbl) = '-'

tbl_output = or_tbl[-n_items,-1] %>% 
  as.data.frame() %>% 
  mutate(question = question_labels[1:(n_items-1)]) %>% 
  dplyr::select(question, everything())

tbl_output %>% 
  write_csv('output/tables/or-pair-output.csv')

#Table of method comparison
data.frame(Method = c('Asym', 'Boot'),
           W = rep(tests$W, 2),
           `W p-value` = Wpvals[c('asymp', 'BT')],
           `W ind` = rep(tests$Wind, 2),
           `W ind p-values` = Wpvals[c('ind-asymp', 'ind-BT')]) %>% 
  write_csv('output/tables/output-table-pisa-W-comparison.csv')