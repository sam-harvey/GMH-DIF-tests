rm(list = ls())
source('libraries.R')

load('output/pisa_purified_results.RData')
df_pisa_aus = read_csv('data/pisa_aus_science_book_1_responses.csv')

data.frame(Method = c('Asym', 'Boot'),
           W = rep(tests$W, 2),
           `W p-value` = Wpvals[c('asymp', 'BT')],
           `W ind` = rep(tests$Wind, 2),
           `W ind p-values` = Wpvals[c('ind-asymp', 'ind-BT')])

load('output/pisa_purified_grouping.RData')

data.frame(Method = c('Asym', 'Boot'),
           W = rep(tests$W, 2),
           `W p-value` = Wpvals[c('asymp', 'BT')],
           `W ind` = rep(tests$Wind, 2),
           `W ind p-values` = Wpvals[c('ind-asymp', 'ind-BT')])


load('output/pisa_reduced_groups.RData')

data.frame(Method = c('Asym', 'Boot'),
           W = rep(tests$W, 2),
           `W p-value` = Wpvals[c('asymp', 'BT')],
           `W ind` = rep(tests$Wind, 2),
           `W ind p-values` = Wpvals[c('ind-asymp', 'ind-BT')])
