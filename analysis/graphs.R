rm(list = ls())
source('libraries.R')
loadfonts(device = "win")
par(family = "LM Roman 10")

df_tbl_1 = read_csv('data/table-1-data.csv', skip=1)
df_tbl_2 = read_csv('data/table-2-data.csv', skip=1)

plot_fn = function(df){
  df %>% 
    select(K, N1, N2, gamma = 4, FH , W_asym = 6) %>% 
    mutate(plot_label = glue("K = {K}, N1 = {N1}, N2 = {N2}, Gamma = {gamma}"),
           plot_axis = factor(FH, levels = unique(.$FH))) %>% 
    ggplot(aes(plot_axis, W_asym, group = plot_label))+
    facet_wrap(~plot_label)+
    geom_line()+
    geom_point()+
    scale_y_continuous(limits = c(0,1))+
    theme_bw()+
    labs(title = 'Power of asymptotic test statistic W_asym',
         y = 'test',
         x = 'sadf')+
    theme(text=element_text(size=16,  family="LM Roman 10"))
}

plot_fn(df_tbl_1)
plot_fn(df_tbl_2)

plot_df = df_tbl_1 %>% 
  rename(gamma = 4) %>% 
  pivot_longer(cols = Asym:Boot_3) %>% 
  mutate(name_position = str_extract(name, '[0-9]+')) %>% 
  mutate(facet_name = case_when(name_position %>% is.na() ~ 'Power~italic(W)',
                                name_position == 1 ~ 'Power~italic(W)[ind]',
                                name_position == 2 ~ 'Coverage~hat(gamma)[j]',
                                name_position == 3 ~ 'Coverage~hat(gamma)[j]~-~hat(gamma)[l]'),
         name_type = str_extract(name, '^[^_]+'))%>% 
  mutate(plot_label = glue("{{K=={K}}}*', '*{{N[1]=={N1}}}*', '*{{N[2]=={N2}}}*', '*{{Gamma=={gamma}}}") %>% as.character(),
         plot_axis = factor(FH, levels = unique(.$FH)))

label_parse <- function(breaks) {
  parse(text = breaks)
}

plot_fn_2 = function(df){
  legend_guide = guide_legend('title', nrow = 4)
  
  df  %>% 
    ggplot(aes(x = plot_axis 
               ,y = value
               # ,colour = plot_label
               ,linetype = plot_label
               ,shape = plot_label
               ,group = plot_label
               ))+
    # facet_grid(name_type~facet_name, scales = 'free_y')+
    facet_grid(name_type~facet_name, 
               scales = 'free_y',
               labeller = labeller(facet_name = label_parsed))+
    geom_line()+
    geom_point()+
    # scale_y_continuous(limits = c(0,1))+
    theme_bw()+
    labs(title = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind]),
         y = NULL,
         x = "FH"
         # ,colour = 'Plot series',
         # linetype = 'Plot series',
         # shape = 'Plot series'
         )+
    theme(text=element_text(size=16,  family="LM Roman 10"),
          legend.position = 'bottom',
          # legend.direction = 'vertical',
          # legend.key.height = unit(0.2, 'npc'),
          legend.text=element_text(size=8))+
    scale_color_manual('title', values=viridisLite::viridis(8), label = label_parse)+
    scale_fill_manual('title', values=viridisLite::viridis(8), label = label_parse)+
    scale_linetype_manual('title', values=seq(0,7), label = label_parse)+
    # scale_shap_discrete(label = label_parse)+label = label_parse
    scale_shape_manual(values=seq(0,7), label = label_parse)+
    guides(
      colour = guide_legend('title', nrow = 4),
      linetype = guide_legend('title', nrow = 4),
      shape = guide_legend('title', nrow = 4)
      # ,fill = guide_legend(override.aes = list(shape = 21))
      )
}

# tikz('power_statistics.tex', width = 6,height = 3.5,pointsize = 12) #define plot name size and font size
# 
# p = plot_fn_2(plot_df)
# p
# 
# dev.off()
# 
# ggsave('table-1-plot.png', p, h = 9, w = 12)
# 
# tikz('power_statistics.tex', width = 6,height = 3.5,pointsize = 12) #define plot name size and font size

plot_df %>% 
  # head(20) %>%
  filter(str_detect(facet_name, 'Power')) %>%
  # mutate(across(where(is.character), ~paste0("$", ., "$"))) %>% 
  plot_fn_2()



#Working
plot_df %>% 
  # head(50) %>% 
  ggplot(aes(x = plot_axis 
             ,y = value
             # ,colour = plot_label
             # ,fill = plot_label
             ,linetype = plot_label
             ,shape = plot_label
             ,group = plot_label
  ))+
  facet_grid(name_type~facet_name, 
             scales = 'free_y',
             labeller = labeller(facet_name = label_parsed))+
  geom_line()+
  geom_point()+
  theme_bw()+
  # theme(text=element_text(size=16,  family="LM Roman 10"),
  #       legend.position = 'bottom',
  #       legend.text=element_text(size=8))+
  # scale_color_manual(values=viridisLite::viridis(8), label = label_parse)+
  # scale_fill_manual(values=viridisLite::viridis(8), label = label_parse)+
  # scale_linetype_manual(values=seq(0,7), label = label_parse)+
  # scale_shape_manual(values=seq(0,7), label = label_parse)
  labs(
    title = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind]),
    y = NULL,
    x = "FH"
    # ,colour = 'Plot series',
    ,linetype = 'Plot series'
    ,shape = 'Plot series'
  )

#Working
plot_df %>% 
  # head(50) %>% 
  ggplot(aes(x = plot_axis 
             ,y = value
             ,colour = plot_label
             # ,fill = plot_label
             ,linetype = plot_label
             ,shape = plot_label
             ,group = plot_label
  ))+
  facet_grid(name_type~facet_name, 
             scales = 'free_y',
             labeller = labeller(facet_name = label_parsed))+
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(text=element_text(size=16,  family="LM Roman 10"),
        legend.position = 'bottom',
        legend.text=element_text(size=8))


plot_df = plot_df %>% 
  arrange(plot_label)

plot_labels = plot_df$plot_label %>% unique %>% sort
plot_labels
map(plot_labels, function(x){
  bquote(x)
})

plot_labs = c(expression({K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}),
expression({K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}) ,
expression({K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}),
expression({K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}) ,
expression({K==20}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}) ,
expression({K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}),
expression({K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}),
expression({K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}))



plot_labs = c(
expression(                {K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}    ),
expression(                {K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}   ),
expression(                {K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}      ),
expression(                {K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}     ),
expression(                {K==20}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}      ),
expression(                {K==20}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}     ),
expression(                {K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}       ),
expression(                {K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}     ))


plot_df %>% 
  # head(50) %>% 
  ggplot(aes(x = plot_axis 
             ,y = value
             ,colour = plot_label
             # ,fill = plot_label
             ,linetype = plot_label
             ,shape = plot_label
             ,group = plot_label
  ))+
  facet_grid(name_type~facet_name, 
             scales = 'free_y',
             labeller = labeller(facet_name = label_parsed))+
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(text=element_text(size=16,  family="LM Roman 10"),
        legend.position = 'bottom',
        legend.text=element_text(size=8))+
  scale_color_manual(values=viridisLite::viridis(8), breaks = plot_labels, labels = plot_labs)+
  scale_linetype_manual(values=seq(1,8), breaks = plot_labels, labels = plot_labs)+
  scale_shape_manual(values=seq(1,8), breaks = plot_labels, labels = plot_labs)+
  labs(
    title = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind])
    ,y = NULL
    ,x = "FH"
    ,colour = NULL
    ,linetype = NULL
    ,shape = NULL
  )


### WORKING ----------

plot_breaks = plot_df$plot_label %>% unique %>% sort

plot_labs_2 = c()

for(i in plot_breaks){
  plot_labs_2 = c(plot_labs_2, parse(text=i))
}

# plot_labs = c(
#   expression(                {K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}    ),
#   expression(                {K==100}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}   ),
#   expression(                {K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}      ),
#   expression(                {K==100}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}     ),
#   expression(                {K==20}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==1}      ),
#   expression(                {K==20}*', '*{N[1]==1000}*', '*{N[2]==200}*', '*{Gamma==10}     ),
#   expression(                {K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==1}       ),
#   expression(                {K==20}*', '*{N[1]==500}*', '*{N[2]==500}*', '*{Gamma==10}     )
# )

plot_labs = plot_labs_2

# plot_labels_parsed = map(plot_labs, as.expression) %>% reduce(function(x, y){x%c%y})

# "%c%" <- function(x, y) parse( text=paste(x, ",", y) )

reduce()

plot_function = function(df,
                         plot_breaks, 
                         plot_labels_parsed, 
                         legend_position = 'bottom',
                         title_value = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind]),
                         x_value = 'FH'){
  df %>% 
    ggplot(aes(x = plot_axis 
               ,y = value
               ,colour = plot_label
               # ,fill = plot_label
               # ,linetype = plot_label
               ,shape = plot_label
               ,group = plot_label
    ))+
    facet_grid(name_type~facet_name, 
               scales = 'free_y',
               labeller = labeller(facet_name = label_parsed))+
    geom_line()+
    geom_point()+
    theme_bw()+
    theme(text=element_text(size=16,  family="LM Roman 10"),
          legend.position = legend_position,
          legend.text=element_text(size=8))+
    # scale_color_manual(values=viridisLite::viridis(8), breaks = plot_breaks, labels = plot_labels_parsed)+
    scale_color_manual(values=brewer_pal(type = 'qual', palette = 2)(8), breaks = plot_breaks, labels = plot_labels_parsed)+
    # scale_linetype_manual(values=seq(1,8), breaks = plot_breaks, labels = plot_labels_parsed)+
    scale_shape_manual(values=seq(1,8), breaks = plot_breaks, labels = plot_labels_parsed)+
    labs(
      title = title_value
      ,y = NULL
      ,x = x_value
      ,colour = NULL
      ,linetype = NULL
      ,shape = NULL
    )
}

plot_function(plot_df, plot_breaks, plot_labels_parsed)

#Full plot
plot_df %>% 
  # head(50) %>% 
  ggplot(aes(x = plot_axis 
             ,y = value
             ,colour = plot_label
             # ,fill = plot_label
             ,linetype = plot_label
             ,shape = plot_label
             ,group = plot_label
  ))+
  facet_grid(name_type~facet_name, 
             scales = 'free_y',
             labeller = labeller(facet_name = label_parsed))+
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(text=element_text(size=16,  family="LM Roman 10"),
        legend.position = 'bottom',
        legend.text=element_text(size=8))+
  # scale_color_manual(values=viridisLite::viridis(8), breaks = plot_labels, labels = plot_labs)+
  scale_color_manual(values=brewer_pal(type = 'qual')(8), breaks = plot_labels, labels = plot_labs)+
  scale_linetype_manual(values=seq(1,8), breaks = plot_labels, labels = plot_labs)+
  scale_shape_manual(values=seq(1,8), breaks = plot_labels, labels = plot_labs)+
  labs(
    title = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind])
    ,y = NULL
    ,x = "FH"
    ,colour = NULL
    ,linetype = NULL
    ,shape = NULL
  )



df_power_stats = plot_df %>% 
  filter(str_detect(facet_name, 'Power'))

p1 = plot_function(df_power_stats, plot_breaks, plot_labels_parsed, legend_position = 'none', x_value = NULL)

df_coverage_stats = plot_df %>% 
  filter(str_detect(facet_name, 'Coverage'))

p2 = plot_function(df_coverage_stats, plot_breaks, plot_labels_parsed, title_value = NULL)

p3 = p1/p2
p3

###### ---------

library(scales)


plot_df %>% 
  # head(50) %>% 
  ggplot(aes(x = plot_axis 
             ,y = value
             ,colour = plot_label
             # ,fill = plot_label
             # ,linetype = plot_label
             ,shape = plot_label
             ,group = plot_label
  )) +
  geom_point()+
  # geom_line()+
  scale_colour_manual(name = "F", 
                     labels = parse_format(),
                     values = 1:8) +
  scale_shape_manual(name = "F", 
                     labels = parse_format(),
                     values = 1:8)

ggplot(data = mtcars, aes(x = mpg, 
                          y = cyl, 
                          color = factor(gear),
                          linetype = factor(gear),
                          shape = factor(gear),
                          group = factor(gear))) +
  geom_point() +
  geom_line()+
  scale_color_manual(name = "F", 
                     labels = parse_format(),
                     values = c("red", "blue", "green")) +
  scale_shape_manual(name = "F", 
                     labels = parse_format(),
                     values = c(16, 15, 18))+
  scale_linetype_manual(name = "F", 
                     labels = parse_format(),
                     values = c(16, 15, 18))

# p = plot_fn_2(plot_df %>% 
#                 filter(str_detect(facet_name, 'Power')) %>% 
#                 mutate(across(where(is.character), ~glue("${.}$"))))
# p
# 
# dev.off()
# 
# tikz('mtcars.tex', width = 6,height = 3.5,pointsize = 12) #define plot name size and font size
# 
# ggplot(mtcars)+
#   geom_point(aes(cyl, mpg))+
#   labs(x = expression(K==20,N1==1000,N2==200,Gamma==1))
# 
# dev.off()
# system("pdflatex mtcars.tex")
# 
# ggplot()+
#   labs(title = expression({k==1}*' and '*k==1*' asdf'))
# 
# ggplot()+
#   labs(title = expression({K==20}*', '*{N1==1000}*', '*{N2==200}*', '*{Gamma==1}))
