rm(list = ls())
source('libraries.R')

#Set publication font
loadfonts(device = "win")
par(family = "LM Roman 10")

df_tbl_1 = read_csv('data/table-1-data.csv', skip=1)
df_tbl_2 = read_csv('data/table-2-data.csv', skip=1)



#' Convert table format as in 
#'
#' @param df_in 
#'
#' @return
prepare_plot_df = function(df_in){
  df_in  %>% 
    rename(gamma = 4) %>% 
    pivot_longer(cols = Asym:Boot_3) %>% 
    mutate(name_position = str_extract(name, '[0-9]+')) %>% 
    #Convert to plotmath eqns to be interpreted in graphs
    mutate(facet_name = case_when(name_position %>% is.na() ~ 'Power~italic(W)',
                                  name_position == 1 ~ 'Power~italic(W)[ind]',
                                  name_position == 2 ~ 'Coverage~hat(gamma)[j]',
                                  name_position == 3 ~ 'Coverage~hat(gamma)[j]~-~hat(gamma)[l]'),
           name_type = str_extract(name, '^[^_]+')) %>% 
    mutate(plot_label = glue("{{K=={K}}}*', '*{{N[1]=={N1}}}*', '*{{N[2]=={N2}}}*', '*{{Gamma=={gamma}}}") %>% as.character(),
           plot_axis = factor(FH, levels = unique(.$FH)))
}

get_plot_series_labels = function(df_in){
  # plot_breaks = df_in$plot_label %>% unique %>% sort
  # plot_labels_parsed = map(plot_breaks, as.expression) %>% reduce(c)
  
  plot_breaks = df_in$plot_label %>% unique %>% sort
  
  plot_labs = c()
  
  for(i in plot_breaks){
    plot_labs = c(plot_labs, parse(text=i))
  }
  
  return(list(plot_breaks, plot_labs))
}

plot_function = function(df_in,
                         plot_breaks, 
                         plot_labels_parsed, 
                         legend_position = 'bottom',
                         title_value = expression('Power and coverage of test statistics '*italic(W)*' and '*italic(W)[ind]),
                         subtitle_value = expression('Under '*{mu[1]==mu[2]}==0),
                         x_value = 'FH',
                         color_aes,
                         linetype_aes,
                         shape_aes,
                         # fill_aes,
                         ylims){
  
  legend_guide = guide_legend(nrow = 4)
  
  df_in %>% 
    ggplot(aes(x = plot_axis 
               ,y = value
               ,colour = !!sym(color_aes)
               # ,fill = !!sym(fill_aes)
               ,linetype = !!sym(linetype_aes)
               ,shape = !!sym(shape_aes)
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
          legend.text=element_text(size=12))+
    scale_color_manual(values=brewer_pal(type = 'qual', palette = 2)(8), breaks = plot_breaks, labels = plot_labels_parsed)+
    scale_linetype_manual(values=seq(1,8), breaks = plot_breaks, labels = plot_labels_parsed)+
    scale_shape_manual(values=seq(1,8), breaks = plot_breaks, labels = plot_labels_parsed)+
    guides(
      colour = legend_guide,
      linetype = legend_guide,
      shape = legend_guide
    )+
    labs(
      title = title_value
      ,subtitle=subtitle_value
      ,y = NULL
      ,x = x_value
      ,colour = NULL
      ,linetype = NULL
      ,shape = NULL
    )+
    scale_y_continuous(limits = ylims)
}

plot_main = function(df_in, 
                     subtitle_value,
                     color_aes = 'plot_label',
                     linetype_aes = 'plot_label',
                     shape_aes = 'plot_label',
                     fill_aes = 'plot_label',
                     ylims){
  plot_df = prepare_plot_df(df_in)
  plot_series_labels = get_plot_series_labels(plot_df)
  
  df_power_stats = plot_df %>% 
    filter(str_detect(facet_name, 'Power'))
  
  p1 = plot_function(df_power_stats, 
                     plot_series_labels[[1]], 
                     plot_series_labels[[2]], 
                     legend_position = 'none',
                     subtitle_value = subtitle_value,
                     color_aes = color_aes,
                     linetype_aes = linetype_aes,
                     shape_aes = shape_aes,
                     # fill_aes = fill_aes,
                     ylims = ylims[[1]],
                     x_value = NULL)
  
  df_coverage_stats = plot_df %>% 
    filter(str_detect(facet_name, 'Coverage'))
  
  p2 = plot_function(df_coverage_stats, 
                     plot_series_labels[[1]], 
                     plot_series_labels[[2]],
                     title_value = NULL,
                     subtitle_value = NULL,
                     color_aes = color_aes,
                     linetype_aes = linetype_aes,
                     shape_aes = shape_aes,
                     # fill_aes = fill_aes,
                     ylims = ylims[[2]])
  
  p3 = p1/p2
  
  return(p3)
}

# df_in = df_tbl_1
# 
# plot_df = prepare_plot_df(df_in)
# 
# plot_breaks = plot_df$plot_label %>% unique %>% sort
# plot_labels_parsed = map(plot_breaks, as.expression) %>% reduce(c)
# 
# 
# df_power_stats = plot_df %>% 
#   filter(str_detect(facet_name, 'Power'))
# 
# p1 = plot_function(df_power_stats, plot_series_labels[[1]], plot_series_labels[[2]], legend_position = 'none', x_value = NULL)
# 
# df_coverage_stats = plot_df %>% 
#   filter(str_detect(facet_name, 'Coverage'))
# 
# p2 = plot_function(df_coverage_stats, plot_series_labels[[1]], plot_series_labels[[2]], title_value = NULL)
# 
# p3 = p1/p2
# 
# return(p3)

# debugonce(plot_function)
# debugonce(plot_main)

plot_height = 12
plot_width = 10

df_tbl_1_plot = plot_main(df_tbl_1, 
                          subtitle_value = expression('Under '*{mu[1]==mu[2]}==0),
                          ylims = list(c(0,1), c(0.93, 0.97)))

df_tbl_2_plot = plot_main(df_tbl_2, 
                          subtitle_value = expression('Under '*{mu[1]==0}*' and '*{mu[2]==1}),
                          ylims = list(c(0,1), c(0.93, 0.97)))

ggsave('output/plots/table-1-plot-all-aes-v2.png',
       df_tbl_1_plot,
       h=plot_height,
       w=plot_width)

ggsave('output/plots/table-2-plot-all-aes-v2.png',
       df_tbl_2_plot,
       h=plot_height,
       w=plot_width)

df_tbl_1_plot = plot_main(df_tbl_1, 
                          subtitle_value = expression('Under '*{mu[1]==mu[2]}==0),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '')

df_tbl_2_plot = plot_main(df_tbl_2, 
                          subtitle_value = expression('Under '*{mu[1]==0}*' and '*{mu[2]==1}),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '')


ggsave('output/plots/table-1-plot-shape-linetype-aes.png',
       df_tbl_1_plot,
       h=plot_height,
       w=plot_width)

ggsave('output/plots/table-2-plot-shape-linetype-aes.png',
       df_tbl_2_plot,
       h=plot_height,
       w=plot_width)

df_tbl_1_plot = plot_main(df_tbl_1, 
                          subtitle_value = expression('Under '*{mu[1]==mu[2]}==0),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '',
                          linetype_aes = '')

df_tbl_2_plot = plot_main(df_tbl_2, 
                          subtitle_value = expression('Under '*{mu[1]==0}*' and '*{mu[2]==1}),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '',
                          linetype_aes = '')

ggsave('output/plots/table-1-plot-shape-aes.png',
       df_tbl_1_plot,
       h=plot_height,
       w=plot_width)

ggsave('output/plots/table-2-plot-shape-aes.png',
       df_tbl_2_plot,
       h=plot_height,
       w=plot_width)

df_tbl_1_plot = plot_main(df_tbl_1, 
                          subtitle_value = expression('Under '*{mu[1]==mu[2]}==0),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '',
                          shape_aes = '')

df_tbl_2_plot = plot_main(df_tbl_2, 
                          subtitle_value = expression('Under '*{mu[1]==0}*' and '*{mu[2]==1}),
                          ylims = list(c(0,1), c(0.94, 0.97)),
                          color_aes = '',
                          shape_aes = '')

ggsave('output/plots/table-1-plot-linetype-aes.png',
       df_tbl_1_plot,
       h=plot_height,
       w=plot_width)

ggsave('output/plots/table-2-plot-linetype-aes.png',
       df_tbl_2_plot,
       h=plot_height,
       w=plot_width)
