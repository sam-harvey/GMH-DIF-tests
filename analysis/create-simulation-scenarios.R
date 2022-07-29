K = c(20, 100)
m = c(50)
m0 = c(10)
gamma = c(1, 10)
FH = c(0, 1, 5
       ,10, 20
)
N_ref = c(1000, 500)
OR = 1.5
sim = 10
mu_delta = c(0,1)
# N_foc=c(200, 500)

df_dist_variables = crossing(K, m, m0, mu_delta, gamma, FH, N_ref, OR, sim) %>% 
  mutate(N_foc = case_when(N_ref == 1000 ~ 200,
                           N_ref == 500 ~ 500)) %>% 
  mutate(simulation_file = glue("data/simulations/Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.RData") %>% 
           as.character())

create_simulation_scenarios = function(K = c(20, 100),
                                       m = c(50),
                                       m0 = c(10),
                                       gamma = c(1, 10),
                                       FH = c(0, 1, 5, 10, 20),
                                       N_ref = c(1000, 500),
                                       OR = 1.5,
                                       sim = 10,
                                       mu_delta = c(0,1)){
  crossing(K, m, m0, mu_delta, gamma, FH, N_ref, OR, sim) %>% 
    mutate(N_foc = case_when(N_ref == 1000 ~ 200,
                             N_ref == 500 ~ 500)) %>% 
    mutate(simulation_file = glue("data/simulations/Sim_GenDif_m{m}_K{K}_Gamma{gamma}_OR{OR}_FH{FH}_Nref{N_ref}_Nfoc{N_foc}_mudelta{mu_delta}_sims{sim}.RData") %>% 
             as.character())
}
