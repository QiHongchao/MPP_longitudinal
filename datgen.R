##Step 0: Dataset preparation
set.seed(seeds_simulation[ns])

##Current data
current <- design_matrix_gen_current(ntp=simulation_scenarios[sc, 2], num_re=2, fe=c(2,1), varb1=0.25, varb2=0.25, varw=1,
                                     num_control = simulation_scenarios[sc, 3], num_treat = simulation_scenarios[sc, 3])
##Historical data
historical <- design_matrix_gen_hc(ntp=simulation_scenarios[sc, 2], num_re=2, fe=c(2,1), varb1=0.25, varb2=0.25, varw=1,
                                   num_control = simulation_scenarios[sc, 3])