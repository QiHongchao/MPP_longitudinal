##Step 0: Dataset preparation
set.seed(seeds_simulation[ns])

##Current data
current <- design_matrix_gen_current(ntp = ntp, group = 2, num_fe = 3, num_re = 2, beta = c(2,1), var_b = rep(0.25, 2), var_e = 1,
                                     num_control = nsubpa, num_treat = nsubpa)
##Historical data
historical <- design_matrix_gen_hc(ntp=ntp, num_fe = 2, num_re=2, beta = c(2,1), var_b = rep(0.25, 2), var_e = 1, 
                                   num_control = nsubpa)