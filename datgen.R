##Step 0: Dataset preparation
set.seed(seeds_simulation[ns])

##Current data
current <- current_data_gen(heterogeneity = simulation_scenarios$heterogeneity[sc], fe = c(2, 1), 
                            trt = simulation_scenarios$trt_eff[sc], ntp=num_measures, ngroup = 2, 
                            varb = c(0.25, 0.25), varw=1, num_control = num_subject, num_treat = num_subject)
##Historical data
historical <- historical_data_gen(heterogeneity = simulation_scenarios$heterogeneity[sc], fe=c(2,1), 
                                  ntp=num_measures, varb = c(0.25, 0.25), varw=1, num_control = num_subject)