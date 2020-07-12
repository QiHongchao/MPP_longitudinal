##Simulation preparation
source("Functions_definition.R")

start_all <- Sys.time()
for (sc in 1:num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, num_simulation)
  
  ##Treatment effect
  trt_eff <- simulation_scenarios[sc, "trt"]
  
  ##Heterogeneity level
  heterogeneity <- simulation_scenarios[sc, "heterogeneity"]
  
  ##Empty vectors to store the result
  power_marg <- bias_marg <- se_marg <- mse_marg <- pow_par_marg <- vector()

for (ns in 1:num_simulation) {
  ##Simulated data generation
  set.seed(seeds_simulation[ns])
  
  ##Current data
  current <- datagen_current(num_tp = num_tp, num_group = num_group, num_fe = num_fe_current, num_re = num_re, 
                             beta = beta, var_b = var_b, var_e = var_e,
                             num_control = nsubpa, num_treat = nsubpa)
  ##Historical data
  historical <- datagen_historical(num_tp = num_tp, num_fe = num_fe_hist, num_re = num_re, 
                                   beta = beta, var_b = var_b, var_e = var_e, 
                                   num_control = nsubpa)
  
  ##Modify the data for stan
  data_marg <- sourceToList("./Marginal MPP/Marginal MPP_stan_data.R")
  
  ##Step 1: The calculation of scaling constant
  ##The order of variables and their index numbers are specified according to 
  ##"parameters block" in the Stan file  (beta, Omega, tau, b, sigmasq).
  beta_max <- data_marg$K0
  Omega_min <- beta_max + 1
  Omega_max <- Omega_min + (data_marg$L)^2 -1
  tau_min <- Omega_max + 1
  tau_max <- tau_min + data_marg$L - 1
  sigmasq_index <- tau_max + 1  
  
  ##Area under the curve and log scaling constant
  logsc <- matrix(NA, num_power, 4)
  logsc[, 1] <- power_all
  logsc[1, 2:3] <- 0
  
  ##The initial values of step 1
  init_s1 <- 0
  
  ##Sampling in step 1
  for (i in 2: num_power) {
    data_marg$power <- power_all[i]
    sample_fixedpow <- sampling(sm1_marg, data = data_marg, chains=1, iter = num_iter_s1, init = init_s1, refresh = 0)
    ##Last sample for a specific power parameter
    sample_power_last <- matrix(unlist(sample_fixedpow@sim$samples[[1]]), nrow=num_iter_s1)[num_iter_s1,]
    init_s1 <- list(list(beta = sample_power_last[1: beta_max],
                         Omega = matrix(sample_power_last[Omega_min: Omega_max], nrow = data_marg$L),
                         tau = sample_power_last[tau_min: tau_max],
                         sigmasq = sample_power_last[sigmasq_index]))
    ##AUC under each interval, AUC for 0 is 0, unweighted
    logsc[i, 2] <- mean(sample_fixedpow@sim$samples[[1]]$loglik_fixedpow[(num_iter_s1*perc_burnin_s1+1):num_iter_s1])
  }
  
  ##Weighted AUC
  logsc[2:num_power, 3] <- logsc[2:num_power, 2]*interval_width
  
  ##Cumulative area under the curve
  logsc[, 4] <- cumsum(logsc[, 3])
  
  ##Only keep power values and corresponding logsc
  logsc <- logsc[, c(1, 4)]
  
  ##Step 2: The sampling of parameters in posterior distribution
  ##Now we have a grid of log scaling constant corresponding to specific power parameter
  
  ##Sampling in step 2
  sample_marg <- sampling(sm2_marg, data = data_marg, chains = num_chains_s2, seed = seeds_simulation[ns], 
                                       refresh = 0, init = "random", iter = num_iter_s2, warmup=num_iter_s2*perc_burnin_s2)
  
  ##Extract the results from the samples
  ##Transform stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_marg)$beta_trt
  
  ##The type I error or statistical power for the treatment effect
  power_marg[ns] <- 
    (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0) 
  
  ##Bias and MSE
  ##No interaction effect
  if (simulation_scenarios[sc, "trt"] == 0) {
    bias_marg[ns] <- mean(beta_trt)
    mse_marg[ns] <- mean(beta_trt^2)
  }
  ##With interaction effect
  if (simulation_scenarios[sc, "trt"] == 1) {
    bias_marg[ns] <- mean(beta_trt - interaction_time_trt)
    mse_marg[ns] <- mean((beta_trt - interaction_time_trt)^2)
  }
  
  ##Standard error for the treatment effect
  se_marg[ns] <- sd(beta_trt)
  
  ##The results of power parameter
  pow_par_marg[ns] <- mean(as.data.frame(sample_marg)$power)
  
  ##Progress 
  print(paste0("scenario: ", sc, ", iteration: ", ns))
}
  write.csv(cbind(bias_marg, se_marg, mse_marg, pow_par_marg, power_marg),
            paste0("res_marg_", sc, ".csv"), row.names = F)
}
end_all <- Sys.time()
end_all - start_all 

##Summary of the results of different scenarios
res_marg <- matrix(NA, num_scenarios, 5)
colnames(res_marg) <- c("power", "bias", "se", "mse", "pow_par_IQR")

for (i in 1:nrow(res_marg)) {
  aa <- read.csv(paste0("res_marg_", i, ".csv"))
  res_marg[i, "power"] <- round(1 - mean(aa$power_marg), 3)
  res_marg[i, "bias"] <- paste0(round(mean(aa$bias_marg), 3), " (", round(quantile(aa$bias_marg, 0.025), 3), ", ",  round(quantile(aa$bias_marg, 0.975), 3),")") 
  res_marg[i, "se"] <- paste0(round(mean(aa$se_marg), 3), " (", round(quantile(aa$se_marg, 0.025), 3), ", ",  round(quantile(aa$se_marg, 0.975), 3),")") 
  res_marg[i, "mse"] <- paste0(round(mean(aa$mse_marg), 3), " (", round(quantile(aa$mse_marg, 0.025), 3), ", ",  round(quantile(aa$mse_marg, 0.975), 3),")")
  ##median and IQR
  res_marg[i, "pow_par_IQR"] <- paste0(round(median(aa$pow_par_marg), 2), 
                                       " (", round(quantile(aa$pow_par_marg, 0.25), 2), ", ",  round(quantile(aa$pow_par_marg, 0.75), 2),")")
}

write.csv(res_marg, "res_marg_summary.csv", row.names = F)
