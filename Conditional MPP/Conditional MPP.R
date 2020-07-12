##Simulation preparation
source("Functions_definition.R")

start_all <- Sys.time()
for (sc in 1:num_scenarios) {
  ##Random seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, num_simulation)
  
  ##Treatment effect
  trt_eff <- simulation_scenarios[sc, "trt"]
  
  ##Heterogeneity level
  heterogeneity <- simulation_scenarios[sc, "heterogeneity"]
  
  ##Empty vectors to store the result
  power_cond <- bias_cond <- se_cond <- mse_cond <- pow_par_cond <- vector()
  
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
    
    ##Step 1: sampling/calculation of average log likelihood
    
    ##Generate data for stan
    data_cond <- sourceToList("./Conditional MPP/Conditional MPP_stan_data.R")
    
    ##The order of variables and their index numbers are specified according to 
    ##"parameters block" in the Stan file  (beta, Omega, tau, b, sigmasq).
    beta_max <- data_cond$K0
    Omega_min <- beta_max + 1
    Omega_max <- Omega_min + (data_cond$L)^2 -1
    tau_min <- Omega_max + 1
    tau_max <- tau_min + data_cond$L - 1
    b0_min <- tau_max + 1
    b0_max <- b0_min + data_cond$J0 * data_cond$L -1
    sigmasq_index <- b0_max + 1  
    
    ##Area under the curve and log scaling constant
    logsc <- matrix(NA, num_power, 4)
    logsc[, 1] <- power_all
    logsc[1, 2:3] <- 0
    
    ##The initial values of step 1
    init_s1 <- 0
    
    ##Sampling in step 1
    for (i in 2: num_power) {
      data_cond$power <- power_all[i]
      sample_fixedpow <- sampling(sm1_cond, data = data_cond, chains=1, iter = num_iter_s1[i], init = init_s1, refresh = 0)
      ##Last sample for a specific power parameter
      sample_power_last <- matrix(unlist(sample_fixedpow@sim$samples[[1]]), nrow=num_iter_s1[i])[num_iter_s1[i],]
      init_s1 <- list(list(beta = sample_power_last[1: beta_max],
                           Omega = matrix(sample_power_last[Omega_min: Omega_max], nrow = data_cond$L),
                           tau = sample_power_last[tau_min: tau_max],
                           b0 = matrix(sample_power_last[b0_min: b0_max], data_cond$J0, data_cond$L), 
                           sigmasq = sample_power_last[sigmasq_index]))
      ##AUC under each interval, AUC for 0 is 0, unweighted
      logsc[i, 2] <- mean(sample_fixedpow@sim$samples[[1]]$loglik_fixedpow[(num_iter_s1[i]*perc_burnin_s1+1):num_iter_s1[i]])
    }
    
    ##Weighted AUC
    logsc[2:num_power, 3] <- logsc[2:num_power, 2] * interval_width
    
    ##Cumulative area under the curve
    logsc[, 4] <- cumsum(logsc[, 3])
    
    ##Only keep power values and corresponding logsc
    logsc <- logsc[, c(1, 4)]
    
    ##Step 2: Sampling for all parameters
    
    ##Warmup phase for step 2
    sample_cond_fixedpow <- sampling(sm2_cond_fixedpow, data = data_cond, chains = 1,
                                       init = "random", seed = 1234,
                                       iter = 100, warmup=50, refresh = 0, cores = detectCores(),
                                       control = list(stepsize = 0.01, adapt_delta = 0.95))
      
    sample_power_last <- as.matrix(as.data.frame(sample_cond_fixedpow))[50,]
      
    init_s2 <- list(beta = sample_power_last[1:data_cond$K0],
                    beta_trt = sample_power_last[data_cond$K],
                    Omega = matrix(sample_power_last[(data_cond$K+1): (data_cond$K+data_cond$L^2)], data_cond$L),
                    tau = sample_power_last[(data_cond$K+data_cond$L^2 + 1):
                                                (data_cond$K+data_cond$L^2 + data_cond$L)],
                    b0 = matrix(sample_power_last[(data_cond$K+data_cond$L^2 + data_cond$L + 1):
                                                      (data_cond$K+data_cond$L^2 + data_cond$L + data_cond$J0*data_cond$L)], data_cond$J0),
                    b = matrix(sample_power_last[(data_cond$K+data_cond$L^2 + data_cond$L + data_cond$J0*data_cond$L + 1):
                                                     (data_cond$K+data_cond$L^2 + data_cond$L + data_cond$J0*data_cond$L + data_cond$J*data_cond$L)], data_cond$J),
                    sigmasq = sample_power_last[data_cond$K+data_cond$L^2 + data_cond$L + data_cond$J0*data_cond$L + data_cond$J*data_cond$L + 1],
                    power  = 0.5
      )
      
      ##Sampling in step 2
      sample_cond <- sampling(sm2_cond, data = data_cond, chains = num_chains_s2,
                                           init = list(init_s2, init_s2, init_s2, init_s2),
                                           iter = num_iter_s2, warmup=num_iter_s2*perc_burnin_s2, refresh = 0, cores = detectCores(),
                                           control = list(stepsize = 0.01, adapt_delta = 0.95))
    
    ##Extract the results from the samples
    ##Transform stanfit to data frames and summarize parameters of interest
    beta_trt <- as.data.frame(sample_cond)$beta_trt

    ##The type I error or statistical power for the treatment effect
    power_cond[ns] <-
      (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0)

    ##Bias and MSE
    ##No interaction effect
    if (simulation_scenarios[sc, "trt"] == 0) {
      bias_cond[ns] <- mean(beta_trt)
      mse_cond[ns] <- mean(beta_trt^2)
    }
    ##With interaction effect
    if (simulation_scenarios[sc, "trt"] == 1) {
      bias_cond[ns] <- mean(beta_trt - interaction_time_trt)
      mse_cond[ns] <- mean((beta_trt - interaction_time_trt)^2)
    }
    
    ##Standard error for the treatment effect
    se_cond[ns] <- sd(beta_trt)

    ##The results of power parameter
    pow_par_cond[ns] <- mean(as.data.frame(sample_cond)$power)

    ##Progress
    print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
    write.csv(cbind(bias_cond, se_cond, mse_cond, pow_par_cond, power_cond),
            paste0("res_cond_", sc, ".csv"), row.names = F)
}
end_all <- Sys.time()
end_all - start_all

##Summary of the results of different scenarios
res_cond <- matrix(NA, num_scenarios, 5)
colnames(res_cond) <- c("power", "bias", "se", "mse", "pow_par_IQR")
for (i in 1:nrow(res_cond)) {
  aa <- read.csv(paste0("res_cond_", i, ".csv"))
  res_cond[i, "power"] <- round(1 - mean(aa$power_cond), 3)
  res_cond[i, "bias"] <- paste0(round(mean(aa$bias_cond), 3), " (", round(quantile(aa$bias_cond, 0.025), 3), ", ",  round(quantile(aa$bias_cond, 0.975), 3),")") 
  res_cond[i, "se"] <- paste0(round(mean(aa$se_cond), 3), " (", round(quantile(aa$se_cond, 0.025), 3), ", ",  round(quantile(aa$se_cond, 0.975), 3),")") 
  res_cond[i, "mse"] <- paste0(round(mean(aa$mse_cond), 3), " (", round(quantile(aa$mse_cond, 0.025), 3), ", ",  round(quantile(aa$mse_cond, 0.975), 3),")")
  res_cond[i, "pow_par_IQR"] <- paste0(round(median(aa$pow_par_cond), 2), 
                                       " (", round(quantile(aa$pow_par_cond, 0.25), 2), ", ",  round(quantile(aa$pow_par_cond, 0.75), 2),")")
}
write.csv(res_cond, "res_cond_summary.csv", row.names = F)
