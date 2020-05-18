##Simulation preparation
source("Functions_definition.R")

start_all <- Sys.time()
for (sc in 1:num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, num_simulation)
  
  ##Step 0: generate simulated data and settings for the following two steps
  ##treatment effect
  trt_eff <- trt_effect(trt = simulation_scenarios[sc, "trt"])
  
  ##heterogeneity level
  heterogeneity <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, "heterogeneity"])
  
  power_cond <- bias_cond <- sd_cond <- mse_cond <- pow_par_cond <- vector()
  
  for (ns in 1:num_simulation) {
    
    ##Generate simulated data
    source("datgen.R")
    
    ##Step 1: sampling/calculation of average log likelihood
    
    ##Generate data for the Stan file of step 1 
    data_cond <- sourceToList("./Conditional MPP/Conditional MPP_stan_data.R")
    
    ##The order of variables and their index numbers are specified according to 
    ##"parameters block" in the Stan file  (beta, Omega, tau, b, sigmasq).
    beta_max <- data_cond$K0
    Omega_min <- beta_max + 1
    Omega_max <- Omega_min + (data_cond$L)^2 -1
    tau_min <- Omega_max + 1
    tau_max <- tau_min + data_cond$L - 1
    b_min <- tau_max + 1
    b_max <- b_min + data_cond$J0 * data_cond$L -1
    sigmasq_index <- b_max + 1  
    
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
                           b = matrix(sample_power_last[b_min: b_max], data_cond$J0, data_cond$L), 
                           sigmasq = sample_power_last[sigmasq_index]))
      ##AUC under each interval, AUC for 0 is 0, unweighted
      logsc[i, 2] <- mean(sample_fixedpow@sim$samples[[1]]$loglik_fixedpow[(num_iter_s1[i]*perc_burnin_s1+1):num_iter_s1[i]])
    }
    ##Weighted AUC
    logsc[2:num_power, 3] <- logsc[2:num_power, 2]*interval_width
    ##Cumulative area under the curve
    logsc[, 4] <- cumsum(logsc[, 3])
    
    ##Step 2: Sampling for all parameters
    ##warmup phase for step 2
    sample_cond_fixedpow <- sampling(sm2_cond_fixedpow, data = data_cond, chains = 1,
                                       init = "random", seed = 1234,
                                       iter = 100, warmup=50, refresh = 0, cores = detectCores(),
                                       control = list(stepsize = 0.01, adapt_delta = 0.95))
      
    sample_power_last <- as.matrix(as.data.frame(sample_cond_fixedpow))[50,]
      
    init_s2 <- list(beta_mutual = sample_power_last[1:data_cond$K0],
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
      ##This is one advantage of HMC, we don't need a lot of iterations any more.
      ##In step 1, the log scaling constant becomes stable after 200 iterations.
      sample_cond <- sampling(sm2_cond, data = data_cond, chains = num_chains_s2,
                                           init = list(init_s2, init_s2, init_s2, init_s2),
                                           iter = num_iter_s2, warmup=num_iter_s2*perc_burnin_s2, refresh = 0, cores = detectCores(),
                                           control = list(stepsize = 0.01, adapt_delta = 0.95))
    
    ##Extract the results from the samples
    ##Transform stanfit to data frames and summarize parameters of interest
    beta_trt <- as.data.frame(sample_cond)$beta_trt

    ##The type I error or statistical power for interaction
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
    ##Standard deviation
    sd_cond[ns] <- sd(beta_trt)

    ##The results of power parameter
    pow_par_cond[ns] <- mean(as.data.frame(sample_cond)$power)

    ##Progress
    print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
    write.csv(cbind(bias_cond, sd_cond, mse_cond, pow_par_cond, power_cond),
            paste0("res_cond_", sc, ".csv"), row.names = F)
}
end_all <- Sys.time()
end_all - start_all

##save the complete results
res_cond <- matrix(NA, num_scenarios, 5)
colnames(res_cond) <- c("power", "bias", "se", "mse", "pow_par_IQR")
for (i in 1:nrow(res_cond)) {
  aa <- read.csv(paste0("res_cond_", i, ".csv"))
  res_cond[i, "power"] <- round(1 - mean(aa$power_cond), 3)
  res_cond[i, "bias"] <- paste0(round(mean(aa$bias_cond), 3), " (", round(quantile(aa$bias_cond, 0.025), 3), ", ",  round(quantile(aa$bias_cond, 0.975), 3),")") 
  res_cond[i, "se"] <- paste0(round(mean(aa$sd_cond), 3), " (", round(quantile(aa$sd_cond, 0.025), 3), ", ",  round(quantile(aa$sd_cond, 0.975), 3),")") 
  res_cond[i, "mse"] <- paste0(round(mean(aa$mse_cond), 3), " (", round(quantile(aa$mse_cond, 0.025), 3), ", ",  round(quantile(aa$mse_cond, 0.975), 3),")")
  res_cond[i, "pow_par_IQR"] <- paste0(round(median(aa$pow_par_cond), 2), 
                                       " (", round(quantile(aa$pow_par_cond, 0.25), 2), ", ",  round(quantile(aa$pow_par_cond, 0.75), 2),")")
}

write.csv(res_cond, "res_cond_summary.csv", row.names = F)
