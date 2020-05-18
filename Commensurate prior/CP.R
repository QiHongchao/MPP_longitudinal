##Simulation preparation
source("Functions_definition.R")

start <- Sys.time()
for (sc in 1:num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, num_simulation)
  
  ##treatment effect
  trt_eff <- trt_effect(trt = simulation_scenarios[sc, "trt"])
  
  ##heterogeneity level
  heterogeneity <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, "heterogeneity"])  
  
  power_beta_trt <- bias_beta_trt <- sd_beta_trt <- mse_beta_trt <- sigmabeta0 <- sigmabeta1 <- vector()
  
  
  for (ns in 1:num_simulation) {
    
    ##Generate simulated data
    source("datgen_HC.R")
    
    data_cp <- sourceToList("./Commensurate prior/CP_stan_data.R")
    
    ##sampling
    sample_cp <- sampling(sm_cp, data = data_cp, chains = num_chains, 
                                seed = seeds_simulation[ns], init = "random",
                                iter = num_iter, warmup=num_iter*perc_burnin, refresh = 0, cores = detectCores(),
                                control = list(stepsize = 0.99, max_treedepth = 15))

  ##Extract the results from the samples
  ##Transform stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_cp)$beta_trt
  
  ##The type I error or statistical power for interaction
  power_beta_trt[ns] <- 
      (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0) 
  
  ##Bias and MSE
  ##No interaction effect
  if (simulation_scenarios[sc, 1] == 0) {
    bias_beta_trt[ns] <- mean(beta_trt)
    mse_beta_trt[ns] <- mean(beta_trt^2)
    }
  ##With interaction effect
  if (simulation_scenarios[sc, 1] == 1) {
    bias_beta_trt[ns] <- mean(beta_trt - interaction_time_trt)
    mse_beta_trt[ns] <- mean((beta_trt - interaction_time_trt)^2)
    }
  ##Standard deviation
  sd_beta_trt[ns] <- sd(beta_trt)
  
  ##The between-study heterogeneity
  sigmabeta0[ns] <- mean(as.data.frame(sample_cp)$"tau_beta[1]")
  sigmabeta1[ns] <- mean(as.data.frame(sample_cp)$"tau_beta[2]")
  
  ##Progress
  print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
  write.csv(cbind(bias_beta_trt, sd_beta_trt, mse_beta_trt, power_beta_trt, sigmabeta0, sigmabeta1),
            paste0("res_cp_", sc, ".csv"), row.names = F)
}
end <- Sys.time()
end - start

##save the results
res_cp <- matrix(NA, num_scenarios, 6)
colnames(res_cp) <- c("power", "bias", "se", "mse", "sigmabeta0", "sigmabeta1")
for (i in 1:nrow(res_cp)) {
  aa <- read.csv(paste0("res_cp_", i, ".csv"))
  res_cp[i, "power"] <- round(1 - mean(aa$power_beta_trt), 3)
  res_cp[i, "bias"] <- paste0(round(mean(aa$bias_beta_trt), 3), " (", round(quantile(aa$bias_beta_trt, 0.025), 3), ", ",  round(quantile(aa$bias_beta_trt, 0.975), 3),")") 
  res_cp[i, "se"] <- paste0(round(mean(aa$sd_beta_trt), 3), " (", round(quantile(aa$sd_beta_trt, 0.025), 3), ", ",  round(quantile(aa$sd_beta_trt, 0.975), 3),")") 
  res_cp[i, "mse"] <- paste0(round(mean(aa$mse_beta_trt), 3), " (", round(quantile(aa$mse_beta_trt, 0.025), 3), ", ",  round(quantile(aa$mse_beta_trt, 0.975), 3),")")
  ##median and IQR for sigmabeta0 and sigmabeta1
  res_cp[i, "sigmabeta0"] <- paste0(round(mean(aa$sigmabeta0), 2), 
                                    " (", round(quantile(aa$sigmabeta0, 0.25), 2), ", ",  round(quantile(aa$sigmabeta0, 0.75), 2),")")
  res_cp[i, "sigmabeta1"] <- paste0(round(mean(aa$sigmabeta1), 2), 
                                    " (", round(quantile(aa$sigmabeta1, 0.25), 2), ", ",  round(quantile(aa$sigmabeta1, 0.75), 2),")")
}

write.csv(res_cp, "res_cp_summary.csv", row.names = F)
