##Simulation preparation
source("Functions_definition.R")

start <- Sys.time()
for (sc in 1:num_scenarios) {
  ##random seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, num_simulation)
  ##Step 0: generate simulated data and settings for the following two steps
  ##treatment effect
  trt_eff <- trt_effect(trt = simulation_scenarios[sc, "trt"])
  
  ##heterogeneity level
  heterogeneity <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, "heterogeneity"]) 
  
  power_beta_trt <- bias_beta_trt <- sd_beta_trt <- mse_beta_trt <- vector()
  
  for (ns in 1:num_simulation) {
    
    ##Generate simulated data
    source("datgen.R")
    
    data_current <- sourceToList("./No borrowing/Current_stan_data.R")
    
    ##sampling
    sample_current <- sampling(sm_current, data = data_current, chains = num_chains, 
                                seed = seeds_simulation[ns], init = "random",
                                iter = num_iter, warmup=num_iter*perc_burnin, refresh = 0, cores = detectCores(),
                                control = list(stepsize = 0.99, max_treedepth = 15))

  ##Extract the results from the samples
  ##Transform stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_current)$"beta[3]"
  
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
  
  ##Progress
  print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
  write.csv(cbind(bias_beta_trt, sd_beta_trt, mse_beta_trt, power_beta_trt),
            paste0("res_current_", sc, ".csv"), row.names = F)
}
end <- Sys.time()
end - start

res_current <- matrix(NA, num_scenarios, 4)
colnames(res_current) <- c("power", "bias", "se", "mse")
for (i in 1:nrow(res_current)) {
  aa <- read.csv(paste0("res_current_", i, ".csv"))
  res_current[i, "power"] <- round(1 - mean(aa$power_beta_trt), 3)
  res_current[i, "bias"] <- paste0(round(mean(aa$bias_beta_trt), 3), " (", round(quantile(aa$bias_beta_trt, 0.025), 3), ", ",  round(quantile(aa$bias_beta_trt, 0.975), 3),")") 
  res_current[i, "se"] <- paste0(round(mean(aa$sd_beta_trt), 3), " (", round(quantile(aa$sd_beta_trt, 0.025), 3), ", ",  round(quantile(aa$sd_beta_trt, 0.975), 3),")") 
  res_current[i, "mse"] <- paste0(round(mean(aa$mse_beta_trt), 3), " (", round(quantile(aa$mse_beta_trt, 0.025), 3), ", ",  round(quantile(aa$mse_beta_trt, 0.975), 3),")")
}
write.csv(res_current, "res_current_summary.csv", row.names = F)