##Simulation preparation
source("../Functions_definition.R")

##Stan file compilation
sm_cp <- stan_model("CP_Stan.stan")

##Samples for different scenarios
sample_cp_scenarios <- list(NULL)

start <- Sys.time()
for (sc in 1:Num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, Num_simulation)
  heterogeneity_trt <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, 4], trt = simulation_scenarios[sc, 1])  
  power_beta_trt <- mean_beta_trt <- sd_beta_trt <- sigma_beta0 <- sigma_beta1 <- vector()
  
  for (ns in 1:Num_simulation) {
    
    ##Generate simulated data
    source("../datgen.R")
    
    data_cp <- sourceToList("CP_Stan_data.R")
    
    init <- list(beta_01_hist = c(2, 1), beta_trt = as.numeric(sc>7)*interaction_time_trt,
                 tau0 = rep(0.5, 2), Omega0 = diag(1, 2), tau1 = rep(0.5, 2), Omega1 = diag(1, 2),
                 sigma0 = 1, sigma1 = 1)
    
    ##sampling
    sample_cp <- sampling(sm_cp, data = data_cp, chains = Num_chains, 
                                init = list(init, init, init, init),
                                iter = Num_iter, warmup=Num_iter*perc_burnin, refresh = 0, cores = Num_chains)

  ##Extract the results from the samples
  ##Transform Stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_cp)$beta_trt
  
  ##The type I error or statistical power for the treatment effect
  power_beta_trt[ns] <- 
    (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0) 
  
  ##Posterior mean of the treatment effect
  mean_beta_trt[ns] <- mean(beta_trt)
  ##Posterior standard deviation of the treatment effect
  sd_beta_trt[ns] <- sd(beta_trt)
  
  ##The between-study heterogeneity
  sigma_beta0[ns] <- mean(as.data.frame(sample_cp)$"tau_beta[1]")
  sigma_beta1[ns] <- mean(as.data.frame(sample_cp)$"tau_beta[2]")
  
  ##Progress
  print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
  sample_cp_scenarios[[sc]] <- cbind(power_beta_trt, mean_beta_trt, sd_beta_trt, sigma_beta0, sigma_beta1)
}
end <- Sys.time()
end - start

##Summary of the results
summary_cp <- data.frame(power = rep(NA, Num_scenarios), bias = NA, se = NA, mse = NA,
                         sigma_beta0 = NA, sigma_beta1 = NA)
for (sc in 1:Num_scenarios) {
  res <- sample_cp_scenarios[[sc]]
  ##Type I error rate/Power
  summary_cp[sc, "power"] <- round((1 - mean(res$power_beta_trt)) * 100, 1)
  ##Bias
  bias <- res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt
  summary_cp[sc, "bias"] <- paste0(round(mean(bias), 3), " (", round(quantile(bias, 0.025), 3), ", ",  round(quantile(bias, 0.975), 3),")") 
  ##SE
  summary_cp[sc, "se"] <- paste0(round(mean(res$sd_beta_trt), 3), " (", round(quantile(res$sd_beta_trt, 0.025), 3), ", ",  round(quantile(res$sd_beta_trt, 0.975), 3),")") 
  ##MSE
  mse <- (res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt)^2
  summary_cp[sc, "mse"] <- paste0(round(mean(mse), 3), " (", round(quantile(mse, 0.025), 3), ", ",  round(quantile(mse, 0.975), 3),")")
  ##Median and IQR for sigma_beta0 and sigma_beta1
  summary_cp[sc, "sigma_beta0"] <- paste0(round(mean(res$sigma_beta0), 2), 
                                    " (", round(quantile(res$sigma_beta0, 0.25), 2), ", ",  round(quantile(res$sigma_beta0, 0.75), 2),")")
  summary_cp[sc, "sigma_beta1"] <- paste0(round(mean(res$sigma_beta1), 2), 
                                    " (", round(quantile(res$sigma_beta1, 0.25), 2), ", ",  round(quantile(res$sigma_beta1, 0.75), 2),")")
}
