##Simulation preparation
source("../Functions_definition.R")

##The candidate power parameters in step 1
interval_width <- 0.02
power_all <- seq(0, 1, interval_width)
num_power <- length(power_all)

##Sampling setting for step 1
Num_iter_s1 <- 100

##Data frame to save log scaling constant
logsc_all <- data.frame(logsc = rep(NA, num_power * Num_simulation))

##Compile Stan files for step 1 and step 2
sm1_mmpp <- stan_model("MMPP_Stan_S1.stan")
sm2_mmpp <- stan_model("MMPP_Stan_S2.stan")

##Empty matrix for the results of the marginal MPP
sample_mmpp_scenarios <- list(NULL)

start_all <- Sys.time()
for (sc in 1:Num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, Num_simulation)
  heterogeneity_trt <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, 4], trt = simulation_scenarios[sc, 1])
  power_beta_trt <- mean_beta_trt <- sd_beta_trt <- pow_par <- vector()

for (ns in 1:Num_simulation) {
  ##Generate simulated data sets
  source("../datgen.R")
  
  ##Step 1: sampling/calculation of average log likelihood
  
  ##Generate data for the Stan file of step 1 
  data_s1 <- sourceToList("MMPP_Stan_data_S1.R")
  
  ##Step 1: The calculation of scaling constant
  ##The order of variables and their index numbers are specified according to 
  ##"parameters block" in the Stan file  (beta, Omega, tau, b, sigma).
  beta_max <- data_s1$K
  Omega_min <- beta_max + 1
  Omega_max <- Omega_min + (data_s1$L)^2 -1
  tau_min <- Omega_max + 1
  tau_max <- tau_min + data_s1$L - 1
  sigma_index <- tau_max + 1  
  
  ##Area under the curve and log scaling constant
  logsc <- matrix(NA, num_power, 4)
  logsc[, 1] <- power_all
  logsc[1, 2:3] <- 0
  
  ##The initial values of step 1
  init_s1 <- list(list(beta = c(2, 1), Omega = diag(1, 2),
                       tau = rep(0.5, 2), sigma = 1))
  
  ##Sampling in step 1
  for (i in 2: num_power) {
    data_s1$power <- power_all[i]
    sample_fixedpow <- sampling(sm1_mmpp, data = data_s1, chains=1, iter = Num_iter_s1, init = init_s1, refresh = 0)
    ##Last sample for a specific power parameter
    sample_power_last <- as.data.frame(sample_fixedpow)[nrow(as.data.frame(sample_fixedpow)),]
    init_s1 <- list(list(beta = sample_power_last[1: beta_max],
                         Omega = matrix(sample_power_last[Omega_min: Omega_max], nrow = data_s1$L),
                         tau = sample_power_last[tau_min: tau_max],
                         sigma = sample_power_last[sigma_index]))
    ##AUC under each interval, AUC for 0 is 0, unweighted
    logsc[i, 2] <- mean(as.data.frame(sample_fixedpow)$loglik_fixedpow)
  }
  ##Weighted AUC
  logsc[2:num_power, 3] <- logsc[2:num_power, 2]*interval_width
  ##Cumulative area under the curve
  logsc[, 4] <- cumsum(logsc[, 3])
  
  ##save logsc
  logsc_all$logsc[((ns - 1) * num_power + 1):(ns * num_power)] <- logsc[,4]
  
  ##Step 2: The sampling of parameters in posterior distribution
  ##Now we have a grid of log scaling constant corresponding to specific power parameter
  ##Data preparation for step 2
  data_s2 <- sourceToList("MMPP_Stan_data_S2.R")
  
  ##initial values for step 2
  init_s2 <- list(beta_mutual = c(2, 1), beta_trt = as.numeric(sc>7)*interaction_time_trt,
                  Omega = diag(1, 2), tau = rep(0.5, 2),
                  sigma = 1)
  
  ##Sampling in step 2
  sample_mmpp <- sampling(sm2_mmpp, data = data_s2, chains = Num_chains, 
                                       init = list(init_s2, init_s2, init_s2, init_s2),
                                       iter = Num_iter, warmup=Num_iter*perc_burnin, refresh = 0, 
                                       cores = Num_chains)
  
  ##Extract the results from the samples
  ##Transform Stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_mmpp)$beta_trt
  
  ##The type I error or statistical power for the treatment effect
  power_beta_trt[ns] <- 
    (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0) 
  
  ##Posterior mean of the treatment effect
  mean_beta_trt[ns] <- mean(beta_trt)
  ##Posterior standard deviation of the treatment effect
  sd_beta_trt[ns] <- sd(beta_trt)
  
  ##Mean power parameter
  pow_par[ns] <- mean(as.data.frame(sample_mmpp)$power)
  
  ##Progress 
  print(paste0("scenario: ", sc, ", iteration: ", ns))
}
  sample_mmpp_scenarios[[sc]] <- cbind(power_beta_trt, mean_beta_trt, sd_beta_trt, pow_par)
}
end_all <- Sys.time()
end_all - start_all 

##Summary of the results
summary_mmpp <- data.frame(power = rep(NA, Num_scenarios), bias = NA, se = NA, mse = NA,
                           pow_par = NA)
for (sc in 1:Num_scenarios) {
  res <- sample_mmpp_scenarios[[sc]]
  ##Type I error rate/Power
  summary_mmpp[sc, "power"] <- round((1 - mean(res$power_beta_trt)) * 100, 1)
  ##Bias
  bias <- res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt
  summary_mmpp[sc, "bias"] <- paste0(round(mean(bias), 3), " (", round(quantile(bias, 0.025), 3), ", ",  round(quantile(bias, 0.975), 3),")") 
  ##SE
  summary_mmpp[sc, "se"] <- paste0(round(mean(res$sd_beta_trt), 3), " (", round(quantile(res$sd_beta_trt, 0.025), 3), ", ",  round(quantile(res$sd_beta_trt, 0.975), 3),")") 
  ##MSE
  mse <- (res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt)^2
  summary_mmpp[sc, "mse"] <- paste0(round(mean(mse), 3), " (", round(quantile(mse, 0.025), 3), ", ",  round(quantile(mse, 0.975), 3),")")
  ##Median and IQR for the power parameter
  summary_mmpp[sc, "pow_par"] <- paste0(round(mean(res$pow_par), 2), 
                                        " (", round(quantile(res$pow_par, 0.25), 2), ", ",  round(quantile(res$pow_par, 0.75), 2),")")
}