##Simulation preparation
source("../Functions_definition.R")

##Stan file compilation
sm_pooling <- stan_model("Pooling_Stan.stan")

##Empty matrix for the results of the marginal MPP
sample_pooling_scenarios <- list(NULL)

start <- Sys.time()
for (sc in 1:Num_scenarios) {
  ##seeds generation
  set.seed(seeds_scenarios[sc])
  seeds_simulation <- sample.int(.Machine$integer.max, Num_simulation)
  heterogeneity_trt <- diff_between_current_hist(heterogeneity = simulation_scenarios[sc, 4], trt = simulation_scenarios[sc, 1])  
  power_beta_trt <- mean_beta_trt <- sd_beta_trt <- vector()
  
  for (ns in 1:Num_simulation) {
    
    ##Generate simulated data
    source("../datgen.R")
    
    ##Special processing for pooling
    ##Change the id in historical control arm
    historical$id <- historical$id + nlevels(as.factor(current$id))
    historical <- cbind(historical[,1:2], 0, historical[,3:4])
    colnames(historical)[3] <- "interaction"
    ##Pool the current data and historical control arm
    pooling_data <- rbind(current, historical)
    
    data_pooling <- sourceToList("Pooling_Stan_data.R")
    
    init <- list(beta = c(2, 1, as.numeric(sc>7)*interaction_time_trt), Omega = diag(1, 2),
                 tau = rep(0.5, 2), sigma = 1)
    
    ##sampling
    sample_pooling <- sampling(sm_pooling, data = data_pooling, chains = Num_chains, 
                               init = list(init, init, init, init),
                               iter = Num_iter, warmup=Num_iter*perc_burnin, refresh = 0, cores = Num_chains)

  ##Extract the results from the samples
  ##Transform Stanfit to data frames and summarize parameters of interest
  beta_trt <- as.data.frame(sample_pooling)$"beta[3]"
  
  ##The type I error or statistical power for the treatment effect
  power_beta_trt[ns] <- 
      (quantile(beta_trt, 0.025) < 0 & quantile(beta_trt, 0.975) > 0) 
  
  ##Posterior mean of the treatment effect
  mean_beta_trt[ns] <- mean(beta_trt)
  ##Posterior standard deviation of the treatment effect
  sd_beta_trt[ns] <- sd(beta_trt)
  
  ##Progress
  print(paste0("scenario: ", sc, ", iteration: ", ns))
  }
  sample_current_scenarios[[sc]] <- cbind(power_beta_trt, mean_beta_trt, sd_beta_trt)
}
end <- Sys.time()
end - start

##Summary of the results
summary_pooling <- data.frame(power = rep(NA, Num_scenarios), bias = NA, se = NA, mse = NA)
for (sc in 1:Num_scenarios) {
  res <- sample_pooling_scenarios[[sc]]
  ##Type I error rate/Power
  summary_pooling[sc, "power"] <- round((1 - mean(res$power_beta_trt)) * 100, 1)
  ##Bias
  bias <- res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt
  summary_pooling[sc, "bias"] <- paste0(round(mean(bias), 3), " (", round(quantile(bias, 0.025), 3), ", ",  round(quantile(bias, 0.975), 3),")") 
  ##SE
  summary_pooling[sc, "se"] <- paste0(round(mean(res$sd_beta_trt), 3), " (", round(quantile(res$sd_beta_trt, 0.025), 3), ", ",  round(quantile(res$sd_beta_trt, 0.975), 3),")") 
  ##MSE
  mse <- (res$mean_beta_trt - as.numeric(sc>7)*interaction_time_trt)^2
  summary_pooling[sc, "mse"] <- paste0(round(mean(mse), 3), " (", round(quantile(mse, 0.025), 3), ", ",  round(quantile(mse, 0.975), 3),")")
}