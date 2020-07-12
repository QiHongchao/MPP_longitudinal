##This code aims to prepare for all methods including 1) invoking R packages, 2) Rstan configuration,
##3) functions definition, 4) setting simulation parameters.

##Packages
library(mvtnorm)
library(rstan)
library(coda)
library(parallel)
library(bayesplot)

##Configuration of RStan
options(mc.cores = detectCores(logical = FALSE))
devAskNewPage(ask = FALSE)

##Define functions in advance
##Func 1: Import data and initial values
sourceToList <- function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}

##Func 2: specify index for each subject
sub_index_specification <- function(data) {
  sub_index <- vector()
  data1 <- data
  data1$index <- 1: nrow(data1)
  for (i in 1: nlevels(as.factor(data1$id))) {
    sub_index[2*i-1] <- min(data1$index[data1$id==i])
    sub_index[2*i] <- max(data1$index[data1$id==i])
  }
  return(sub_index)
}

##Func4: define the function to generate simulated data
##Interaction effect between time and treatment
interaction_time_trt <- 0.36

##Number of time points, current/historical fixed effects, random effects, subjects per arm
num_tp <- 6
num_group <- 2
num_fe_current <- 3
num_fe_hist <- 2
num_re <- 2
beta <- c(2, 1)
var_b <- c(0.25, 0.25)
var_e <- 1
nsubpa <- 100

##Generation of current data
datagen_current<-function(num_tp, num_group, num_fe, num_re, beta, var_b, var_e, num_control, num_treat){
  ##with treatment effect or not
  if (trt_eff == 0) {trt <- 0}
  if (trt_eff == 1) {trt <- interaction_time_trt}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,num_tp*num_group,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-rep(seq(0, 1, 1/(num_tp-1)),num_group) ##time point start from 0, equally spaced, time interval is 1
  dm_fe[,3]<-dm_fe[,2]*c(rep(0,num_tp),rep(1,num_tp))
  
  ##Design matrix for random effects
  dm_re<-matrix(NA,num_tp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(num_tp-1))
  
  ##Control arm
  control<-rmvnorm(num_control,dm_fe[1:num_tp,]%*%c(beta, trt), dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e,num_tp))
  dm_fe_control <- data.frame(int = rep(1, num_tp * num_control), 
                              time = rep(seq(0, 1, 1/(num_tp-1)),num_control),
                              interaction = rep(0, num_tp * num_control))
  ##Response for the control arm
  for (i in 1:num_control){
    dm_fe_control[(num_tp * (i-1) + 1):(num_tp * i), "response"] <- control[i,]
  }
  
  ##Treatment arm
  treat <- rmvnorm(num_treat,dm_fe[(num_tp+1):(num_tp*2),]%*%c(beta, trt), dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e,num_tp))
  dm_fe_treat <- data.frame(int = rep(1, num_tp * num_treat), 
                            time = rep(seq(0, 1, 1/(num_tp-1)),num_treat),
                            interaction = rep(seq(0, 1, 1/(num_tp-1)),num_treat))
  ##Response for the treatment arm
  for (i in 1:num_treat){
    dm_fe_treat[(num_tp*(i-1)+1):(num_tp*i), "response"] <- treat[i,]
  }
  
  ##Combine the control and treatment arms
  current <- rbind(dm_fe_control,dm_fe_treat)
  current[, "id"] <- sort(rep(1:(num_control + num_treat), num_tp))
  
  ##Return the dataset
  return(current)
}

##Levels of heterogeneity: 0:"No", 1: "Low+RI", 2: "Moderate+RI", 3: "High+RI",
##4:"Low+RIS", 5: "Moderate+RIS", 6: "High+RIS"
datagen_historical<-function(num_tp, num_fe, num_re, beta, var_b, var_e, num_control){
  ##Covariance matrix of the shared fixed effects
  if (heterogeneity == 0) {cov_study <- matrix(rep(0, 4), 2)}
  if (heterogeneity == 1) {cov_study <- matrix(c(0.01, 0, 0, 0), 2)}
  if (heterogeneity == 2) {cov_study <- matrix(c(0.09, 0, 0, 0), 2)}
  if (heterogeneity == 3) {cov_study <- matrix(c(0.16, 0, 0, 0), 2)}
  if (heterogeneity == 4) {cov_study <- diag(0.01, 2)}
  if (heterogeneity == 5) {cov_study <- diag(0.09, 2)}
  if (heterogeneity == 6) {cov_study <- diag(0.16, 2)}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,num_tp,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-seq(0, 1, 1/(num_tp-1)) ##time point start from 0, equally spaced, time interval is 1/5
  
  ##Design matrix for random effects
  dm_re<-matrix(NA,num_tp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(num_tp-1))
  
  ##Control arm
  control <- rmvnorm(num_control,dm_fe[1:num_tp,]%*%beta + dm_fe[1:num_tp,]%*%as.numeric(rmvnorm(1, rep(0,ncol(dm_fe)), cov_study)), 
                     dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e, num_tp))
  dm_fe_control <- data.frame(int = rep(1, num_tp * num_control), 
                              time = rep(seq(0, 1, 1/(num_tp-1)),num_control))
  ##Response for the historical control
  for (i in 1:num_control){
    dm_fe_control[(num_tp*(i-1)+1):(num_tp*i), "response"] <- control[i,]
  }
  
  ##Historical id
  dm_fe_control[, "id"] <- sort(rep(1:num_control,num_tp))
  
  ##Return the historical controls
  return(dm_fe_control)
}

##simulation scenarios generation
num_scenarios <- 14
simulation_scenarios <- data.frame(trt = rep(0:1, each = num_scenarios/2),
                                   heterogeneity = rep(0:(num_scenarios/2 - 1), 2))

##Simulation settings for CP, no borrowing and pooling
num_iter <- 500
num_chains <- 4
perc_burnin <- 0.5
num_simulation <- 500

##The candidate power parameters in step 1 of the MPP
interval_width <- 0.02
power_all <- seq(0, 1, interval_width)
num_power <- length(power_all)

##Sampling setting for the MPP
num_iter_s1 <- 100
perc_burnin_s1 <- 0.2
num_chains_s2 <- 4
num_iter_s2 <- 500
perc_burnin_s2 <- 0.5

##Compile stan files for the methods
##The marginal MPP
sm1_marg <- stan_model("./Marginal MPP/Marginal MPP_stan_s1.stan")
sm2_marg <- stan_model("./Marginal MPP/Marginal MPP_stan_s2.stan")

##The conditional MPP
sm1_cond <- stan_model("./Conditional MPP/Conditional MPP_stan_s1.stan")
sm2_cond_fixedpow <- stan_model("./Conditional MPP/Conditional/MPP_stan_warmup_s1.stan")
sm2_cond <- stan_model("./Conditional MPP/Conditional MPP_stan_s2.stan")

##The commensurate prior
sm_cp <- stan_model("./Commensurate prior/CP_stan.stan")

##No borrowing
sm_current <- stan_model("./No borrowing/Current_stan.stan")

##Pooling
sm_pooling <- stan_model("./Pooling/Pooling_stan.stan")

##The generation of random seeds for each scenario
set.seed(1234)
seeds_scenarios <- sample.int(.Machine$integer.max, num_scenarios)