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
##Func1: Import data and initial values
sourceToList <- function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}

##Func2: Heterogeneity between current and historical datasets and w/ or w/o treatment effect
trt_effect <- function(trt) {
  if (trt == 0) {
    trt <- "No"
  } else {
    trt <- "Yes"
  }
  return(trt)
}

diff_between_current_hist <- function(heterogeneity) {
  if (heterogeneity==0) {
    heterogeneity<-"No"
  } 
  if (heterogeneity==1) {
    heterogeneity<-"Low+RI"
  }
  if (heterogeneity==2) {
    heterogeneity<-"Moderate+RI"
  }
  if (heterogeneity==3) {
    heterogeneity<-"High+RI"
  }
  if (heterogeneity==4) {
    heterogeneity<-"Low+RIS"
  }
  if (heterogeneity==5) {
    heterogeneity<-"Moderate+RIS"
  }
  if (heterogeneity==6) {
    heterogeneity<-"High+RIS"
  }
  return(heterogeneity)
}

##Func3: specify index for each subject
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
##In oder to simplify the case, I simulate an RCT with equal number of patients in treatment and control arm, respectively.
##Design matrix: we assume current and historical study have the same design matrix.
##Design matrix for fixed effects and random effects
design_matrix_gen_current<-function(ntp, group, num_fe, num_re, beta, var_b, var_e, num_control, num_treat){
  ##with treatment effect or not
  if (trt_eff == "No") {trt <- 0}
  if (trt_eff == "Yes") {trt <- interaction_time_trt}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,ntp*group,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-rep(seq(0, 1, 1/(ntp-1)),group) ##time point start from 0, equally spaced, time interval is 1
  dm_fe[,3]<-dm_fe[,2]*c(rep(0,ntp),rep(1,ntp))
  
  ##Design matrix for random effects
  dm_re<-matrix(NA,ntp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(ntp-1))
  
  ##Control arm
  control<-rmvnorm(num_control,dm_fe[1:ntp,]%*%c(beta, trt), dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e,ntp))
  dm_fe_control<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_control[,1]<-1
  dm_fe_control[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_control)
  dm_fe_control[,3]<-0
  for (i in 1:num_control){
    dm_fe_control[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-control[i,]
  }
  ##Treatment arm
  treat<-rmvnorm(num_treat,dm_fe[(ntp+1):(ntp*2),]%*%c(beta, trt), dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e,ntp))
  dm_fe_treat<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_treat[,1]<-1
  dm_fe_treat[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_treat)
  dm_fe_treat[,3]<-dm_fe_treat[,2]
  for (i in 1:num_treat){
    dm_fe_treat[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-treat[i,]
  }
  ##Response
  response<-rbind(dm_fe_control,dm_fe_treat)
  response[,num_fe+2]<-sort(rep(1:(num_control+num_treat),ntp))
  colnames(response)<-c("int","time","interaction","response","id")
  ##Return the dataset
  return(as.data.frame(response))
}

design_matrix_gen_hc<-function(ntp, num_fe, num_re, beta, var_b, var_e, num_control){
  ##standard deviation of time trend
  if (heterogeneity == "No") {cov_study <- matrix(rep(0, 4), 2)}
  if (heterogeneity == "Low+RI") {cov_study <- matrix(c(0.01, 0, 0, 0), 2)}
  if (heterogeneity == "Moderate+RI") {cov_study <- matrix(c(0.09, 0, 0, 0), 2)}
  if (heterogeneity == "High+RI") {cov_study <- matrix(c(0.16, 0, 0, 0), 2)}
  if (heterogeneity == "Low+RIS") {cov_study <- diag(0.01, 2)}
  if (heterogeneity == "Moderate+RIS") {cov_study <- diag(0.09, 2)}
  if (heterogeneity == "High+RIS") {cov_study <- diag(0.16, 2)}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,ntp,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-seq(0, 1, 1/(ntp-1)) ##time point start from 0, equally spaced, time interval is 1/5
  
  ##Design matrix for random effects
  dm_re<-matrix(NA,ntp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(ntp-1))
  
  ##Control arm
  control<-rmvnorm(num_control,dm_fe[1:ntp,]%*%beta + dm_fe[1:ntp,]%*%as.numeric(rmvnorm(1, rep(0,ncol(dm_fe)), cov_study)), 
                   dm_re%*%(diag(var_b))%*%t(dm_re) + diag(var_e, ntp))
  dm_fe_control<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_control[,1]<-1
  dm_fe_control[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_control)
  for (i in 1:num_control){
    dm_fe_control[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-control[i,]
  }
  ##Response
  response<-dm_fe_control
  response[,num_fe+2]<-sort(rep(1:num_control,ntp))
  colnames(response)<-c("int","time","response","id")
  ##Return the dataset
  return(as.data.frame(response))
}

##simulation scenarios generation
num_scenarios <- 14
ntp <- 6
nsubpa <- 100
simulation_scenarios <- data.frame(trt = rep(0:1, each = num_scenarios/2), ntp = ntp, nsub = nsubpa,
                                   heterogeneity = rep(0:(num_scenarios/2 - 1), 2))

##Simulation settings
num_iter <- 500
num_chains <- 4
perc_burnin <- 0.5
num_simulation <- 500
set.seed(1234)
seeds_scenarios <- sample.int(.Machine$integer.max, num_scenarios)

##Interaction effect between time and treatment
interaction_time_trt <- 0.36

##The candidate power parameters in step 1: 0 to 1 by interval width
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
sm2_cond_fixedpow <- stan_model("./Conditional MPP/Conditional MPP_stan_warmup_s1.stan")
sm2_cond <- stan_model("./Conditional MPP/Conditional MPP_stan_s2.stan")
##The commensurate prior
sm_cp <- stan_model("./Commensurate prior/CP_stan.stan")
##No borrowing
sm_current <- stan_model("./No borrowing/Current_stan.stan")
##Pooling
sm_pooling <- stan_model("./Pooling/Pooling_stan.stan")