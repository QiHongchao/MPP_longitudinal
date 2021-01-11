##This code aims to prepare for all methods including 1) invoking R packages, 2) Rstan configuration,
##3) functions definition, 4) setting simulation parameters.

##Packages
library(mvtnorm)
library(lme4)
library(rstan)
library(coda)
library(parallel)
library(bayesplot)

##Configuration of RStan
options(mc.cores = parallel::detectCores(logical = FALSE))
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
diff_between_current_hist <- function(heterogeneity, trt) {
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
  if (trt == 0) {
    trt <- "No"
  }
  if (trt == 1) {
    trt <- "Yes"
  }
  return(c(heterogeneity, trt))
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
source("datgen_func.R")

##All the 14 scenarios for commensurate prior method
Num_scenarios <- 14
simulation_scenarios <- matrix(NA, Num_scenarios, 4)
simulation_scenarios[, 1] <- c(rep(0, Num_scenarios/2), rep(1, Num_scenarios/2))  ##without or with treatment effect
simulation_scenarios[, 2] <- 6  ##number of repeated measurements
simulation_scenarios[, 3] <- 100  ##number of subjects per arm
simulation_scenarios[, 4] <- rep(0:6, nlevels(factor(simulation_scenarios[, 1])))  ##levels of heterogeneity
colnames(simulation_scenarios) <- c("trt", "num_measures", "num_subject", "heterogeneity")

##Simulation settings
Num_chains <- 4
perc_burnin <- 0.5
Num_iter <- 2000
Num_simulation <- 500
set.seed(1104)
seeds_scenarios <- sample.int(.Machine$integer.max, Num_scenarios)

##Interaction effect between time and treatment
interaction_time_trt <- 0.36