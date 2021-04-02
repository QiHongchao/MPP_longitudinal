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

##Func 3: define the function to generate simulated data
source("datgen_func.R")

##Simulation scenarios
num_measures <- 6
num_subject <- 100
##Interaction effect between time and treatment
interaction_time_trt <- 0.36
simulation_scenarios <- expand.grid(trt_eff = c(0, 1) * interaction_time_trt, heterogeneity = 0:6)
simulation_scenarios <- simulation_scenarios[order(simulation_scenarios$trt_eff),]
Num_scenarios <- nrow(simulation_scenarios)

##Simulation settings
Num_chains <- 4
perc_burnin <- 0.5
Num_iter <- 2000
Num_simulation <- 500
set.seed(1104)
seeds_scenarios <- sample.int(.Machine$integer.max, Num_scenarios)