//data block
data {
  int<lower=0> N0;              // num observations
  int<lower=1> K0;              // num ind predictors
  int<lower=1> J;              // num subjects
  int<lower=1> L;              // num subject-specific predictors
  int<lower=1> sub_index0[2*J] ;      //lower and upper bounds for each subject
  matrix[N0, K0] x0;               // individual predictors
  matrix[N0, L] z0;          // subject-specific predictors
  vector[N0] y0;                 // outcomes
  real<lower=0> power;         //power parameter
}
//parameter block
parameters {
  //fixed effects
  vector[K0] beta;           // indiv coeffs by group
  //covariance structure for random effects
  corr_matrix[L] Omega;        // prior correlation of random effects
  vector<lower=0>[L] tau;      // prior scale of random effects
  //error term
  real<lower=0> sigmasq;         // variance of error term
}
//model block
model {
  //prior for fixed effects
  beta ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  //priors for D matrix
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(1);
  //prior for error term variance
  sigmasq ~ inv_gamma(0.01, 0.01);
  //sampling
  for (j in 1:J)
    target += power * 
    multi_normal_lpdf(y0[sub_index0[2*j-1] : sub_index0[2*j]] |  
    x0[sub_index0[2*j-1] : sub_index0[2*j]] * beta, 
    z0[sub_index0[2*j-1] : sub_index0[2*j]]*quad_form_diag(Omega, tau)*z0[sub_index0[2*j-1] : sub_index0[2*j]]' + 
    diag_matrix(rep_vector(sigmasq, sub_index0[2*j]-sub_index0[2*j-1]+1))); 
}
//generated quantities block
generated quantities {
  vector[J] loglik_fixedpow_sub;
  real loglik_fixedpow;
  for (j in 1:J)
  loglik_fixedpow_sub[j]=multi_normal_lpdf(y0[sub_index0[2*j-1] : sub_index0[2*j]] |  
    x0[sub_index0[2*j-1] : sub_index0[2*j]] * beta, 
    z0[sub_index0[2*j-1] : sub_index0[2*j]]*quad_form_diag(Omega, tau)*z0[sub_index0[2*j-1] : sub_index0[2*j]]' + 
    diag_matrix(rep_vector(sigmasq, sub_index0[2*j]-sub_index0[2*j-1]+1))); 
  loglik_fixedpow = sum(loglik_fixedpow_sub);
}
