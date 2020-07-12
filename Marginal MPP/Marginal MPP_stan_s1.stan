//data block
data {
  int<lower=0> N0;                      // number of historical observations
  int<lower=1> K0;                      // number of fixed effects
  int<lower=1> J0;                      // number of historical subjects 
  int<lower=1> L;                       // number of random effects
  int<lower=1> sub_index0[2*J0];       // index numbers for historical subjects
  matrix[N0, K0] x0;                    // individual predictors
  matrix[N0, L] z0;                     // subject-specific predictors
  vector[N0] y0;                        // outcomes
  real<lower=0> power;                  // power parameter
}
//parameter block
parameters {
  vector[K0] beta;                      // fixed effects
  corr_matrix[L] Omega;                 // correlation matrix of random effects
  vector<lower=0>[L] tau;               // standard deviations of random effects
  real<lower=0> sigmasq;                // error variance
}
//model block
model {
  //prior for fixed effects
  beta ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  //priors for the correlation matrix of the random effects
  Omega ~ lkj_corr(1);
  //prior for standard deviations of the random effects
  tau ~ normal(0, 1);
  //prior for sigmasq
  sigmasq ~ inv_gamma(0.01, 0.01);
  //sampling
  for (j in 1:J0) {
    int sub_first;
    int sub_last;
    sub_first = sub_index0[2*j-1];
    sub_last = sub_index0[2*j];
    target += power * 
    multi_normal_lpdf(y0[sub_first : sub_last] |  
    x0[sub_first : sub_last] * beta, 
    z0[sub_first : sub_last] * quad_form_diag(Omega, tau) * z0[sub_first : sub_last]' + 
    diag_matrix(rep_vector(sigmasq, sub_last - sub_first + 1)));
  }
}
//generated quantities block
generated quantities {
  vector[J0] loglik_fixedpow_sub;
  real loglik_fixedpow;
  for (j in 1:J0) {
    int sub_first;
    int sub_last;
    sub_first = sub_index0[2*j-1];
    sub_last = sub_index0[2*j];
    loglik_fixedpow_sub[j] = multi_normal_lpdf(y0[sub_first : sub_last] |  
    x0[sub_first : sub_last] * beta, 
    z0[sub_first : sub_last] * quad_form_diag(Omega, tau) * z0[sub_first : sub_last]' + 
    diag_matrix(rep_vector(sigmasq, sub_last - sub_first + 1))); 
  }
  loglik_fixedpow = sum(loglik_fixedpow_sub);
}
