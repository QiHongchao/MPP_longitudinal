//data block
data {
  int<lower=0> N0;              // num observations
  int<lower=1> K0;              // num ind predictors
  int<lower=1> J0;              // num subjects
  int<lower=1> L;              // num subject-specific predictors
  int<lower=1> sub_index0[2*J0] ;      //lower and upper bounds for each subject
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
  //random effects
  matrix[J0, L] b;                 // vector of random effects
  //error term
  real<lower=0> sigmasq;         // error variance
}
//model block
model {
  //prior for beta
  beta ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  //priors for D matrix
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(1);
  //prior for sigmasq
  sigmasq ~ inv_gamma(0.01, 0.01);
  //sampling 
  //Alternative likelihood, it seems take longer, but if we use this, we don't need transformed parameters block.
  for (j in 1:J0)
  target += power * (multi_normal_lpdf(y0[sub_index0[2*j-1] : sub_index0[2*j]] |
  x0[sub_index0[2*j-1] : sub_index0[2*j]] * beta + z0[sub_index0[2*j-1] : sub_index0[2*j]] * to_vector(b[j, ]),
  diag_matrix(rep_vector(sigmasq, sub_index0[2*j]-sub_index0[2*j-1]+1)))) +
  multi_normal_lpdf(b[j, ] | rep_vector(0, L), quad_form_diag(Omega, tau));
  //In old code, I also raised the distribution of b to power, now I have corrected the error (2019-08-05)
}
generated quantities{
  vector[J0] loglik_fixedpow_sub;
  real loglik_fixedpow;
  for (j in 1:J0)
  loglik_fixedpow_sub[j]= multi_normal_lpdf(y0[sub_index0[2*j-1] : sub_index0[2*j]] |
  x0[sub_index0[2*j-1] : sub_index0[2*j]] * beta + z0[sub_index0[2*j-1] : sub_index0[2*j]] * to_vector(b[j, ]),
  diag_matrix(rep_vector(sigmasq, sub_index0[2*j]-sub_index0[2*j-1]+1)));
  loglik_fixedpow = sum(loglik_fixedpow_sub);
}
