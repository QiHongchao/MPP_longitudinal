//data block
data {
  int<lower=0> N;              // number of observations
  int<lower=1> K;              // number of fixed effects
  int<lower=1> J;              // number of subjects
  int<lower=1> L;              // number of random effects
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds of index for each subject
  matrix[N, K] x;               // design matrix of the fixed effects
  matrix[N, L] z;          // design matrix of the random effects
  vector[N] y;                 // response, historical control
  real<lower=0> power;         // fixed power parameter
  matrix[K, K] var_beta;   // g*inv(X'X) element for the g-prior
}
//parameter block
parameters {
  //fixed effects
  vector[K] beta;         
  //covariance structure for random effects
  corr_matrix[L] Omega;        // prior for the correlation matrix of random effects
  vector<lower=0>[L] tau;      // prior for scales of random effects
  //random effects
  matrix[J, L] b;                 
  //error term
  real<lower=0> sigma;         // sd of error term
}
transformed parameters {
  real sigmasq;
  matrix[K, K] var_beta_g;
  sigmasq = sigma^2;
  var_beta_g = sigmasq * var_beta;
}
//model block
model {
  //g-prior for beta
  beta ~ multi_normal(rep_vector(0, K), var_beta_g);
  //priors for G matrix
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(1);
  //prior for sigma
  sigma ~ normal(0, 2);
  //likelihood downweighted with fixed power
  for (j in 1:J)
  target += power * (multi_normal_lpdf(y[sub_index[2*j-1] : sub_index[2*j]] |
  x[sub_index[2*j-1] : sub_index[2*j]] * beta + z[sub_index[2*j-1] : sub_index[2*j]] * to_vector(b[j, ]),
  diag_matrix(rep_vector(sigmasq, sub_index[2*j]-sub_index[2*j-1]+1)))) +
  multi_normal_lpdf(b[j, ] | rep_vector(0, L), quad_form_diag(Omega, tau));
}
//generated quantities block, log likelihood given the samples, to calculate log scaling constant
generated quantities{
  vector[J] loglik_fixedpow_sub;
  real loglik_fixedpow;
  for (j in 1:J)
  loglik_fixedpow_sub[j]= multi_normal_lpdf(y[sub_index[2*j-1] : sub_index[2*j]] |
  x[sub_index[2*j-1] : sub_index[2*j]] * beta + z[sub_index[2*j-1] : sub_index[2*j]] * to_vector(b[j, ]),
  diag_matrix(rep_vector(sigmasq, sub_index[2*j]-sub_index[2*j-1]+1)));
  loglik_fixedpow = sum(loglik_fixedpow_sub);
}
