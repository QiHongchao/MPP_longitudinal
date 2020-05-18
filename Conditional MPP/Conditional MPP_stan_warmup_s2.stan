//function block
functions {
  //construction of current likelihood and power prior
    real likelihood_prior_calc(int J, int L, int[] sub_index, vector y, matrix x, matrix z, vector beta, real power, matrix b, real sigmasq, matrix Omega, vector tau) {
    //power prior distribution
    vector[J] likelihood_prior_sub;
    real likelihood_prior;
    int sub_index_lower;
    int sub_index_upper;
    for (j in 1: J) {
      sub_index_lower = sub_index[2*j-1];
      sub_index_upper = sub_index[2*j];
      likelihood_prior_sub[j] = power * (multi_normal_lpdf( y[sub_index_lower : sub_index_upper] | 
      x[sub_index_lower : sub_index_upper] * beta + z[sub_index_lower : sub_index_upper] * to_vector(b[j, ]),
      diag_matrix(rep_vector(sigmasq, sub_index_upper-sub_index_lower+1)))) +
      multi_normal_lpdf(b[j, ] | rep_vector(0, L), quad_form_diag(Omega, tau));
    }
    likelihood_prior = sum(likelihood_prior_sub);
    return likelihood_prior;
  }
}
//In this block, we want to define a linear interpolation function to calculate log scaling constant.
//data block
data {
  int<lower=0> N;              // num observations, current data
  int<lower=0> N0;             // num observations, historical data
  int<lower=1> K;              // num ind predictors, current data
  int<lower=1> K0;              // num ind predictors, historical data
  int<lower=1> J;              // num subjects, current data
  int<lower=1> J0;             // num subjects, historical data
  int<lower=1> L;              // num subject-specific predictors, both current and historical data
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds for each subject, current data
  int<lower=1> sub_index0[2*J0] ;      //lower and upper bounds for each subject, historical data
  matrix[N, K] x;               // individual predictors, current data
  matrix[N, L] z;          // subject-specific predictors, current data
  vector[N] y;             // response, current data
  matrix[N0, K0] x0;               // individual predictors, historical data
  matrix[N0, L] z0;          // subject-specific predictors, historical data
  vector[N0] y0;             // response, historical data
}
//parameter block
parameters {
  //fixed effects
  vector[K0] beta_mutual;           // regression coefficients for control and treatment group
  real beta_trt;            // regression coefficients only for current treatment group
  //covariance structure for random effects
  corr_matrix[L] Omega;
  vector<lower=0>[L] tau;
  //random effects
  matrix[J0, L] b0;               // historical random effects
  matrix[J, L] b;                 //current random effects
  //error term
  real<lower=0> sigmasq;         // attention! this is the variance of the error term
}
//model block
model {
  //prior for beta
  beta_mutual ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  beta_trt ~ normal(0, 10);
  //priors for D matrix
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(1);
  //prior for sigmasq
  sigmasq ~ inv_gamma(0.01, 0.01);
  //construct the likelihood
  target += 
    //current likelihood
  likelihood_prior_calc(J, L, sub_index, y, x, z, append_row(beta_mutual, beta_trt), 1, b, sigmasq, Omega, tau) + 
  //power prior distribution without normalizing
  likelihood_prior_calc(J0, L, sub_index0, y0, x0, z0, beta_mutual, 0.5, b0, sigmasq, Omega, tau);
}
