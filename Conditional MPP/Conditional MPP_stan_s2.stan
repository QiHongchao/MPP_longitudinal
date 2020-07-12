//function block
functions {
  //function for linear interpolation of the log scaling constant
  real linear_interpolation (matrix logsc, real power, int num_power) {
  // max fixed power less than the proposal and logsc
  real x1; real y1;
  // min fixed power larger than the proposal and logsc
  real x2; real y2;
  // row number of x1
  int np;
  // logsc for the proposal
  real logsc_new;
  // locate the proposal 
  for (i in 1: num_power) {
    if (logsc[i, 1] < power) np = i;
  }
  // logsc, column 1 is power value, column 2 is log scaling constant
  x1 = logsc[np, 1];
  y1 = logsc[np, 2];
  x2 = logsc[np+1, 1];
  y2 = logsc[np+1, 2];
  // new log scaling constant
  logsc_new = (y2-y1)/(x2-x1) * (power-x1) + y1;
  return logsc_new;
  }
  // construction of current likelihood and power prior
  real likelihood_prior_calc(int J, int L, int[] sub_index, vector y, matrix x, 
       vector beta, matrix z, matrix b, matrix Omega, vector tau, real sigmasq, 
       real power) {
  vector[J] likelihood_prior_sub;
  real likelihood_prior;
  for (j in 1: J) {
    int sub_first;
    int sub_last;
    sub_first = sub_index[2*j-1];
    sub_last = sub_index[2*j];
    likelihood_prior_sub[j] = power * 
      (multi_normal_lpdf(y[sub_first : sub_last] | 
      x[sub_first : sub_last] * beta + 
      z[sub_first : sub_last] * to_vector(b[j, ]),
      diag_matrix(rep_vector(sigmasq, sub_last - sub_first+1)))) +
      multi_normal_lpdf(b[j, ] | rep_vector(0, L), quad_form_diag(Omega, tau));
    }
    likelihood_prior = sum(likelihood_prior_sub);
    return likelihood_prior;
  }
}
//data block
data {
  int<lower=0> N;                     // number of current observations
  int<lower=0> N0;                    // number of historical observations
  int<lower=1> K;                     // number of current fixed effects
  int<lower=1> K0;                    // number of historical fixed effects
  int<lower=1> J;                     // number of current subjects
  int<lower=1> J0;                    // number of historical subjects
  int<lower=1> L;                     // number of random effects
  int<lower=1> sub_index[2*J];        // index numbers for current subjects
  int<lower=1> sub_index0[2*J0];      // index numbers for historical subjects
  matrix[N, K] x;                     // design matrix for current fixed effects
  matrix[N, L] z;                     // design matrix for current random effects
  vector[N] y;                        // current response
  matrix[N0, K0] x0;                  // design matrix for historical fixed effects
  matrix[N0, L] z0;                   // design matrix for historical random effects
  vector[N0] y0;                      // historical response
  int<lower=0> num_power;             // number of fixed power values
  matrix[num_power, 2] logsc;         // power values and logsc
}
//parameter block
parameters {
  vector[K0] beta;                    // mutual fixed effects
  real beta_trt;                      // treatment effect
  corr_matrix[L] Omega;               // correlation matrix of random effects
  vector<lower=0>[L] tau;             // standard deviations of random effects 
  matrix[J0, L] b0;                   // historical random effects
  matrix[J, L] b;                     // current random effects
  real<lower=0> sigmasq;              // error variance
  real<lower=0, upper=1> power;       // power parameter
}
//model block
model {
  // prior for mutual fixed effects
  beta ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  // prior for treatment effect
  beta_trt ~ normal(0, 10);
  // priors for the correlation matrix of the random effects
  Omega ~ lkj_corr(1);
  // prior for standard deviations of the random effects
  tau ~ normal(0, 1);
  // prior for sigmasq
  sigmasq ~ inv_gamma(0.01, 0.01);
  // uniform prior for the power parameter
  power ~ beta(1, 1);
  // sampling
  target += 
  // current likelihood
  likelihood_prior_calc(J, L, sub_index, y, x, append_row(beta, beta_trt), 
                        z, b, Omega, tau, sigmasq, 1) + 
  // power prior distribution without normalizing
  likelihood_prior_calc(J0, L, sub_index0, y0, x0, beta, 
                        z0, b0, Omega, tau, sigmasq, power) - 
  // log scaling constant
  linear_interpolation(logsc, power, num_power);  
}
