//marginal likelihood raised to power, step 2
//function block
functions {
  //linear interpolation of log scaling constant
  real linear_interpolation(matrix logsc, real power, int num_power) {
    // power value smaller than candidate and its corresponding logsc
    real x1; real y1;
    // power value larger than candidate and its corresponding logsc
    real x2; real y2;
    //row number of x1 and x2
    int np;
    //linear interpolation
    real logsc_new;

    for (i in 1: num_power) {
      if (logsc[i, 1] < power) np=i;
    }
    //x and y, column 1 is power value, column 2 is log scaling constant
    x1 = logsc[np, 1];
    y1 = logsc[np, 2];
    x2 = logsc[np+1, 1];
    y2 = logsc[np+1, 2];
    //new log scaling constant
    logsc_new = (y2-y1)/(x2-x1)*(power-x1)+y1;
    return logsc_new;
  }
  //construction of power prior distribution or current likelihood (power=1, data is current data), sub_index0 is a vector of integers,
  //its definition refers to "Dimensionality Declaration" in "User-Defined Functions" of the reference manual.
   real prior_likelihood_calc(int J, int[] sub_index, vector y, matrix x, matrix z, vector beta, real power, matrix Omega, vector tau, real sigmasq) {
   //power prior distribution
   vector[J] prior_likelihood_sub;
   real prior_likelihood;
   int sub_index_lower;
   int sub_index_upper;
   for (j in 1: J) {
    sub_index_lower = sub_index[2*j-1];
    sub_index_upper = sub_index[2*j];
    prior_likelihood_sub[j] = power *
    multi_normal_lpdf( y[sub_index_lower : sub_index_upper] | x[sub_index_lower : sub_index_upper] * beta,
    z[sub_index_lower : sub_index_upper]*quad_form_diag(Omega, tau)*z[sub_index_lower : sub_index_upper]'+
    diag_matrix(rep_vector(sigmasq, sub_index_upper-sub_index_lower+1)));
   }
   prior_likelihood = sum(prior_likelihood_sub);
   return prior_likelihood;
 }
}
//In this block, we want to define a linear interpolation function to calculate log scaling constant.
//data block
data {
  int<lower=0> N;              // num observations, current data
  int<lower=0> N0;             // num observations, historical data
  int<lower=1> K;              // num ind predictors, current data
  int<lower=1> K0;             // num ind predictors, historical data
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
  int<lower=0> num_power;    // num power values
  matrix[num_power, 2] logsc;  // matrix of power values and their corresponding log scaling constant
}
//parameter block
parameters {
  //fixed effects
  vector[K0] beta_mutual;           // regression coeffs for control and treatment group
  real beta_trt;                    // regression coeff only for current treatment group
   //covariance structure for random effects
  corr_matrix[L] Omega;        // prior correlation of random effects
  vector<lower=0>[L] tau;      // prior scale of random effects
  //error term
  real<lower=0> sigmasq;         // attention! this is the variance of the error term
  // Power parameter
  real<lower=0, upper=1> power;
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
  //distribution of power parameter
  power ~ beta(1, 1);
  //construct the likelihood
  target += 
  prior_likelihood_calc(J, sub_index, y, x, z, append_row(beta_mutual, beta_trt), 1, Omega, tau, sigmasq) +  //current likelihood
            prior_likelihood_calc(J0, sub_index0, y0, x0, z0, beta_mutual, power, Omega, tau, sigmasq) -  //power prior distribution without normalizing
            linear_interpolation(logsc, power, num_power);  //log scaling constant
}
