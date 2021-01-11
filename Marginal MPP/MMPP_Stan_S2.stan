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
  //construction of power prior distribution or current likelihood (power=1, data is current data)
   real prior_likelihood_calc(int J, int[] sub_index, vector y, matrix x, matrix z, vector beta, real power, matrix Omega, vector tau, real sigmasq) {
   //power prior distribution
   vector[J] prior_likelihood_sub;
   real prior_likelihood;
   int sub_index_lower;
   int sub_index_upper;
   for (j in 1: J) {
    sub_index_lower = sub_index[2*j-1];
    sub_index_upper = sub_index[2*j];
    prior_likelihood_sub[j] = power * multi_normal_lpdf(y[sub_index_lower : sub_index_upper] | x[sub_index_lower : sub_index_upper] * beta,
    z[sub_index_lower : sub_index_upper]*quad_form_diag(Omega, tau)*z[sub_index_lower : sub_index_upper]'+
    diag_matrix(rep_vector(sigmasq, sub_index_upper-sub_index_lower+1)));
   }
   prior_likelihood = sum(prior_likelihood_sub);
   return prior_likelihood;
 }
}
//data block
data {
  int<lower=0> N;              // number of observations, current data
  int<lower=0> N0;             // number of observations, historical control
  int<lower=1> K;              // number of fixed effects, current data
  int<lower=1> K0;             // number of fixed effects, historical control
  int<lower=1> J;              // number of subjects, current data
  int<lower=1> J0;             // number of subjects, historical control
  int<lower=1> L;              // number of random effects, both current data and historical control
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds of index for each subject, current data
  int<lower=1> sub_index0[2*J0] ;      //lower and upper bounds of index for each subject, historical control
  matrix[N, K] x;               // design matrix of the fixed effects, current data
  matrix[N, L] z;          // design matrix of the random effects, current data
  vector[N] y;             // response, current data
  matrix[N0, K0] x0;               // design matrix of the fixed effects, historical control
  matrix[N0, L] z0;          // design matrix of the random effects, historical control
  vector[N0] y0;             // response, historical control
  int<lower=0> num_power;    // number of power values
  matrix[num_power, 2] logsc;  // matrix of power values and their corresponding log scaling constant
  matrix[K0, K0] var_betamutual;   // g*inv(X'X) element for the g-prior
  real<lower=0> var_betatrt;   // g*inv(X'X) element for the g-prior
}
//parameter block
parameters {
  //fixed effects
  vector[K0] beta_mutual;           // mutual regression coefficients for current data and historical control
  real beta_trt;                    // treatment effect
  //covariance structure for random effects
  corr_matrix[L] Omega;        // prior for the correlation matrix of random effects
  vector<lower=0>[L] tau;      // prior for scales of random effects
  //error term
  real<lower=0> sigma;         // sd of error term
  // Power parameter
  real<lower=0, upper=1> power;
}
transformed parameters {
  real sigmasq;
  matrix[K0, K0] var_betamutual_g;
  real sd_betatrt_g;
  sigmasq = sigma^2;
  var_betamutual_g = sigmasq * var_betamutual;
  sd_betatrt_g = sqrt(sigmasq * var_betatrt);
}
//model block
model {
  //g-prior for beta
  beta_mutual ~ multi_normal(rep_vector(0, K0), var_betamutual_g);
  //g-prior for treatment effect
  beta_trt ~ normal(0, sd_betatrt_g);
  //priors for G matrix
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(1);
  //prior for sigma
  sigma ~ normal(0, 2);
  //prior for power parameter
  power ~ beta(1, 1);
  //construct the likelihood
  target += 
  prior_likelihood_calc(J, sub_index, y, x, z, append_row(beta_mutual, beta_trt), 1, Omega, tau, sigmasq) +  //current likelihood
            prior_likelihood_calc(J0, sub_index0, y0, x0, z0, beta_mutual, power, Omega, tau, sigmasq) -  //power prior without normalizing
            linear_interpolation(logsc, power, num_power);  //log scaling constant
}
