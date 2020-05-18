//function block
functions {
  real likelihood_calc(int J, int[] sub_index, vector y, matrix x, matrix z, vector beta, matrix Omega, vector tau, real sigmasq) {
    //likelihood
    vector[J] likelihood_sub;
    real likelihood;
    int sub_index_lower;
    int sub_index_upper;
    for (j in 1: J) {
      sub_index_lower = sub_index[2*j-1];
      sub_index_upper = sub_index[2*j];
      likelihood_sub[j] =
        multi_normal_lpdf( y[sub_index_lower : sub_index_upper] | x[sub_index_lower : sub_index_upper] * beta,
                           z[sub_index_lower : sub_index_upper]*quad_form_diag(Omega, tau)*z[sub_index_lower : sub_index_upper]' +
    diag_matrix(rep_vector(sigmasq, sub_index_upper-sub_index_lower+1)));
   }
   likelihood = sum(likelihood_sub);
   return likelihood;
 }
}
//data block
data {
  int<lower=0> N;              // num observations, current data
  int<lower=0> N0;             // num observations, historical control arm
  int<lower=1> K;              // num ind predictors, current data
  int<lower=1> K0;              // num ind predictors, historical control arm
  int<lower=1> J;              // num subjects, current data
  int<lower=1> J0;             // num subjects, historical control arm
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
  //fixed effects for historical control arm
  vector[K0] beta_01_hist;           // regression coefficients for historical control arm
  //fixed effects for current study
  vector[K0] beta_01_current;           // regression coefficients for current data intercept and time effect
  real beta_trt;            // regression coefficient only for current treatment group
  //variance-covriance matrix of beta_01_current given beta_01_hist
  corr_matrix[K0] Omega_beta;        // prior correlation of current beta01 given historical beta01
  vector<lower=0>[K0] tau_beta;      // prior scale of current beta01 given historical beta01
  //covariance structure for random effects in historical control arm
  corr_matrix[L] Omega0;        // prior correlation of random effects
  vector<lower=0>[L] tau0;      // prior scale of random effects
  //covariance structure for random effects in current study
  corr_matrix[L] Omega1;        // prior correlation of random effects
  vector<lower=0>[L] tau1;      // prior scale of random effects
  //error term for historical control arm
  real<lower=0> sigmasq0;         // attention! this is the variance of the error term
  //error term for current study
  real<lower=0> sigmasq1;         // attention! this is the variance of the error term
}
//model block
model {
  //prior for beta in historical control armand current data
  beta_01_hist ~ multi_normal(rep_vector(0, K0), diag_matrix(rep_vector(100, K0)));
  beta_01_current ~ multi_normal(beta_01_hist, quad_form_diag(Omega_beta, tau_beta));
  beta_trt ~ normal(0, 10);
  //prior for Omega_beta and tau_beta
  tau_beta ~ normal(0, 1);
  Omega_beta ~ lkj_corr(1);
  //priors for D matrix in historical control arm
  tau0 ~ normal(0, 1); //maybe we can also specify a half normal normal(0, 1)
  Omega0 ~ lkj_corr(1); //larger parameter, stronger prior for identity matrix
  //priors for D matrix in current study
  tau1 ~ normal(0, 1); //maybe we can also specify a half normal normal(0, 1)
  Omega1 ~ lkj_corr(1); //larger parameter, stronger prior for identity matrix
  //prior for sigmasq in historical control arm
  sigmasq0 ~ inv_gamma(0.01, 0.01);
  //prior for sigmasq in current study
  sigmasq1 ~ inv_gamma(0.01, 0.01);
  //construct the likelihood
  target += 
  //current likelihood
  likelihood_calc(J, sub_index, y, x, z, append_row(beta_01_current, beta_trt), Omega1, tau1, sigmasq1) + 
  //historical likelihood
  likelihood_calc(J0, sub_index0, y0, x0, z0, beta_01_hist, Omega0, tau0, sigmasq0);  
}
