//function block, calculate the likelihood
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
        multi_normal_lpdf(y[sub_index_lower : sub_index_upper] | x[sub_index_lower : sub_index_upper] * beta,
                          z[sub_index_lower : sub_index_upper]*quad_form_diag(Omega, tau)*z[sub_index_lower : sub_index_upper]' +
                          diag_matrix(rep_vector(sigmasq, sub_index_upper-sub_index_lower+1)));
   }
   likelihood = sum(likelihood_sub);
   return likelihood;
 }
}
//data block
data {
  int<lower=0> N;              // number of observations, current data
  int<lower=0> N0;             // number of observations, historical control
  int<lower=1> K;              // number of fixed effects, current data
  int<lower=1> K0;              // number of fixed effects, historical control
  int<lower=1> J;              // number of subjects, current data
  int<lower=1> J0;             // number of subjects, historical control
  int<lower=1> L;              // number of random effects, both current and historical data
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds of index for each subject, current data
  int<lower=1> sub_index0[2*J0] ;      //lower and upper bounds of index for each subject, historical control
  matrix[N, K] x;               // design matrix of the fixed effects, current data
  matrix[N, L] z;          // design matrix of the random effects, current data
  vector[N] y;             // response, current data
  matrix[N0, K0] x0;               // design matrix of the fixed effects, historical control
  matrix[N0, L] z0;          // design matrix of the random effects, historical control
  vector[N0] y0;             // response, historical control
  matrix[K0, K0] var_beta01;   // g*inv(X'X) element for the g-prior of the intercept and time effect
  real<lower=0> var_betatrt;   // g*inv(X'X) element for the g-prior of treatment effect
}
//parameter block
parameters {
  //fixed effects for historical control
  vector[K0] beta_01_hist;           // the intercept and time effect for historical control
  //fixed effects for current study
  vector[K0] beta_01_current;           // the intercept and time effect for current data 
  real beta_trt;            // the treatment effect for current data
  //variance-covriance matrix of beta_01_current given beta_01_hist
  corr_matrix[K0] Omega_beta;        // prior correlation of current beta01 given historical beta01
  vector<lower=0>[K0] tau_beta;      // prior scale of current beta01 given historical beta01
  //covariance structure for random effects in historical control
  corr_matrix[L] Omega0;        // prior correlation of random effects
  vector<lower=0>[L] tau0;      // prior scale of random effects
  //covariance structure for random effects in current study
  corr_matrix[L] Omega1;        // prior correlation of random effects
  vector<lower=0>[L] tau1;      // prior scale of random effects
  //error term for historical control
  real<lower=0> sigma0;        // sd of error term
  //error term for current study
  real<lower=0> sigma1;         // sd of error term
}
transformed parameters {
  real sigmasq0;
  real sigmasq1;
  matrix[K0, K0] var_beta01_g;
  real sd_betatrt_g;
  sigmasq0 = sigma0^2;
  sigmasq1 = sigma1^2;
  var_beta01_g = sigmasq0 * var_beta01;
  sd_betatrt_g = sqrt(sigmasq1 * var_betatrt);
}
//model block
model {
  //g-prior for beta in historical control
  beta_01_hist ~ multi_normal(rep_vector(0, K0), var_beta01_g);
  beta_01_current ~ multi_normal(beta_01_hist, quad_form_diag(Omega_beta, tau_beta));
  //g-prior for the treatment effect
  beta_trt ~ normal(0, sd_betatrt_g);
  //prior for Omega_beta and tau_beta
  tau_beta ~ normal(0, 1);
  Omega_beta ~ lkj_corr(1);
  //priors for G matrix in historical control 
  tau0 ~ normal(0, 1); 
  Omega0 ~ lkj_corr(1);
  //priors for G matrix in current study
  tau1 ~ normal(0, 1); 
  Omega1 ~ lkj_corr(1);
  //prior for sigma in historical control 
  sigma0 ~ normal(0, 2);
  //prior for sigma in current study
  sigma1 ~ normal(0, 2);
  //construct the likelihood
  target += 
  //current likelihood
  likelihood_calc(J, sub_index, y, x, z, append_row(beta_01_current, beta_trt), Omega1, tau1, sigmasq1) + 
  //historical likelihood
  likelihood_calc(J0, sub_index0, y0, x0, z0, beta_01_hist, Omega0, tau0, sigmasq0);  
}
