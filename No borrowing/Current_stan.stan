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
  int<lower=1> K;              // number of fixed effects
  int<lower=1> J;              // number of subjects, current data
  int<lower=1> L;              // number of random effects
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds of index for each subject
  matrix[N, K] x;               // design matrix of the fixed effects
  matrix[N, L] z;          // design matrix of the random effects
  vector[N] y;             // response
  matrix[K, K] var_beta;   // g*inv(X'X) element for the g-prior
}
//parameter block
parameters {
  //fixed effects
  vector[K] beta;         
  //covariance structure for random effects in current data
  corr_matrix[L] Omega;        // prior for the correlation matrix of random effects
  vector<lower=0>[L] tau;      // prior for scales of random effects
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
  //construct the likelihood
  target += 
  likelihood_calc(J, sub_index, y, x, z, beta, Omega, tau, sigmasq);
}
