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
  int<lower=0> N;              // num observations, pooled data
  int<lower=1> K;              // num ind predictors, pooled data
  int<lower=1> J;              // num subjects, pooled data
  int<lower=1> L;              // num subject-specific predictors, pooled data
  int<lower=1> sub_index[2*J] ;      //lower and upper bounds for each subject, pooled data
  matrix[N, K] x;               // individual predictors, pooled data
  matrix[N, L] z;          // subject-specific predictors, pooled data
  vector[N] y;             // response, pooled data
}
//parameter block
parameters {
  //fixed effects
  vector[K] beta;         
  //covariance structure for random effects in pooled data
  corr_matrix[L] Omega;        // prior correlation of random effects
  vector<lower=0>[L] tau;      // prior scale of random effects
  //error term for pooled study
  real<lower=0> sigmasq;         // attention! this is the variance of the error term
}
//model block
model {
  //prior for beta in pooled data
  beta ~ multi_normal(rep_vector(0, K), diag_matrix(rep_vector(100, K)));
  //priors for D matrix
  tau ~ normal(0, 1); //maybe we can also specify a half normal normal(0, 1)
  Omega ~ lkj_corr(1); //larger parameter, stronger prior for identity matrix
  //prior for sigmasq 
  sigmasq ~ inv_gamma(0.01, 0.01);
  //construct the likelihood
  target += 
  //pooled likelihood
  likelihood_calc(J, sub_index, y, x, z, beta, Omega, tau, sigmasq);
}
