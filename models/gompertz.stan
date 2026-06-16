//
// Latent Gompertz function modelled with multiplicative noise 
//
data {
  int<lower=1> N_train;  // number of training data
  int<lower=0> T_train;  // number of unique train time points
  array[T_train] real train_time;  // train time
  int<lower=1, upper=T_train> time_train_idx[N_train];  // train time idxs
  int<lower=0> N_b_train;  // number of blanks in training
  int<lower=0, upper=N_train> blank_train_idx[N_b_train];  // blank indices
  int<lower=0> N_auris_train;  // number of non blanks in training
  int<lower=0, upper=N_train> auris_train_idx[N_auris_train];  // non blank indices
  vector[N_train] y_obs_train; // input data of OD values
  
  int<lower=0> N_test;  // number of predictions
  int<lower=0> T_test;  // number of unique test time points
  array[T_test] real test_time;  // test time
  int<lower=1, upper=T_test> time_test_idx[N_test];  // test time idxs
  int<lower=0> N_b_test;  // number of blanks in testing
  int<lower=0, upper=N_test> blank_test_idx[N_b_test];  // blank indices
  int<lower=0> N_auris_test;  // number of non blanks in testing
  int<lower=0, upper=N_test> auris_test_idx[N_auris_test];  // non blank indices
  vector[N_test] y_obs_test; // input data of OD test values

  int<lower=1, upper=N_train+N_test> y_train_idx[N_train];  // training idxs
  int<lower=1, upper=N_train+N_test> y_test_idx[N_test];  // testing idxs
  
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N_log_lik = N_test*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  real<lower=0> IC = (2.5)*1e2; // initial inoculum 
}

parameters {
  real<lower=0> y0;  // initial inoculum parameter
  real<lower=log10(IC)> K_tilde;  // log10 carrying capacity 
  real<lower=log10(IC)> delta_tilde;  // log10 linear transform parameter (scale)
  real<lower=0> beta;  // initial growth rate
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> sigma;  // scale of observed noise
}

transformed parameters {
  vector[N_auris_train] f;
  real<lower=IC> K = pow(10, K_tilde);  // carrying capacity 
  real<lower=IC> delta = pow(10, delta_tilde);  // linear transform parameter (scale)
  for (i in 1:N_auris_train) f[i] = y0*exp(log(K/y0)*(1-exp(-(beta*train_time[time_train_idx[auris_train_idx[i]]]))));  // gompertz function  
}

model {
  y0 ~ lognormal(log(IC), 1);
  sigma ~ normal(0, 0.5);
  K_tilde ~ normal(9, 2);
  beta ~ std_normal();
  basal ~ lognormal(0, 1);
  delta_tilde ~ cauchy(log10(IC), 1);

  if (include_likelihood) {
    y_obs_train[auris_train_idx] ~ lognormal(log(basal + f/delta), sigma);
    y_obs_train[blank_train_idx] ~ lognormal(log(basal), sigma);
  }
}

generated quantities {
  vector[N_auris_test] f_pred;
  real y_tot[N_train+N_test];
  real y_rep[N_train];
  real y_pred[N_test];
  vector[N_log_lik] log_lik;  // log likelihood
  
  y_rep[auris_train_idx] = lognormal_rng(log(basal + f/delta), sigma);
  for (i in blank_train_idx) y_rep[i] = lognormal_rng(log(basal), sigma);
  for (i in 1:N_auris_test) f_pred[i] = y0*exp(log(K/y0)*(1-exp(-(beta*test_time[time_test_idx[auris_test_idx[i]]]))));  // gompertz function  
  y_pred[auris_test_idx] = lognormal_rng(log(basal + f_pred/delta), sigma);
  for (i in blank_test_idx) y_pred[i] = lognormal_rng(log(basal), sigma);

  // calculate log likelihood
  if (include_likelihood) {
    for (i in blank_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal), sigma);
    for (i in 1:N_auris_test)
      log_lik[auris_test_idx[i]] = lognormal_lpdf(y_obs_test[auris_test_idx[i]] | log(basal + f_pred[i]/delta), sigma);
  }
  y_tot[y_train_idx] = y_rep;
  y_tot[y_test_idx] = y_pred;
}
