//
// Latent Edwards function modelled with multiplicative noise 
//

functions {
  vector edwards(real t,
                 vector y,
                 real beta,
                 real L) {
    vector[1] dydt;
    dydt[1] = (beta)*y[1]*exp(-(y[1]/L));  // Edwards ODE
    return dydt;
  }
}

data {
  int<lower=1> N_train;  // number of training data
  int<lower=0> T_train;  // number of unique train time points
  array[T_train] real train_time;  // training time 
  int<lower=1, upper=T_train> time_train_idx[N_train];  // training time idxs
  int<lower=0> N_b_train;  // number of blanks in training
  int<lower=0, upper=N_train> blank_train_idx[N_b_train];  // blank indices
  int<lower=0> N_auris_train;  // number of non blanks in training
  int<lower=0, upper=N_train> auris_train_idx[N_auris_train];  // non blank indices
  vector[N_train] y_obs_train; // input data of OD values
  
  int<lower=0> N_test;  // number of predictions
  int<lower=0> T_test;  // number of unique test time points
  array[T_test] real test_time;  // testing time
  int<lower=1, upper=T_test> time_test_idx[N_test];  // testing time idxs
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
  real<lower=0> y0;  // initial value for auris in wells
  real L_tilde;  // log10 growth-impedance constant
  real<lower=log10(IC)> delta_tilde;  // log10 linear transform parameter (scale)
  real<lower=0> beta;  // growth rate
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> sigma;  // scale of observed noise
}

transformed parameters {
  real<lower=0> L = pow(10, L_tilde); // growth-impedance constant
  real<lower=IC> delta = pow(10, delta_tilde);  // linear transform parameter (scale)
  array[T_train-1] vector[1] mu_hat = ode_rk45(edwards, to_vector({y0}), train_time[1], train_time[2:T_train], beta, L);  // solve edwards ODE
  array[T_train] vector[1] mu;
  mu[1, 1] = y0;
  mu[2:T_train, 1] = mu_hat[, 1];
}

model {
  y0 ~ lognormal(log(IC), 1);
  sigma ~ normal(0, 0.5);
  L_tilde ~ normal(log10(IC), 2);
  beta ~ std_normal();
  basal ~ lognormal(0, 1);
  delta_tilde ~ cauchy(log10(IC), 1);

  if (include_likelihood) {
    y_obs_train[auris_train_idx] ~ lognormal(log(basal + to_vector(mu[time_train_idx[auris_train_idx], 1])/delta), sigma);
    y_obs_train[blank_train_idx] ~ lognormal(log(basal), sigma);
  }
}

generated quantities {
  array[T_test] vector[1] mu_pred;
  real y_tot[N_train+N_test];
  real y_rep[N_train];
  real y_pred[N_test];
  vector[N_log_lik] log_lik;  // log likelihood

  y_rep[auris_train_idx] = lognormal_rng(log(basal + to_vector(mu[time_train_idx[auris_train_idx], 1])/delta), sigma);
  for (i in blank_train_idx) y_rep[i] = lognormal_rng(log(basal), sigma);
  array[T_test-1] vector[1] mu_hat_pred = ode_rk45(edwards, to_vector({y0}), test_time[1], test_time[2:T_test], beta, L);  // solve edwards ODE
  mu_pred[1, 1] = y0;
  mu_pred[2:T_test, 1] = mu_hat_pred[, 1];
  y_pred[auris_test_idx] = lognormal_rng(log(basal + to_vector(mu_pred[time_test_idx[auris_test_idx], 1])/delta), sigma);
  for (i in blank_test_idx) y_pred[i] = lognormal_rng(log(basal), sigma);

  // calculate log likelihood
  if (include_likelihood) {
    for (i in blank_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal), sigma);
    for (i in auris_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal + mu_pred[time_test_idx[i], 1]/delta), sigma);
  }
  y_tot[y_train_idx] = y_rep;
  y_tot[y_test_idx] = y_pred;
}
