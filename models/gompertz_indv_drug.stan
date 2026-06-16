//
// Latent gompertz function modelled with multiplicative noise with individual drug action added
//

functions {
  vector gompertz_drug(real t,
                       vector y,
                       real y0,
                       real beta_gamma_max_time,
                       real K_epsilon,
                       real kappa_star) {
    vector[1] dydt;
    dydt[1] = kappa_star*y0*exp((log(K_epsilon/y0)-(kappa_star/beta_gamma_max_time))*(1-exp(-beta_gamma_max_time*t)));  // ODE of dead fungi scaled by max time
    return dydt;
  }
}

data {
  int<lower=0> D;  // number of drugs (with no drug added)
  int<lower=1> D_c;  // number of drug conditions 
  int<lower=0, upper=1> X[D_c, D];  // matrix of drug condition x drug indicators
  
  int<lower=1> N_train;  // number of obs
  int<lower=0> T_train;  // number of unique train time points
  array[T_train] real train_time;  // train time
  int<lower=1, upper=T_train> time_train_idx[N_train];  // train time idxs
  real max_train_time;  //  max time 
  int<lower=0> N_b_train;  // number of blanks in training
  int<lower=0, upper=N_train> blank_train_idx[N_b_train];  // blank indices
  int<lower=0> N_auris_train;  // number of data points for wells with c. auris in training
  int<lower=0, upper=N_train> auris_train_idx[N_auris_train];  // auris indices in training
  vector[N_train] y_obs_train; // input data of OD values
  int<lower=1, upper=D_c> drug_cond_train_idx[N_train];  // indices of drug conditions in training
  
  int<lower=0> N_test;  // number of predictions
  int<lower=0> T_test;  // number of unique test time points
  array[T_test] real test_time;  // test times
  int<lower=1, upper=T_test> time_test_idx[N_test];  // test time idxs
  int<lower=0> N_b_test;  // number of blanks in testing
  int<lower=0, upper=N_test> blank_test_idx[N_b_test];  // blank indices
  int<lower=0> N_auris_test;  // number of data points for wells with c. auris in testing
  int<lower=0, upper=N_test> auris_test_idx[N_auris_test];  // auris idxs in testing
  vector[N_test] y_obs_test; // input data of OD test values
  int<lower=1, upper=D_c> drug_cond_test_idx[N_test]; // indices of drug conditions in testing

  int<lower=1, upper=N_train+N_test> y_train_idx[N_train];  // training indices
  int<lower=1, upper=N_train+N_test> y_test_idx[N_test];  // testing indices
  
  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N_log_lik = N_test*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  real<lower=0> IC = (2.5)*1e2;  // initial inoculum
  array[T_train] real t_star_train = to_array_1d(to_vector(train_time)/max_train_time);  // transformed train time points
  array[T_test] real t_star_test = to_array_1d(to_vector(test_time)/max_train_time);  // transformed test time points
}

parameters {
  real<lower=0> y0;  // initial inoculum parameter
  real<lower=log10(IC)> K_tilde;  // log10 carrying capacity 
  real<lower=log10(IC)> delta_tilde;  // log10 linear transform parameter (scale)
  real<lower=0> beta;  // initial growth rate
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> sigma;  // scale of observed noise
  // params for each D-1 drugs
  vector<lower=0>[D-1] kappa_star;  // killing
  vector<lower=0>[D-1] gamma;  // inhibition of growth rate
  vector<lower=0>[D-1] epsilon;   // enhancement of carrying capacity reduction
}

transformed parameters {
  matrix[T_train-1, D_c] f_v;  // viable fungi
  real<lower=IC> K = pow(10, K_tilde);  // carrying capacity
  real<lower=IC> delta = pow(10, delta_tilde);  // linear transform parameter (scale)
  array[T_train-1, D_c] vector[1] mu_hat;
  array[T_train, D_c] vector[1] mu;
  for (i in 1:D_c) {
    row_vector[D] X_i = to_row_vector(X[i, ]);
    real beta_gamma_max_time = (beta*max_train_time)/(1+X_i*append_row(0.0, gamma));
    real K_epsilon = K/(1+X_i*append_row(0.0, epsilon));
    real kappa_tilde = X_i*append_row(0.0, kappa_star);
    f_v[, i] = y0*exp((log(K_epsilon/y0)-(kappa_tilde/beta_gamma_max_time))*(1-exp(-beta_gamma_max_time*to_vector(t_star_train[2:T_train]))));  // gompertz scaled by max time for viable fungi
    mu_hat[, i, ] = ode_rk45(gompertz_drug, to_vector({0.0}), t_star_train[1], t_star_train[2:T_train], y0, beta_gamma_max_time, K_epsilon, kappa_tilde);  // solve ODE for dead fungi
    mu[1, i, 1] = 0;
    f_v[1, i] = y0;
    mu[2:T_train, i, 1] = to_array_1d(to_vector(f_v[, i]) + to_vector(mu_hat[, i, 1]));  // viable and dead fungi added
  }
}

model {
  y0 ~ lognormal(log(IC), 1);
  sigma ~ normal(0, 0.5);
  K_tilde ~ normal(9, 2);
  beta ~ std_normal();
  kappa_star ~ std_normal();
  gamma ~ std_normal();
  epsilon ~ std_normal();
  basal ~ lognormal(0, 1);
  delta_tilde ~ cauchy(log10(IC), 1);

  if (include_likelihood) {
    for (i in auris_train_idx) y_obs_train[i] ~ lognormal(log(basal + mu[time_train_idx[i], drug_cond_train_idx[i], 1]/delta), sigma);
    y_obs_train[blank_train_idx] ~ lognormal(log(basal), sigma);
  }
}

generated quantities {
  vector<lower=0>[D-1] kappa = kappa_star/max_train_time; // untransformed killing rate
  matrix[T_test-1, D_c] f_v_test;  // predicted viable fungi
  array[T_test-1, D_c] vector[1] mu_hat_pred;
  array[T_test, D_c] vector[1] mu_pred;
  real y_tot[N_train+N_test];
  real y_rep[N_train];
  real y_pred[N_test];
  vector[N_log_lik] log_lik;  // log likelihood

  for (i in auris_train_idx) y_rep[i] = lognormal_rng(log(basal + mu[time_train_idx[i], drug_cond_train_idx[i], 1]/delta), sigma);
  for (i in blank_train_idx) y_rep[i] = lognormal_rng(log(basal), sigma);
  for (i in 1:D_c) {  // predictions for CV stratified by rep
    row_vector[D] X_i = to_row_vector(X[i, ]);
    real beta_gamma_max_time = (beta*max_train_time)/(1+X_i*append_row(0.0, gamma));
    real K_epsilon = K/(1+X_i*append_row(0.0, epsilon));
    real kappa_tilde = X_i*append_row(0.0, kappa_star);
    f_v_test[, i] = y0*exp((log(K_epsilon/y0)-(kappa_tilde/beta_gamma_max_time))*(1-exp(-beta_gamma_max_time*to_vector(t_star_test[2:T_test]))));  // gompertz transformed by max time
    mu_hat_pred[, i, ] = ode_rk45(gompertz_drug, to_vector({0.0}), t_star_test[1], t_star_test[2:T_test], y0, beta_gamma_max_time, K_epsilon, kappa_tilde);  // ODE for dead fungi
    mu_pred[1, i, 1] = 0;
    f_v_test[1, i] = y0;
    mu_pred[2:T_test, i, 1] = to_array_1d(to_vector(f_v_test[, i]) + to_vector(mu_hat_pred[, i, 1]));  // predictions for viable and dead fungi added
  }
  for (i in auris_test_idx) y_pred[i] = lognormal_rng(log(basal + mu_pred[time_test_idx[i], drug_cond_test_idx[i], 1]/delta), sigma);
  for (i in blank_test_idx) y_pred[i] = lognormal_rng(log(basal), sigma);
  // calculate log likelihood
  if (include_likelihood) {
    for (i in blank_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal), sigma);
    for (i in auris_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal + mu_pred[time_test_idx[i], drug_cond_test_idx[i], 1]/delta), sigma);
  }
  y_tot[y_train_idx] = y_rep;
  y_tot[y_test_idx] = y_pred;
}
