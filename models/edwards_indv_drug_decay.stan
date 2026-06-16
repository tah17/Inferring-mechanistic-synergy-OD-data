//
// Latent Edwards function modelled with multiplicative noise and individual drug action added with added decay of dead fungi
//

functions {
  vector edwards_drug(real t,
                      vector y,
                      real beta,
                      real L,
                      real delta_d_star,
                      real kappa_star,
                      real gamma,
                      real epsilon,
                      real max_time) {
    vector[2] dydt;
    dydt[1] = max_time*(beta/(1+gamma))*y[1]*exp(-((y[1]*(1+epsilon))/L)) - kappa_star*y[1];   // ODE for viable fungi scaled by max_time
    dydt[2] = kappa_star*y[1] - delta_d_star*y[2];   // ODE for dead fungi scaled by max_time and added decay rate
    return dydt;
  }
}

data {
  int<lower=0> D;  // number of drugs + 1 (+1 includes no drug option)
  int<lower=1> D_c;  // number of drug conditions 
  int<lower=0, upper=1> X[D_c, D];  // binary matrix of drug conditions x drugs

  int<lower=1> N_train;  // number of training obs
  int<lower=0> T_train;  // number of unique train time points
  array[T_train] real train_time;  //  train time points
  int<lower=1, upper=T_train> time_train_idx[N_train];  //  train time indices
  real max_train_time;  //  max time 
  int<lower=0> N_b_train;  // number of blanks in training
  int<lower=0, upper=N_train> blank_train_idx[N_b_train];  // blank indices
  int<lower=0> N_auris_train;  // number of non-blanks in training
  int<lower=0, upper=N_train> auris_train_idx[N_auris_train];  // non-blank indices
  vector[N_train] y_obs_train; // input data of OD values
  int<lower=1, upper=D_c> drug_cond_train_idx[N_train];  // indices of drug conditions in training
  
  int<lower=0> N_test;  // number of predictions
  int<lower=0> T_test;  // number of unique test time points
  array[T_test] real test_time;  // test time points
  int<lower=1, upper=T_test> time_test_idx[N_test];  // test time indices
  int<lower=0> N_b_test;  // number of blanks in testing
  int<lower=0, upper=N_test> blank_test_idx[N_b_test];  // blank indices
  int<lower=0> N_auris_test;  // number of non-blanks in testing
  int<lower=0, upper=N_test> auris_test_idx[N_auris_test];  // non-blank testing indices
  vector[N_test] y_obs_test; // input data of OD test values
  int<lower=1, upper=D_c> drug_cond_test_idx[N_test];  // indices of drug conditions in testing

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
  real<lower=0> y0;  // initial value for auris in wells
  real L_tilde;  // log10 growth impedance constant
  real<lower=log10(IC)> delta_tilde;  // log10 linear transform parameter (scale)
  real<lower=0> beta;  // growth rate
  real<lower=0> basal;  // linear transform parameter (offset)
  real<lower=0> sigma;  // scale of observed noise
  real<lower=0> delta_d_star;  // decay rate of dead fungi
  // params for each D-1 drugs
  vector<lower=0>[D-1] kappa_star;  // killing
  vector<lower=0>[D-1] gamma;  // inhibition of growth rate
  vector<lower=0>[D-1] epsilon;   // enhancement of growth impedance
}

transformed parameters {
  real<lower=0> L = pow(10, L_tilde);  // growth impedance const  
  real<lower=IC> delta = pow(10, delta_tilde);  // linear transform parameter (scale)
  array[T_train-1, D_c] vector[2] mu_hat;
  array[T_train, D_c] vector[1] mu;
  for (i in 1:D_c) {
    row_vector[D] X_i = to_row_vector(X[i, ]);
    mu_hat[, i, ] = ode_rk45(edwards_drug, to_vector({y0, 0.0}), t_star_train[1], t_star_train[2:T_train], beta, L, delta_d_star, X_i*append_row(0.0, kappa_star), X_i*append_row(0.0, gamma), X_i*append_row(0.0, epsilon), max_train_time);  // solve ODE
    mu[1, i, 1] = y0;
    mu[2:T_train, i, 1] = to_array_1d(to_vector(mu_hat[, i, 1]) + to_vector(mu_hat[, i, 2]));   // viable and dead fungi added
  }
}

model {
  y0 ~ lognormal(log(IC), 1);
  sigma ~ normal(0, 0.5);
  L_tilde ~ normal(log10(IC), 2);
  beta ~ std_normal();
  delta_d_star ~ std_normal();
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
  real<lower=0> delta_d = delta_d_star/max_train_time; // decay rate in hrs^(-1)
  vector<lower=0>[D-1] kappa = kappa_star/max_train_time; // killing in hrs^(-1)
  array[T_test-1, D_c] vector[2] mu_hat_pred;
  array[T_test, D_c] vector[1] mu_pred;
  real y_tot[N_train+N_test];
  real y_rep[N_train];
  real y_pred[N_test];
  vector[N_log_lik] log_lik;  // log likelihood
  
  for (i in auris_train_idx) y_rep[i] = lognormal_rng(log(basal + mu[time_train_idx[i], drug_cond_train_idx[i], 1]/delta), sigma);
  for (i in blank_train_idx) y_rep[i] = lognormal_rng(log(basal), sigma);
  for (i in 1:D_c) {  // predictions for CV stratified by rep
    row_vector[D] X_i = to_row_vector(X[i, ]);
    mu_hat_pred[, i, ] = ode_rk45(edwards_drug, to_vector({y0, 0.0}), t_star_test[1], t_star_test[2:T_test], beta, L, delta_d_star, X_i*append_row(0.0, kappa_star), X_i*append_row(0.0, gamma), X_i*append_row(0.0, epsilon), max_train_time);  // solve ODE
    mu_pred[1, i, 1] = y0;
    mu_pred[2:T_test, i, 1] = to_array_1d(to_vector(mu_hat_pred[, i, 1]) + to_vector(mu_hat_pred[, i, 2]));
  }
  for (i in auris_test_idx) y_pred[i] = lognormal_rng(log(basal + mu_pred[time_test_idx[i], drug_cond_test_idx[i], 1]/delta), sigma);
  for (i in blank_test_idx) y_pred[i] = lognormal_rng(log(basal), sigma);
  // calculate log likelihood
  if (include_likelihood) {
    for (i in blank_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal), sigma);
    for (i in auris_test_idx)
      log_lik[i] = lognormal_lpdf(y_obs_test[i] | log(basal + mu_pred[time_test_idx[i],  drug_cond_test_idx[i], 1]/delta), sigma);
  }
  y_tot[y_train_idx] = y_rep;
  y_tot[y_test_idx] = y_pred;
}
