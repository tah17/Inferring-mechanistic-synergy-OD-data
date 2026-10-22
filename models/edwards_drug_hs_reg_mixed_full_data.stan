//
// Latent Edwards function modelled with multiplicative noise and added drug action + mixed effects on growth rate, L and drug action params
//

functions {
  vector edwards_drug(real t,
                      vector y,
                      real beta,
                      real L,
                      real kappa_star,
                      real gamma,
                      real epsilon,
                      real max_time) {
    vector[2] dydt;
    dydt[1] = max_time*(beta/(1+gamma))*y[1]*exp(-((y[1]*(1+epsilon))/L)) - kappa_star*y[1];  // ODE for viable fungi scaled by max_time
    dydt[2] = kappa_star*y[1];   // ODE for dead fungi scaled by max_time
    return dydt;
  }
}

data {
  int<lower=0> D;  // number of drugs + 1 (+1 includes no drug option)
  int<lower=0> D_i;  // number of drug interactions 
  int<lower=1> D_c;  // number of drug conditions 
  int<lower=0, upper=1> X[D_c, D+D_i];  // binary matrix of drug conditions x drugs

  int<lower=1> N;  // number of obs
  int<lower=0> T;  // number of unique time points
  array[T] real time;  // time points
  int<lower=1, upper=T> time_idx[N];  // time indices
  real max_time;  //  max time 
  int<lower=0> R;  // number reps
  int<lower=0> R_auris;  // number auris reps
  int<lower=1, upper=R> rep_idx[N];  // rep indices 
  int<lower=0> N_b;  // number of blanks
  int<lower=0, upper=N> blank_idx[N_b];  // blank indices
  int<lower=0> N_auris;  // number of non-blanks 
  int<lower=0, upper=N> auris_idx[N_auris];   // non-blank indices
  vector[N] y_obs; // input data of OD values
  int<lower=1, upper=R> auris_in_reps[R_auris];  // indices of auris reps 
  int<lower=1, upper=D_c> drug_cond_for_rep[R_auris];  // drug cond indices for auris reps

  int<lower=0, upper=1> include_likelihood;  // if the likelihood is included in the model or not, e.g. set to 0 during a prior predictive check
}

transformed data {
  int<lower=0> N_log_lik = N*include_likelihood;  // sets the data size to 0 if the likelihood is not in the model
  real<lower=0> IC = (2.5)*1e2;  // initial inoculum
  real<lower=0> slab_scale = 2;
  real<lower=0> scale_global = 1;
  array[T] real t_star = to_array_1d(to_vector(time)/max_time);  // transformed time points
}

parameters {
  real<lower=0> y0;   // initial value for auris in wells
  vector[R_auris] L_tilde_eta;  // non cenetred log10 growth impedance constant
  real log10_mu_L_tilde;  // mean log10 growth impedance constant
  real<lower=0> sigma_L_tilde;  // sd of log10 growth impedance constant
  real<lower=log10(IC)> delta_tilde;  // log10 linear transform parameter (scale)
  real<lower=0> mu_beta;  // mean growth rate
  real<lower=0> sigma_log_beta;  // sd of logged growth rate
  vector[R_auris] log_beta_eta;  // non centred growth rate per rep
  real<lower=0> basal;   // offset per rep
  real<lower=0> sigma;  // scale of observed noise
  // population params for each D-1 drugs and D_i drug combinations
  vector<lower=0>[D+D_i-1] mu_kappa_star; // killing
  vector<lower=0>[D+D_i-1] mu_gamma;   // inhibition of growth rate
  vector<lower=0>[D+D_i-1] mu_epsilon;  // enhancement of growth impedance
  real<lower=0> sigma_log_drug;  // scale of drug noise
  // non_centered params for each D-1 drugs and D_i drug combinations
  vector[R_auris-3] kappa_star_eta; // killing ( - 3 RPMI reps, this works due to the order of the data where the RPMI reps are last)
  vector[R_auris-3] gamma_eta;   // inhibition of growth rate
  vector[R_auris-3] epsilon_eta;  // enhancement of growth impedance
  // regularised horseshoe params (https://doi.org/10.1214/17-EJS1337SI)
  real<lower=0> caux;
  vector<lower=0>[(D+D_i-1)*3] lambda;  // local regularisation for each drug
  real<lower=0> tau;  // global regularisation
}

transformed parameters {
  real<lower=0> c = slab_scale*sqrt(caux);
  vector<lower=0>[(D+D_i-1)*3] lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda)));
  vector<lower=0>[R_auris-3] kappa_star; // killing ( - 3 RPMI reps, this works due to the order of the data where the RPMI reps are last)
  vector<lower=0>[R_auris-3] gamma;   // inhibition of growth rate
  vector<lower=0>[R_auris-3] epsilon;  // enhancement of growth impedance
  vector[R_auris] log_beta = log(mu_beta) + sigma_log_beta*log_beta_eta;   // log growth rate per rep
  vector<lower=0>[R_auris] beta = exp(log_beta);   // growth rate per rep
  vector[R_auris] L_tilde = log10_mu_L_tilde + sigma_L_tilde*L_tilde_eta;   // log10 growth impedance constant
  vector<lower=0>[R_auris] L = pow(10, L_tilde);  // growth impedance const
  real<lower=IC> delta = pow(10, delta_tilde);  // linear transform parameter (scale)
  array[T-1, R] vector[2] mu_hat = rep_array(rep_vector(0, 2), T-1, R);
  array[T, R] vector[1] mu = rep_array(rep_vector(0, 1), T, R);
  for (i in 1:R_auris) {
    row_vector[D+D_i] X_i = to_row_vector(X[drug_cond_for_rep[i], ]);
    if (drug_cond_for_rep[i] == 2) {
      mu_hat[, auris_in_reps[i], ] = ode_rk45(edwards_drug, to_vector({y0, 0.0}), t_star[1], t_star[2:T], beta[i], L[i], 0, 0, 0, max_time);  // solves no drug Edwards ODE
    } else {
      # non-centered hierarchy
      real log_kappa_star = log(X_i*append_row(0.0, mu_kappa_star)) + sigma_log_drug*kappa_star_eta[i];
      real log_gamma = log(X_i*append_row(0.0, mu_gamma)) + sigma_log_drug*gamma_eta[i];
      real log_epsilon = log(X_i*append_row(0.0, mu_epsilon)) + sigma_log_drug*epsilon_eta[i];
      
      kappa_star[i] = exp(log_kappa_star);
      gamma[i] = exp(log_gamma);
      epsilon[i] = exp(log_epsilon);

      mu_hat[, auris_in_reps[i], ] = ode_rk45(edwards_drug, to_vector({y0, 0.0}), t_star[1], t_star[2:T], beta[i], L[i], kappa_star[i], gamma[i], epsilon[i], max_time);  // solves Edwards ODE
    }
    mu[1, auris_in_reps[i], 1] = y0;
    mu[2:T, auris_in_reps[i], 1] = to_array_1d(to_vector(mu_hat[, auris_in_reps[i], 1]) + to_vector(mu_hat[, auris_in_reps[i], 2]));   // viable and dead fungi added
  }
}

model {
  y0 ~ lognormal(log(IC), 1);
  sigma ~ normal(0, 0.5);
  log10_mu_L_tilde ~ normal(log10(IC), 2);
  L_tilde_eta ~ std_normal();
  sigma_L_tilde ~ normal(0, 0.5);
  mu_beta ~ std_normal();
  sigma_log_beta ~ normal(0, 0.5);
  log_beta_eta ~ std_normal();
  delta_tilde ~ cauchy(log10(IC), 1);
  kappa_star_eta ~ std_normal();
  gamma_eta ~ std_normal();
  epsilon_eta ~ std_normal();
  sigma_log_drug ~ normal(0, 0.5);
  basal ~ lognormal(0, 1);
  // regularised horseshoe
  tau ~ cauchy(0, scale_global*sigma);
  caux ~ inv_gamma(4*0.5, 4*0.5); // default value in packages is usually slab_df=4
  lambda ~ cauchy(0, 1);
  for (i in 1:(D+D_i-1)) {
    mu_kappa_star[i] ~ normal(0, lambda_tilde[i]*tau);
    mu_gamma[i] ~ normal(0, lambda_tilde[i+(D+D_i-1)]*tau);
    mu_epsilon[i] ~ normal(0, lambda_tilde[i+(2*(D+D_i-1))]*tau);
  }

  if (include_likelihood) {
    for (i in auris_idx) y_obs[i] ~ lognormal(log(basal + mu[time_idx[i], rep_idx[i], 1]/delta), sigma);
    y_obs[blank_idx] ~ lognormal(log(basal), sigma);
  }
}

generated quantities {
  vector<lower=0>[R_auris-3] kappa = kappa_star/max_time; // killing rate in hrs^(-1)
  vector<lower=0>[D+D_i-1] mu_kappa = mu_kappa_star/max_time; // killing rate in hrs^(-1)
  real y_tot[N];
  vector[N_log_lik] log_lik;  // log likelihood
  
  for (i in auris_idx) y_tot[i] = lognormal_rng(log(basal + mu[time_idx[i], rep_idx[i], 1]/delta), sigma); 
  for (i in blank_idx) y_tot[i] = lognormal_rng(log(basal), sigma);
  // calculate log likelihood
  if (include_likelihood) {
    for (i in blank_idx)
      log_lik[i] = lognormal_lpdf(y_obs[i] | log(basal), sigma);
    for (i in auris_idx)
      log_lik[i] = lognormal_lpdf(y_obs[i] | log(basal + mu[time_idx[i], rep_idx[i], 1]/delta), sigma);
  }
}
