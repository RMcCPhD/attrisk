// Clamp p to prevent exact 0 or 1s, allows STAN to sample in stable regions
functions {
  real clamp_p(real x) {
    return fmin(fmax(x, 1e-12), 1 - 1e-12);
  }
}

data {
  int<lower=1> I;                       // number arms
  int<lower=1> K;                       // number of covariates
  int<lower=1> J;                       // number of sobol participants
  array[I] int<lower=0> r;              // observed attrition counts
  array[I] int<lower=1> N;              // number of arm participants
  vector[I] dt;                         // total follow-up time in years
  
  // Pre-CITL draws
  vector[I] mu1_draw;                   // intercept of log-rate
  vector[K] beta_draw;                  // covariate betas
  matrix[I, K] X;                       // covariate matrix
  vector[I] lp;                         // linear predictor of shape
  real beta_rate;                       // rate beta from linear model of shape
  vector[I] error_draw;                 // residual errors passed to shape
  
  // CITL inputs
  int<lower=1> pK;                      // number of sobol covariates
  matrix[I * J, pK] Z;                  // stacked matrix of Sobol sequences
  vector[pK] pB;                        // pooled Cox log-HRs
}

transformed data {
  vector[I * J] hr_all = exp(Z * pB);   // participant-level HR
}

parameters {
  real gamma_mu;                        // global calibration
  vector[I] gamma_raw;                  // arm-specific calibration
  real<lower=0> gamma_sigma;            // shrinkage scale
}

transformed parameters {
  vector[I] gamma = gamma_mu + gamma_sigma * gamma_raw;
  vector[I] log_b_pre = mu1_draw + X * beta_draw;             // pre-CITL log-rate
  vector[I] log_b_cal = log_b_pre + gamma;                    // calibrated log-rate
  vector[I] b_cal = exp(log_b_cal);                           // rate
  vector[I] a_cal = log_b_cal * beta_rate + lp + error_draw;  // shape
  
  // Compute arm-level average participant risk
  vector[I] p_i;
  for (i in 1:I) {
    int s = (i - 1) * J + 1;
    vector[J] hr_i = segment(hr_all, s, J);
    real a_i = a_cal[i];
    real b_i = b_cal[i];
    real t_i = dt[i];
    
    // Safely compute Gompertz constant used to compute risk
    real C_i = (fabs(a_i) < 1e-4) ? (b_i * t_i) : (b_i / a_i) * expm1(a_i * t_i);
    p_i[i] = clamp_p(1 - mean(exp(-hr_i * C_i)));
  }
}

model {
  // Normal priors for CITL gamma
  gamma_mu ~ normal(0, 1);
  gamma_raw ~ normal(0, 1);
  gamma_sigma ~ normal(0, 1);
    
  // Pass average risk to likelihood
  target += binomial_lpmf(r | N, p_i);
    
}

generated quantities {
  vector[I] delta = log_b_cal - log_b_pre;
  vector[I] haz_mult = exp(delta);
  
  // PPC and log-likelihood
  array[I] int r_new;
  vector[I] log_lik;
  for (i in 1:I) {
    r_new[i]  = binomial_rng(N[i], p_i[i]);
    log_lik[i]= binomial_lpmf(r[i] | N[i], p_i[i]);
  }
}

