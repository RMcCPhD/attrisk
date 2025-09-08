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
  vector[I] gamma_raw;                  // arm-specific calibration
  real<lower=0> gamma_sigma;            // shrinkage scale
}

transformed parameters {
  // Non-centred and sum-to-zero gamma
  vector[I] gamma = gamma_sigma * gamma_raw;
  vector[I] gamma_i = gamma - mean(gamma);
  
  vector[I] log_b_pre = mu1_draw + X * beta_draw;             // pre-CITL log-rate
  vector[I] log_b_cal = log_b_pre + gamma_i;                  // calibrated log-rate
  vector[I] b_cal = exp(log_b_cal);                           // rate
  vector[I] a_cal = log_b_cal * beta_rate + lp + error_draw;  // shape
}

model {
  // Normal priors for CITL gamma
  gamma_raw ~ normal(0, 1);
  gamma_sigma ~ normal(0, 1);
  
  // Compute arm-level average risk from participant HRs
  for (i in 1:I) {
    int s = (i - 1) * J + 1;
    vector[J] hr_i = segment(hr_all, s, J);
    
    real a_i = a_cal[i];
    real b_i = b_cal[i];
    real t_i = dt[i];
    
    // Safely compute Gompertz constant used to compute risk
    real C_i = (fabs(a_i) < 1e-4) ? (b_i * t_i) : (b_i / a_i) * expm1(a_i * t_i);
    
    // Average participant-level risk per arm
    real p_i = clamp_p(1 - mean(exp(-hr_i * C_i)));
    
    // Pass average risk to likelihood
    target += binomial_lpmf(r[i] | N[i], p_i);
  }
}

generated quantities {
  vector[I] delta = log_b_cal - log_b_pre;
  vector[I] haz_mult = exp(delta);
}

