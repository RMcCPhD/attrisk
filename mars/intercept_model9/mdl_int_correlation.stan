functions {
  // Original Gompertz function
  vector gompertz_prob(vector b, vector a, vector dt) {
    int I = rows(b);
    vector[I] p;
    for (i in 1:I) {
      if (fabs(a[i]) < 1e-4) {
        p[i] = 1 - exp(-b[i] * dt[i]);
      } else {
        p[i] = 1 - exp(-b[i] / a[i] * (exp(a[i] * dt[i]) - 1));
      }
    }
    return p;
  }
  
  // Alternative: parameterize directly in terms of cumulative probability
  vector logit_prob(vector logit_p_base, vector shape_effect, vector dt) {
    int I = rows(logit_p_base);
    vector[I] p;
    for (i in 1:I) {
      // Use a flexible shape that can capture Gompertz-like behavior
      p[i] = inv_logit(logit_p_base[i] + shape_effect[i] * log(dt[i] + 1));
    }
    return p;
  }
  
  // Alternative: parameterize in terms of initial rate and acceleration
  vector rate_accel_prob(vector log_rate_init, vector acceleration, vector dt) {
    int I = rows(log_rate_init);
    vector[I] p;
    for (i in 1:I) {
      // Rate increases exponentially with time: rate(t) = rate_init * exp(accel * t)
      real rate_t = exp(log_rate_init[i] + acceleration[i] * dt[i]);
      p[i] = 1 - exp(-rate_t * dt[i]);
    }
    return p;
  }
}

data {
  int<lower=0> I;                // number of arms
  array[I] int<lower=0> r;       // observed attrition
  array[I] int<lower=0> N;       // total enrolled
  vector[I] dt;                  // follow-up time in years
  int<lower=0> K;                // number of covariates
  matrix[I, K] X;                // design matrix of covariates
}

parameters {
  real mu_log_rate;              // mean log-rate
  real<lower=0> sigma_log_rate;  // SD of log-rate
  vector[I] z_log_rate;          // standard normal for log-rate
  
  // Shape as function of rate and covariates
  real alpha;                    // intercept for shape
  real beta_rate;                // coefficient for rate effect on shape
  vector[K] beta_cov;            // coefficients for covariates
  real<lower=0> sigma_shape;     // residual SD for shape
  vector[I] z_shape;             // standard normal for shape residuals
}

transformed parameters {
  vector[I] log_rate;            // log-rate for each arm
  vector[I] rate;                // rate for each arm
  vector[I] shape;               // shape for each arm
  vector[I] p;                   // attrition probability

  // Generate log-rates
  log_rate = mu_log_rate + sigma_log_rate * z_log_rate;
  rate = exp(log_rate);
  
  // Generate shapes as function of rate and covariates
  for (i in 1:I) {
    shape[i] = alpha + beta_rate * log_rate[i] + X[i] * beta_cov + sigma_shape * z_shape[i];
  }
  
  p = gompertz_prob(rate, shape, dt);
}

model {
  // Priors for rate parameters
  mu_log_rate ~ normal(-1.0, 1.0);
  sigma_log_rate ~ normal(0.5, 0.3);
  
  // Priors for shape parameters
  alpha ~ normal(-1.0, 1.0);           // shape intercept
  beta_rate ~ normal(0, 0.5);          // rate effect on shape
  beta_cov ~ normal(0, 0.5);           // covariate effects
  sigma_shape ~ normal(0.5, 0.3);      // shape residual SD
  
  // Latent variables
  z_log_rate ~ normal(0, 1);
  z_shape ~ normal(0, 1);
  
  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;
  
  // Parameter summaries
  real mean_log_rate = mean(log_rate);
  real mean_shape = mean(shape);
  real mean_rate = mean(rate);
  
  // Model fit statistics
  real deviance = -2 * sum(log_lik);
  
  // Parameter ranges
  real min_rate = min(rate);
  real max_rate = max(rate);
  real min_shape = min(shape);
  real max_shape = max(shape);
  
  // Rate-shape relationship
  real corr_rate_shape = corr(log_rate, shape);
  
  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }
}