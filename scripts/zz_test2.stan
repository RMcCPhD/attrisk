
data {
  int<lower=1> I; // Number of trials
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial
  vector[I] var1; // Covariate
}

parameters {
  real mu1; // Log-baseline hazard rate
  vector[I] u; // Random effect for rate
  real beta; // Coefficient
  real a_raw; // Raw shape
  real alpha; // Shape intercept
}

transformed parameters {
  vector[I] a; // Shape per trial
  vector[I] b; // Rate per trial
  vector[I] p; // Attrition probability per trial

  for (i in 1:I) {
    a[i] = fmin(fmax(exp(a_raw + alpha + beta * var1[i]), 1e-3), 5); // Place covariate on shape
    b[i] = exp(mu1 + u[i]); // Combine baseline hazard rate with random effect
    
    real exp_term = fmin(exp(a[i] * dt[i]) - 1, 1000);
    real rate_ratio = fmin(b[i] / a[i], 100);
    
    p[i] = 1 - exp(-rate_ratio * exp_term); // Get p as Gompertz CDF
  }
}

model {
  beta ~ normal(0, 1);
  mu1 ~ normal(0, 1);
  u ~ normal(0, 0.2);
  alpha ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  
  // likelihood
  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new; // Predicted count
  vector[I] p_new; // Predicted probability of attrition
  vector[I] rate; // Baseline hazard rate per trial
  
  for (i in 1:I) {
    real exp_term = exp(a[i] * dt[i]) - 1;
    real rate_ratio = b[i] / a[i];
    
    p_new[i] = 1 - exp(-rate_ratio * exp_term);
    r_new[i] = binomial_rng(N[i], p_new[i]);
    rate[i] = exp(mu1 + u[i]);
  }
}
