data {
  int<lower=0> I; // Number of trials
  int<lower=0> K; // Number of covariates
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in years
  vector[I] arm; // Treatment arm (1 = placebo, 0 = active)
  matrix[I, K] X; // Covariates
}

parameters {
  real mu1; // Intercept
  vector[I] a_raw; // Raw shape parameter
  real<lower=0> sigma_a; // SD of shape variation across trials
  
  real beta_arm;
  vector[K] z;
  real<lower=0> tau; // Global shrinkage factor
  vector<lower=0>[K] lambda; // Local shrinkage factor
}

transformed parameters {
  vector[I] a = a_raw * sigma_a; // Trial-specific shape
  vector[K] beta = z .* tau .* lambda; // Coefficient with shrinkage factors applied
  
  vector[I] b = exp(mu1 + arm * beta_arm + X * beta); // Rate
  vector[I] p;

  for (i in 1:I) {
    real exponent = a[i] * dt[i];
    real time_part = (fabs(a[i]) < 1e-6)
      ? dt[i]
      : expm1(exponent) / a[i];
    
    p[i] = 1 - exp(-b[i] * time_part); // Gompertz CDF
  }
}

model {
  mu1 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  sigma_a ~ normal(0,1);
  beta_arm ~ normal(0, 1);
  z ~ normal(0, 1);
  
  lambda ~ normal(0, 1); // Local scale
  tau ~ normal(0, 1); // Global scale

  r ~ binomial(N, p);
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;
  
  real shrinkage_effective_params = sum(lambda ./ sqrt(1 + lambda.^2));

  for (i in 1:I) {
    real exponent = a[i] * dt[i];
    real time_part = (fabs(a[i]) < 1e-6) // Safely approximate a if close to 0
      ? dt[i] // Reduce to exponential if 0
      : expm1(exponent) / a[i]; // Otherwise use full expression
    p_new[i] = 1 - exp(-b[i] * time_part);
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }

  real rate = exp(mu1); // Baseline hazard
}
