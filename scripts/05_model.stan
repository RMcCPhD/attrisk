data {
  int<lower=0> I; // Number of trials
  int<lower=0> K; // Number of covariates
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial
  vector[I] arm; // Treatment arm (1 = placebo, 0 = active)
  matrix[I, K] X; // Covariates
}

parameters {
  real mu1; // Intercept
  real a_raw; // Raw shape parameter
  real beta_arm;
  
  vector[K] z;
  real<lower=0> tau; // Global shrinkage factor
  vector<lower=0>[K] lambda; // Local shrinkage factor
}

transformed parameters {
  real a = a_raw;
  
  vector[K] beta = z .* tau .* lambda;
  
  vector[I] b = exp(mu1 + arm * beta_arm + X * beta); // Rate
  vector[I] p;

  for (i in 1:I) {
    real exponent = a * dt[i];
    real time_part = (fabs(a) < 1e-6)
      ? dt[i]
      : expm1(exponent) / a;
    
    p[i] = 1 - exp(-b[i] * time_part); // Gompertz CDF
  }
}

model {
  mu1 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  beta_arm ~ normal(0, 1);
  z ~ normal(0, 1);
  
  lambda ~ student_t(3, 0, 1); // Local scale
  tau ~ student_t(3, 0, 1); // Global scale

  r ~ binomial(N, p);
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;
  
  real shrinkage_effective_params = sum(lambda ./ sqrt(1 + lambda.^2));

  for (i in 1:I) {
    real exponent = a * dt[i];
    real time_part = (fabs(a) < 1e-6)
      ? dt[i]
      : expm1(exponent) / a;
    p_new[i] = 1 - exp(-b[i] * time_part);
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }

  real rate = exp(mu1); // Baseline hazard
}
