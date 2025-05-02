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
  real mu1_mean; // Global baseline log-hazard
  vector[I] mu1_raw; // Raw baseline log-hazard
  real<lower=0> sigma_mu1; // Baseline log-hazard variation

  real a; // Fixed shape parameter across all trials

  real beta_arm;
  vector[K] z;
  real<lower=0> tau; // Global shrinkage factor
  vector<lower=0>[K] lambda; // Local shrinkage factor
}

transformed parameters {
  vector[I] mu1 = mu1_mean + mu1_raw * sigma_mu1;
  vector[K] beta = z .* tau .* lambda;

  vector[I] b;
  vector[I] p;

  for (i in 1:I) {
    b[i] = exp(mu1[i] + arm[i] * beta_arm + dot_product(X[i], beta));
    real exponent = a * dt[i];
    real time_part = (fabs(a) < 1e-6)
      ? dt[i]
      : expm1(exponent) / a;
    p[i] = 1 - exp(-b[i] * time_part);
  }
}

model {
  mu1_mean ~ normal(0, 1);
  mu1_raw ~ normal(0, 1);
  sigma_mu1 ~ normal(0, 1);

  a ~ normal(0, 1); // Fixed shape parameter

  beta_arm ~ normal(0, 1);
  z ~ normal(0, 1);
  lambda ~ normal(0, 1);
  tau ~ normal(0, 1);

  r ~ binomial(N, p);
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;
  vector[I] rate = exp(mu1);
  real mean_rate = mean(rate);
  real sd_rate = sd(rate);
  real shrinkage_effective_params = sum(lambda ./ sqrt(1 + lambda.^2));

  for (i in 1:I) {
    real
