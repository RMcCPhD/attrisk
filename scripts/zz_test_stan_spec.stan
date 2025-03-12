
data {
  int<lower=0> N; // number of data points
  int<lower=0> r[N]; // observed counts
  int<lower=0> n[N]; // number of trials
  vector[N] dt; // time intervals
  vector[N] var1; // covar1 - numeric
}

parameters {
  real mu1; // intercept
  real<lower=0> tau_beta; // scale parameter for beta (Laplace prior, shrinkage)
  real beta_raw; // raw covar1 coefficient
  real a_raw; // raw parameter for 'a' (shape)
}

transformed parameters {
  real a; // shape
  vector[N] b; // b (rate)
  vector[N] p; // probability of attrition
  real beta = beta_raw * tau_beta; // apply shrinkage
  
  a = a_raw == 0 ? 0.0001 : a_raw; // ensure 'a' is not exactly zero
  for (i in 1:N) {
    b[i] = exp(mu1 + var1[i] * beta);
    p[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
  }
}

model {
  // priors
  mu1 ~ normal(0, 1);
  tau_beta ~ exponential(1); // Shrinkage prior
  beta_raw ~ normal(0, 1); // Normal prior for covar1
  a_raw ~ normal(0, 1); // Normal prior for shape
  
  // likelihood
  for (i in 1:N) {
    r[i] ~ binomial(n[i], p[i]);
  }
}

generated quantities {
  vector[N] r_new; // predicted count
  vector[N] p_new; // predicted probability of attrition
  
  for (i in 1:N) {
    p_new[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
    r_new[i] = binomial_rng(n[i], p_new[i]);
  }
  real rate = exp(mu1); // transformed rate
}
