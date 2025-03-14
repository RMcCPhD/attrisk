
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
  real a; // shape
  vector[I] b; // rate
  vector[I] p; // probability of attrition
  
  a = a_raw == 0 ? 0.0001 : a_raw; // ensure 'a' is not exactly zero
  for (i in 1:I) {
    b[i] = exp(mu1 + var1[i] * beta + u[i]);
    p[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
  }
}

model {
  // priors
  mu1 ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  a_raw ~ normal(0, 1);
  u ~ normal(0, 0.5);
  
  // likelihood
  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new; // predicted count
  vector[I] p_new; // predicted probability of attrition
  
  for (i in 1:I) {
    p_new[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }
  real rate = exp(mu1); // transformed rate
}
