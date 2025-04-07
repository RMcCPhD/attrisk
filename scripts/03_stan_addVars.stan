
data {
  int<lower=0> I; // Number of trials
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial
  
  // Covariates
  vector[I] arm; // Treatment arm (0 = placebo, 1 = active)
  vector[I] var1; // Number of arms (integer)
  vector[I] var2; // Sample size (scaled)
  vector[I] var3; // Number of sites (scaled)
}

parameters {
  real mu1; // intercept
  real a_raw; // raw parameter for 'a'
  
  // Betas
  real beta_arm;
  real beta1;
  real beta2;
  real beta3;
}

transformed parameters {
  real a; // Shape
  vector[I] b; // Rate
  vector[I] p; // Probability of attrition
  
  a = a_raw == 0 ? 0.0001 : a_raw; // ensure 'a' is not exactly zero
    
  for (i in 1:I) {
    b[i] = exp(
      mu1 + beta_arm * arm[i] + beta1 * var1[i] + beta2 * var2[i] +
      beta3 * var3[i]
    );
    p[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
  }
}

model {
  // priors
  mu1 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  beta_arm ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  
  // likelihood
  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;
  
  for (i in 1:I) {
    p_new[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }
  real rate = exp(mu1); // Baseline rate
}
