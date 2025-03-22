
data {
  int<lower=1> I; // Number of trials
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial
  vector[I] var1; // Covariate
  vector[I] arm; // Treatment arm (0 = placebo, 1 = active)
}

parameters {
  real b_mu; // Log-baseline hazard rate
  vector[I] u; // Random intercept for rate
  real a_raw; // Shape intercept
  real beta; // Covariate coefficient
  real gamma; // Treatment arm effect
}

transformed parameters {
  vector[I] b; // Rate
  vector[I] a; // Shape
  vector[I] p; // Probability of attrition
  
  for (i in 1:I) {
    a[i] = a_raw + beta * var1[i] + gamma * arm[i];
    b[i] = exp(b_mu + beta * var1[i] + gamma * arm[i] + u[i]);
    p[i] = 1 - exp(-b[i]/a[i] * (exp(a[i] * dt[i]) - 1));
    
    // Check values
    // if (i <= 5) {
    //   print("a_raw = ", a_raw, ", var1[i] = ", var1[i], "beta = ", beta);
    //   print("a = ", a[i], ", b[i] = ", b[i]);
    //   print("p[i] = ", p[i]);
    // }
  }
}

model {
  b_mu ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  beta ~ normal(0, 1);
  u ~ normal(0, 1);
  gamma ~ normal(0, 1);
  
  // Likelihood
  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new; // Predicted count
  vector[I] p_new; // Predicted probability of attrition
  vector[I] rate; // Baseline hazard rate
  
  for (i in 1:I) {
    p_new[i] = 1 - exp(-b[i]/a[i] * (exp(a[i] * dt[i]) - 1));
    r_new[i] = binomial_rng(N[i], p_new[i]);
    rate[i] = exp(b_mu);
  }
}
