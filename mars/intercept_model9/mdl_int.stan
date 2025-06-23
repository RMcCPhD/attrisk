functions {
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
}

data {
  int<lower=0> I;                // number of arms
  array[I] int<lower=0> r;       // observed attrition
  array[I] int<lower=0> N;       // total enrolled
  vector[I] dt;                  // follow-up time
}

parameters {
  vector[I] z;                   // latent normal variable
}

transformed parameters {
  vector[I] mu1;                 // log-rate
  vector[I] b;                   // rate
  vector[I] a;                   // shape
  vector[I] p;                   // attrition probability
  
  real mu1_mean = -5.237;
  real mu1_sd = 0.071;
  
  real a_intercept = -0.118;
  real a_slope = -0.02;

  for (i in 1:I) {
    mu1[i] = mu1_mean + mu1_sd * z[i]; 
    b[i] = exp(mu1[i]);
    a[i] = a_intercept + a_slope * (mu1[i] - mu1_mean);
  }

  p = gompertz_prob(b, a, dt);
}

model {
  z ~ normal(0, 1);               // standard normal latent
  r ~ binomial(N, p);             // likelihood
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

  real mean_log_rate = mean(mu1);
  real mean_shape = mean(a);
}
