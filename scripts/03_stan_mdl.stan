
data {
  int<lower=0> I; // Number of trials
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial
  
  // Covariates
  vector[I] arm; // Treatment arm (0 = placebo, 1 = active)
  vector[I] var1; // Number of arms (integer)
  vector[I] var2; // Sample size (scaled)
  vector[I] var3; // Number of regions (integer; 0 - <= 2, 1 - >2)
  vector[I] var4; // Number of facilities (integer; 1 - <=50, 2 - >50)
  vector[I] var5; // Phase 3 (integer; 1 - yes, 0 - no)
  vector[I] var6; // Phase 4 (integer; 1 - yes, 0 - no)
  // var7-12 integer; 1 - yes, 0 - no
  vector[I] var7; // Has sites in region of the americas
  vector[I] var8; // In european region
  vector[I] var9; // In western pacific region
  vector[I] var10; // In african region
  vector[I] var11; // In south-east asia region
  vector[I] var12; // In eastern mediterranean
  vector[I] var13; // Double-blinded
  vector[I] var14; // Open label
  vector[I] var15; // Quadruple-blinded
  vector[I] var16; // Triple-blinded
  vector[I] var17; // Categorised duration (integer; 0 - <=6mths, 1 - >6mths)
  vector[I] var18; // Categorised n POMS (integer; 0 - 1 POM, 1 - >1 POM)
}

parameters {
  real mu1; // intercept
  real a_raw; // Raw shape parameter
  
  // Betas
  real beta_arm;
  vector[18] beta; // Covariate betas
  real<lower=0> tau; // Shrinkage factor
}

transformed parameters {
  // Safely transform a if near 0
  real a;
  real epsilon = 1e-4;
  a = (a_raw < 0) ? a_raw - epsilon : (a_raw > 0 ? a_raw + epsilon : 0);
  
  vector[I] b; // Rate
  vector[I] p; // Probability of attrition
  
  for (i in 1:I) {
    b[i] = exp(
      mu1 + beta_arm * arm[i] + beta[1] * var1[i] + beta[2] * var2[i] +
      beta[3] * var3[i] + beta[4] * var4[i] + beta[5] * var5[i] +
      beta[6] * var6[i] + beta[7] * var7[i] + beta[8] * var8[i] +
      beta[9] * var9[i] + beta[10] * var10[i] + beta[11] * var11[i] +
      beta[12] * var12[i] + beta[13] * var13[i] + beta[14] * var14[i] +
      beta[15] * var15[i] + beta[16] * var16[i] + beta[17] * var17[i] +
      beta[18] * var18[i]
    );
    p[i] = 1 - exp(-b[i] * (expm1(a * dt[i]) / a));
  }
}

model {
  // priors
  mu1 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  beta_arm ~ normal(0, 1);
  beta ~ normal(0, tau); // Shared shrinkage scale
  tau ~ cauchy(0, 0.5); // Weakly informative shrinkage factor
  
  // likelihood
  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;
  
  for (i in 1:I) {
    p_new[i] = 1 - exp(-b[i] * (expm1(a * dt[i]) / a));
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }
  real rate = exp(mu1); // Baseline rate
}
