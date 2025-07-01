
### Description of export request #################################################
This export request contains the outputs of fitting intercept-only time-to-event models for trial attrition.

We used data from 72 trials for type 2 diabetes from the following sponsors:
- Astrazeneca (1)
- Boehringer (26)
- Eli lilly (8)
- GSK (8)
- Janssen (8)
- Sanofi (11)
- Takeda (10)

### Description of models #########################################################

The outcome was attrition, defined as participant non-completion of a trial for any reason (e.g. adverse events, lost to follow-up) before the end of protocol-defined follow-up

We fitted intercept-only models using the flexsurv R package. Models were fitted using the Gompertz distribution, specified as:
flexsurvreg(Surv(years, event) ~ 1, dist = "gompertz", data = df)

We will use the outputs of the intercept-only models to estimate Gompertz shape conditioned on trial-level covariates using data from clinicaltrials.gov

### Contents of export ############################################################

1. tidy_res.csv
These are tidied model results (coefficients and standard errors for parameters). 
One row per trial (72 in total).
* nct_id - ClinicalTrials.gov identifier
* shape_est, shape_se - Coefficient and standard error for Gompertz shape
* rate_est, rate_se - Coefficient and standard error for Gompertz rate

2. tidy_cov2cor.csv
This contains a long-format dataset of the correlation between Gompertz model parameters. This was obtained by converting the variance-covariance matrix of parameters in each trial to correlation, and transformation into a long-format dataset.
Four rows per trial (corresponding to the 2x2 correlation matrix for each trial).
* nct_id
* row, col - Gompertz model parameters in each parameter correlation matrix
* value - correlation

3. 00_functions.r - Script calling functions written for the project
4. 00_config.r - Configuration arguments for the project
5. 10_flexsurv_int.r - Script to fit intercept-only models and get model outputs