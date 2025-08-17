
## 1. Description of export request ###############################################

This export request contains model outputs from an analysis evaluating the effect 
of baseline participant characteristics on attrition risk in 74 randomised 
controlled trials for type 2 diabetes mellitus.

The 74 trials are from seven sponsors:
- AstraZeneca: 1 (1.35%)
- Boehringer Ingelheim: 28 (37.8%)
- Eli lilly: 8 (10.8%)
- GlaxoSmithKline: 8 (10.8%)
- Janssen: 8 (10.8%)
- Sanofi: 11 (14.9%)
- Takeda: 10 (13.5%)

## 2. Analysis description ########################################################

The outcome was attrition, defined as participant non-completion of a trial for any
reason (e.g. adverse events, lost to follow-up)

The goal of the analysis was to identify baseline participant characteristics that
predict attrition risk among participants receiving a placebo or an active comparator

To do this, we fitted Cox proportional hazard models to each trial:
coxph(time, event ~ characteristics)

Variables we included in the model were:
* time       | Duration of randomised participation in weeks
* event      | Attrition status (0 = completer, 1 = attrition)
* trt_class  | Antidiabetic treatment class as per ATC classification (see below)
* age10      | Age, measured in decades (i.e. divided by 10)
* bmi        | Body mass index (kg/m2)
* chol_hdl   | HDL cholesterol (mmol/L)
* chol_total | Total cholesterol (mmol/L)
* comorbs    | Number of comorbidities
* durn       | Duration of type 2 diabetes mellitus in years
* egfr       | Estimated glomerular filtration rate (ml/min/1.73m2)
* hba1c      | Glycosylated haemoglobin (%)
* n_atc      | Number of concomitant medications
* race       | Race (collapsed into Asian, Black, White and Other)
* sbp        | Systolic blood pressure (mmHg)
* sex        | Participant sex (0 = female, 1 = male)
* smoke      | Smoking status (0 = Never/Former/Not reported, 1 = Current)
* trttype    | Treatment type (placebo, active comparator)

## Contents of export request #####################################################

# Main folder
1. 01_res.csv - Model coefficients and standard errors, one row per term/trial.
* nct_id - ClinicalTrials.gov identifier
* term - Covariate term (i.e. treatment class, age, sex, interactions)
* estimate - Model log-hazard estimates
* std.error - Log-standard errors
* p.value - P value

2. 02_diag.csv - Diagnostics for fitted models, one row per model/trial.
* nct_id
* statistic.log - Likelihood ratio test statistic
* p.value.log - P-value associated with likelihood ratio test
* statistic.sc - Log-rank test statistic
* p.value.sc - P-value associated with log-rank test
* statistic.wald - Wald test statistic
* p.value.wald - P-value associated with Wald test
* statistic.robust - Robust likelihood ratio test statistic
* p.value.robust - P-value associated with robust likelihood ratio test
* r.squared - R-squared describing proportion of variation explained by model
* r.squared.max - Maximum R-squared
* concordance - Harrel's concordance (C-index)
* std.error.concordance - Standard error of C-index
* logLik - Log-likelihood
* AIC - Akaike information criterion
* BIC - Bayesian information criterion
* nobs - Number of observations used in fitted model

3. 03_cor.csv
Correlation matrix for model coefficients (not the correlation between variables).
5740 rows, one row per pair of terms (i.e. each correlation) per trial.
This was obtained by taking upper triangle for the variance-covariance matrix, 
from each model, converting to correlation and reshaping into a long data.
* nct_id
* row, col - first and second coefficient terms
* value - correlation

4. 04_prop_hazards.pdf
Plots to evaluate the proportional hazards assumption of Cox models. The plots were
created by calculating Schoenfeld residuals and plotting the median log-hazard and confidence intervals over time. The plots are grouped by covariates and assembled into a grid, with one plot per trial.

5. 05_non_linearity.pdf
Plots to evaluate non-linearity among continuous covariates fitted in Cox models.
These plots were created by obtaining Martingale residuals from fitted models and
estimating a smoothed loess curve (with confidence intervals). This visualises the relationship between a covariate and residuals. The plots are assembled in the same format as the proportional hazard plots (i.e. grid, one plot per trial, grouped by covariates). Residual points are not included in the plots to prevent re-calculation of individual values.

# Scripts
00_config.R - Configuration (i.e. disable scientific notation, clear cache)
00_functions.R - Functions written for the analysis
00_lookups.R - Prepares metadata used in the analysis
01_prepare_data.R - Data preparation for analysis
02_check_t2e.R - Check time-to-event data before analysis
03_comorbs.R - Extract comorbidities for analysis
04_conmeds.R - Extact concomitant medications for analysis
05_add_comorbs_conmeds.R - Add comorbidities and concomitant medications to dataset
06_prepare_cox.R - Create model-ready data
07_impute_within_trial.R - Perform imputation to address partially missing data
08_impute_hier.R - Imputation to address completely missing variables
09a_cox.R - Fit cox models and save outputs
09b_cox_checks.R - Model checks (proportional hazards, non-linearity)
09c_cox_ma.R - Test meta-analysis to check pooled estimates
10_export.R - Prepare model outputs and checks for export from Vivli