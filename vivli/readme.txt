
This export request contains trial-level model outputs and summaries from an analysis evaluating the effect of baseline characteristics on attrition risk in 76 randomised controlled trials for type 2 diabetes mellitus.

This export is an updated version of a previous export requested 2025-08-14, with the following additions:
* Two new predictors added to the modelling: number of comorbidities, and number of concomitant medications

The trials involved in this analysis were obtained from seven sponsors:
* AstraZeneca: 2 
* Boehringer Ingelheim: 27 
* Eli-Lilly: 8 
* GlaxoSmithKline: 8 
* Janssen: 9
* Sanofi: 12
* Takeda: 10

File structure
Export/
| - Models/
| - Plots/
| - Scripts/
| - Summaries/
| - readme.txt

---------------------------------------------------

# Description of analysis

* Outcome: Attrition (participant non-completion of a trial for any reason, e.g. adverse events)
* Aim: Determine the effect of baseline participant characteristics on attrition risk

This analysis fitted Cox proportional hazard models to the control arms of each of the 76 trials
* Model: coxph(time, event ~ characteristics)
* Covariates are centred on means and scaled by standard deviations
* e.g. (age10 - age10mean) / age10sd

The following variables were specified for each Cox model:
* time       | Duration of follow-up in weeks
* event      | Attrition status (0 = completer, 1 = attrition)
* trt_class  | Non-experimental comparator (0 = placebo, 1 = active comparator)
* age10      | Age, measured in decades (i.e. divided by 10)
* bmi        | Body mass index (kg/m2)
* chol_hdl   | HDL cholesterol (mmol/L)
* chol_total | Total cholesterol (mmol/L)
* comorbs    | Number of comorbidities
* durn       | Duration of type 2 diabetes mellitus in years
* egfr       | Estimated glomerular filtration rate (ml/min/1.73m2)
* hba1c      | Glycosylated haemoglobin (%)
* conmeds    | Number of concomitant medications
* race       | Race (collapsed into Asian, Black, White and Other)
* sbp        | Systolic blood pressure (mmHg)
* sex        | Participant sex (0 = female, 1 = male)
* smoke      | Smoking status (0 = Never/Former, 1 = Current)
* trttype    | Treatment type (placebo, active comparator)

---------------------------------------------------

# Export/Models/ - Outputs for the trial-level Cox models

Models/
| - export_covcor.csv
| - export_diag.csv
| - export_res.csv

1. Models/export_covcor.csv

Description:
* Upper triangle of the coefficient correlation matrices from fitted Cox models
* These are not correlations between covariates
* Produced by converting variance-covariance matrices to correlation matrices
* Only the upper triangle of correlation matrices retained to reduce file size
* Reshaped into long format for export

Variables:
* nct_id | ClinicalTrials.gov identifier
* row    | First coefficient term
* col    | Second coefficient term
* value  | Correlation between coefficient terms, rounded to 3 decimal places

===========================

2. Models/export_diag.csv

Description:
* Diagnostics for fitted Cox models
* One row per trial

ROUNDING

Variables:
* nct_id                | ClinicalTrials.gov identifier
* statistic.log         | Likelihood ratio test statistic
* p.value.log           | P-value associated with likelihood ratio test
* statistic.sc          | Log-rank test statistic
* p.value.sc            | P-value associated with log-rank test
* statistic.wald        | Wald test statistic
* p.value.wald          | P-value associated with Wald test
* statistic.robust      | Robust likelihood ratio test statistic
* p.value.robust        | P-value associated with robust likelihood ratio test
* r.squared             | R-squared describing proportion of variation explained by model
* r.squared.max         | Maximum R-squared
* concordance           | Harrel's concordance (C-index)
* std.error.concordance | Standard error of C-index
* logLik                | Log-likelihood
* AIC                   | Akaike information criterion
* BIC                   | Bayesian information criterion
* nobs                  | Number of observations used in fitted model

===========================

3. Models/export_res.csv

Description:
* Model coefficients, standard errors, and significance tests
* One row per term and trial

Variables:
* nct_id    | ClinicalTrials.gov identifier
* term      | Coefficient term
* estimate  | Log-hazard estimate, rounded to 3 decimal places
* std.error | Log-standard error, rounded to 3 decimal places
* p.value   | P value, rounded to 3 decimal places

---------------------------------------------------

# Export/Plots/ - Plots supporting methodological approaches

Plots/
| - cox_non_linearity_loess.pdf
| - cox_ph_schoenfeld.pdf

# 1. Plots/cox_non_linearity_loess.pdf

Description:
* Evaluates non-linearity among continuous covariates in fitted Cox models
* Created by calculating Martingale residuals from fitted models
* Residuals used to estimate a smoothed loess curve with confidence intervals
* Plots grouped by covariates and assembled into a grid (one cell per trial)
* Residuals were removed to preserve privacy

===========================

# 2. Plots/cox_ph_schoenfeld.pdf

Description:
* Evaluates the proportional hazard assumption of fitted Cox models
* Created by calculating Schoenfeld residuals
* Plots Schoenfeld residuals against median log-hazard with confidence intervals over time
* Plots grouped by covariates and assembled into a grid (one cell per trial)
* Residuals were removed to preserve privacy

---------------------------------------------------

# Export/Scripts/ - R scripts

Description
* Scripts written to prepare data and perform the analysis

Scripts/
| - 00_config.R
| - 00_functions.R
| - 00_lookups.R
| - 01_prepare_data.R
| - 02_check_t2e.R
| - 03_prepare_cox.R
| - 04_impute_within_trial.R
| - 05_impute_hier.R
| - 06_cox.R
| - 07_pool.R
| - 08_checks.R
| - 09_export.R

Overview:
* 00_config.R              | Setup/Configuration for RStudio environment
* 00_functions.R           | Functions used for data preparation and analysis
* 00_lookups.R             | Generates lookups (outputs not exported) used in workflow
* 01_prepare_data.R        | Prepares analysis data
* 02_check_t2e.R           | Checks time-to-event data
* 03_prepare_cox.R         | Finalises data preparation
* 04_impute_within_trial.R | Within-trial missing data imputation
* 05_impute_hier.R         | Across-trial missing data imputation
* 06_cox.R                 | Fit Cox models
* 07_pool.R                | Check pooled estimates (outputs not exported)
* 08_checks.R              | Check proportional hazard assumption and non-linearity
* 09_export.R              | Prepare export request

---------------------------------------------------

# Export/Summaries/ - Summary statistics

Summaries/
| - missingness_imputation_stages.csv
| - overall_summary.csv
| - trial_level_summary.csv
| - covariate_mu_d.csv

1. Summaries/missingness_imputation_stages.csv

Description:
* Overall summary of missingness addressed for each covariate across imputation stages
* One row per covariate, percentages rounded to 2 decimal places

Variables:
* term      | Coefficient term
* na_start  | Missingness (n [%]) before any imputation
* na_stg1   | Remaining missingness (n [%]) after within-trial multiple imputation
* na_stg2   | Confirmation of no missingness (n [%]) remaining after across-trial imputation

===========================

2. Summaries/overall_summary.csv

Description:
* Overall summary of study characteristics
* Overall summary of baseline characteristics not previously analysed
* Baseline characteristics are number of comorbidities and concomitant medications

Variables:
* name      | Characteristic summarised
* value     | Summary statistic

Summaries:
* n_trials    | Number of trials
* n_subj      | Number of randomised participants
* attr        | Total number of attrition events
* risk        | Overall attrition risk as median [min, max], rounded to 2 decimal places
* pt_yr_total | Total person-time (years), rounded to 2 decimal places
* followup    | Follow-up in weeks as median [min, max]
* comorbs     | Number of comorbidities as median [min, max]
* conmeds     | Number of concomitant medications as median [min, max]

===========================

3. Summaries/trial_level_summary.csv

Description:
* Trial-level summary of attrition
* One row per arm in trial

Variables:
* nct_id      | ClinicalTrials.gov identifier
* trtcls5     | Arm treatment class
* trttype     | Arm treatment type (either placebo or active comparator)
* n           | Number of randomised participants
* attr        | Number of non-completers
* attr_p      | Attrition (%) rounded to 2 decimal places
* pt_yr_total | Arm-level total person-time (years), rounded to 2 decimal places

===========================

4. Summaries/covariate_mu_d.csv

Description:
* Mean and standard deviation across all trials for covariates
* One row per covariate, values rounded to 2 decimal places

Variables:
* term   | Covariate name
* mu     | Overall covariate mean (across all trials)
* d      | Overall standard deviation (across all trials)   


