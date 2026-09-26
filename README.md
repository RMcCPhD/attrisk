# Predicting aggregate-level and participant-level attrition risk in newer glucose-lowering drug trials: Model development and internal validation

This repository contains the analysis code, non-individual-level data, metadata, and selected outputs supporting the manuscript *Predicting aggregate-level and participant-level attrition risk in newer glucose-lowering drug trials: Model development and internal validation.*

The study develops a framework for predicting control arm attrition risk in phase III/IV randomised controlled trials of newer glucose-lowering drugs for type 2 diabetes mellitus. The framework comprises:

- a Bayesian design-stage model linking aggregate control arm attrition counts to trial design characteristics and follow-up duration through a Gompertz cumulative distribution function;
- internal validation of the design-stage model using ten-fold cross-validation;
- a conduct-stage model incorporating participant-level associations obtained from trial-specific Cox models and multivariate random-effects meta-analysis;
- numerical integration over generated participant profiles and a log-rate intercept adjustment that preserves design-stage predicted risk at each control arm's follow-up; and
- sensitivity analyses of the integration grid size and illustrative participant-level risk trajectories within a fixed trial design.

The analysis includes 342 trials with 72,692 control arm participants, of which 76 trials contributed individual participant data (IPD) for 19,023 participants. The conduct-stage alignment checks assess agreement with design-stage predictions; the conduct-stage model was not internally validated, and the complete framework has not been externally validated.

The associated AttRisk Shiny interface, described in Chapter 7, is maintained separately in [`attrisk_shiny`](https://github.com/RMcCPhD/attrisk_shiny).

## Reproducibility and data access

The public materials include aggregate data and metadata in `data/`, `created_metadata/`, and `scratch_data/`, together with the disclosure-reviewed exports in `vivli/`. Selected tabular summaries and figures are retained in `outputs/` so that results can be inspected without rerunning the models.

The raw IPD are not included and cannot be shared publicly. They were accessed through the Vivli secure research environment. The scripts in `vivli/Scripts/` document data preparation, imputation, and trial-specific Cox modelling in that environment. Rerunning those scripts requires authorised access to the source IPD. Approved non-individual-level exports, including model coefficients, coefficient correlations, diagnostics, and baseline summaries, are provided for downstream analyses. Generated participant profiles in the public workflow are numerical integration points, not records of actual participants.

The design-stage model and ten-fold cross-validation were run on a high-performance computing (HPC) system. Both sets of materials are supplied in `hpc/`, including the R scripts, Stan model definitions, staged aggregate-level `stan_data.rds` inputs, and SLURM submission scripts. These allow the fitting stages to be rerun without first reconstructing the upstream data-preparation workflow, subject to the setup adjustments described below. Fitting can also be performed on a suitably configured local machine, but is computationally intensive.

Large fitted objects and draw-level intermediates are not given in this checkout. Downstream scripts require the following objects to be regenerated:

- the full design-stage fit, stored locally as `outputs/design_stage/model/design_model.rds`;
- the ten cross-validation fits and associated fold metadata under `hpc/design_cross_val/`; and
- conduct-stage fitted objects, imputed baseline summaries, generated profiles, and intercept-adjustment draws.

The retained CSV files are summaries rather than substitutes for posterior draws. Full reproduction from the preparation scripts additionally requires the external inputs listed below. Several scripts retain paths from the original working environment; this repository is an analysis record with rerun materials, rather than an automated end-to-end pipeline.

## Relationship to the preceding attrition and network meta-analyses

The initial trial set comes from [`attrition_differential_t2d`](https://github.com/RMcCPhD/attrition_differential_t2d); `scripts/01a_prepare_agg.R` generates `processed_data/tidied_agg.rds`. Copy that generated file to `data/tidied_agg.rds` in this repository before running `scripts/01_prepare_agd.R` as this is not supplied here.

The baseline summary preparation and profile generation also use materials from the preceding network meta-analysis and related data-preparation projects. The scripts assume the following sibling-folder arrangement:

``` text
parent_directory/
|-- attrisk/
|-- attrition_differential_t2d/
|-- nma_agesex_public/
|-- nma_more_covariates/
|-- additional_covariates/
`-- project_data/
```

The external files referenced by the scripts are:

| Source folder | Required files | Used by |
|------------------------|------------------------|------------------------|
| `nma_agesex_public/` | `Scratch_data/agg_ipd_hba1c.Rds` | `08a_sobol_prepare_ipd.R` |
| `nma_more_covariates/` | Covariate-summary files matching `Data/*covariates*`, plus `Data/ipd_age_sex.csv`, `Data/categorical_variables_percentage.csv`, and `Data/hba1c_base_summary_ipd.csv` | `08a_sobol_prepare_ipd.R` |
| `nma_more_covariates/` | `Scratch_data/correlation_summary.rds` | `08c_sobol_generate.R` |
| `additional_covariates/` | `Outputs/baseline_characteristics_trial_wide.csv` | `08b_sobol_prepare_agd.R` |
| `project_data/` | `data/aact/detailed_descriptions.csv` | `02_explore_gompertz.R` |

The preceding public network meta-analysis is available at [`nma_agesex_public`](https://github.com/Type2DiabetesSystematicReview/nma_agesex_public). Some dependencies are generated intermediate files, and the other source folders are not bundled with this repository. Confirm that the required files are available, or adapt the input paths to equivalent authorised exports, before running the dependent stages. A fresh clone of this repository alone does not provide all inputs for the full preparation and conduct-stage workflow.

## Repository structure

### `created_metadata/`

An antidiabetic treatment-class lookup (`atc.csv`) used to harmonise control arm treatments during data preparation.

### `data/`

Extracted trial eligibility criteria and trial-level Gompertz parameter summaries used in the exploratory analyses. The initial aggregate attrition dataset, `tidied_agg.rds`, must be obtained from the preceding analysis as described above.

### `hpc/`

Materials for the computationally intensive design-stage model and internal validation. `hpc/design_model/` contains the full-model fitting bundle, and `hpc/design_cross_val/` contains the ten-fold cross-validation bundle. Each includes a staged Stan dataset, model definition, R script, and SLURM submission script. Cluster-specific settings and the preparation steps below must be reviewed before submission.

### `outputs/`

Selected machine-readable results and figures. `design_stage/` contains exploratory Gompertz results, posterior summaries, development performance, and cross-validation summaries. `conduct_stage/` contains baseline summaries, pooled participant-level associations, intercept-adjustment diagnostics, integration-grid sensitivity results, and illustrative profile predictions. `pmsampsize/` records the exploratory sample-size assessment, and `trial_summaries/` contains descriptive summaries. Additional fitted objects, intermediate files, tables, and figures are generated here when scripts are run.

### `processed_data/`

The destination for analysis-ready objects generated by the preparation scripts, including the prepared aggregate dataset, design matrix, exploratory shape-rate relationship, and Stan input list. These intermediates are not supplied in the checkout; create this directory before running the preparation workflow.

### `scratch_data/`

Checked trial characteristics, manually completed missing information, covariate labels, and distribution choices used during preparation and reporting. These are analysis inputs and records of processing decisions, rather than disposable temporary files.

### `scripts/`

Data preparation, model specification, fitting, validation, conduct-stage construction, and reporting code. The intended sequence and dependencies are summarised below.

### `vivli/`

Disclosure-reviewed exports from the secure research environment. `Models/` contains Cox model coefficients, diagnostics, and the upper triangles of coefficient-correlation matrices used to reconstruct within-trial variance-covariance matrices. `Summaries/` contains trial summaries, covariate scaling information, and missingness summaries. `Cox_ph_nonlinearity/` contains exported diagnostic plots. `Scripts/` records the restricted IPD workflow, and `readme.txt` describes the exports.

## Analysis workflow

Run preparation and post-processing scripts from the repository root so that relative paths resolve correctly. The standalone fitting scripts instead expect a staged working directory, as described in the HPC instructions. Script filenames in the table are relative to `scripts/`.

| Stage | Script | Purpose |
|------------------------|------------------------|------------------------|
| Setup | `00_config.R` | Clears the R workspace and sets global options. |
| Setup | `00_functions_packages_etc.R` | Loads packages and defines shared preparation, prediction, and validation functions. |
| Sample-size assessment | `00_protocol_ss.R` | Explores the number of candidate predictors supported under alternative sample-size assumptions. |
| Data preparation | `01_prepare_agd.R` | Selects and combines control arms, prepares trial design characteristics and eligibility criteria, and constructs the design matrix. |
| Exploratory analysis | `02_explore_gompertz.R` | Examines trial-level Gompertz parameters and estimates the shape-log-rate relationship used by the design-stage model. |
| Stan preparation | `03_prepare_stan.R` | Combines the prepared data, design matrix, and exploratory estimates into `processed_data/stan_data.rds`. |
| Model specification | `04a_design_model.stan` | Defines the Bayesian binomial-Gompertz model, arm-specific intercepts, and coefficient shrinkage. |
| Design-stage fitting | `04b_design_fit.R` | Fits the design-stage model; the staged HPC equivalent is `hpc/design_model/design_model_full.R`. |
| Design-stage checks | `04c_design_checks.R` | Produces sampling diagnostics, posterior predictive checks, and development performance summaries. |
| Design-stage reporting | `04d_design_sum.R` | Summarises posterior parameters and trial design associations. |
| Internal validation | `05a_design_validate.R` | Refits the model across ten held-out folds; the staged HPC equivalent is `hpc/design_cross_val/design_validate.R`. |
| Validation reporting | `05b_design_validate_output.R` | Reconstructs held-out predictions and summarises error, calibration, and concordance on count and risk scales. |
| Conduct-stage associations | `07_conduct_ma.R` | Reconstructs within-trial covariance matrices and pools Cox coefficients using multivariate random-effects meta-analysis. |
| IPD baseline preparation | `08a_sobol_prepare_ipd.R` | Prepares exported IPD baseline summaries for generating participant profiles. Requires the prepared aggregate dataset and external summaries. |
| Aggregate baseline preparation | `08b_sobol_prepare_agd.R` | Combines baseline summaries, imputes missing aggregate characteristics across 20 datasets, and creates reference distributions. Run after `08a_sobol_prepare_ipd.R`. |
| Profile generation | `08c_sobol_generate.R` | Generates 64 Sobol integration profiles per control arm per imputation and computes participant-level predictors. Requires pooled coefficients from `07_conduct_ma.R` and completed summaries from `08b_sobol_prepare_agd.R`. |
| Conduct-stage adjustment | `09a_citl_all.R` | Calculates draw-specific log-rate intercept adjustments to align participant-averaged risk with design-stage predictions. Requires the design-stage fit and generated profiles. |
| Sensitivity analysis | `09b_sobol_grid_stability.R` | Compares integration grids of 64, 128, 256, and 512 profiles. Run after the primary intercept adjustment. |
| Illustrative predictions | `10_illustrative_profiles.R` | Generates risk trajectories for hypothetical participant profiles within a fixed 52-week trial design. Requires the design-stage fit, profiles, and intercept-adjustment draws. |
| Reporting | `11_sum_trials.R` | Produces trial, attrition, follow-up, baseline, and Cox model descriptive summaries. Run after aggregate and IPD baseline preparation. |

### Preparing a fresh checkout

Create the `processed_data/` directory and any missing output sub directories before running scripts that write to them. In addition to obtaining the external inputs, reconcile these retained working-environment references:

- In `01_prepare_agd.R`, change `vivli/v2/Models/` and `vivli/v2/Summaries/` to the supplied `vivli/Models/` and `vivli/Summaries/` directories.
- In `04d_design_sum.R`, point the fitted-model import at `outputs/design_stage/model/design_model.rds` instead of `outputs/design/design_model.rds`.
- In `07_conduct_ma.R`, the model is fitted and then replaced by an import of `outputs/conduct_stage/ipd_meta_analysis/fit.rds`, which is not supplied. For a fresh fit, save the newly fitted `ma_fit` there and omit the subsequent cached-fit import.
- In `08c_sobol_generate.R`, replace the legacy `outputs/conduct/ma/out/ma_tidy.rds` import with the result produced by `07_conduct_ma.R`, which currently writes `outputs/conduct_stage/ipd_meta_analysis/ma_tidy.csv`. Use `read_csv()` for that CSV, or explicitly save and read an equivalent RDS object.

After these adjustments and once the required inputs are available, the preparation sequence is:

``` sh
Rscript scripts/01_prepare_agd.R
Rscript scripts/02_explore_gompertz.R
Rscript scripts/03_prepare_stan.R
```

The staged HPC inputs provide an alternative starting point for design-stage fitting and cross-validation. They do not replace the labelled aggregate data and baseline summaries needed by all downstream reporting and conduct-stage scripts.

### Running the design-stage model and cross-validation on a cluster

The two `hpc/` directories contain identical staged `stan_data.rds` inputs. The full model and each validation fold use four chains, 1,000 warmup iterations and 3,000 post-warmup iterations per chain, with `adapt_delta = 0.99` and `max_treedepth = 15`. Cross-validation assigns the collapsed control-arm rows to ten folds and refits the model for each held-out fold.

Before submitting on another system:

1.  Review the inputs and R scripts in each bundle. Run each job from its own bundle directory. If rebuilding the inputs, copy the same regenerated `processed_data/stan_data.rds` into both bundles.

2.  In `hpc/design_model/design_model_full.R`, change `stan_model("04a_design_model.stan")` to `stan_model("design_model.stan")` to match the supplied filename.

3.  In `hpc/design_cross_val/design_validate.R`, add the following immediately after `out_dir <- "folds"`. The current script uses `folds_models` without defining it and assumes both output directories exist:

    ``` r
    folds_models <- "folds_models"
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(folds_models, recursive = TRUE, showWarnings = FALSE)
    ```

4.  Adapt the account, partition, R module, library configuration, time limit, memory, and CPU request in each `slurm_2026.sh`. Install the required R packages and Stan toolchain on the target system, and ensure R's parallel-core setting matches the allocated CPUs.

5.  Submit `sbatch slurm_2026.sh` from each respective bundle directory. Alternatively, run `Rscript design_model_full.R` and `Rscript design_validate.R` from their respective directories on a suitably configured local machine after making the same adjustments.

6.  Copy the full-model output `design_model.rds` to `outputs/design_stage/model/design_model.rds`. Retain the cross-validation `folds/` and `folds_models/` directories under `hpc/design_cross_val/`, matching the paths used by `05b_design_validate_output.R` and its helper functions.

7.  Place the matching Stan input at `processed_data/stan_data.rds` for local post-processing. Run `04c_design_checks.R`, `04d_design_sum.R`, and `05b_design_validate_output.R` once their other inputs and output directories are ready.

The supplied SLURM scripts record the original submission setup and are not expected to run unchanged on another cluster. The standalone fitting scripts in `scripts/` likewise require input staging and directory setup; they are not root-level pipeline entry points.

## Software requirements

Chapter 6 reports R 4.6.0 and Stan for the analysis. The supplied HPC scripts load the cluster's `apps/R` module without pinning a version. Exact package versions are not locked in this repository.

The public preparation, modelling, and reporting scripts use:

- `tidyverse`
- `broom`
- `rstan`
- `bayesplot`
- `tidybayes`
- `posterior`
- `markdown`
- `flexsurv`
- `metafor`
- `multinma`
- `mice`
- `pmsampsize`
- `cowplot`
- `Matrix`, `glue`, and `scales`

Stan fitting requires a working C++ toolchain appropriate to the operating system. Profile generation uses `multinma::add_integration()`. The secure-environment scripts have additional IPD preparation, imputation, and survival-analysis dependencies recorded in `vivli/Scripts/`.

## Statistical analysis plan

The pre-specified study protocol is maintained separately in [`attrisk_protocol`](https://github.com/RMcCPhD/attrisk_protocol). Chapter 6 of the thesis describes the implemented framework, and Appendix 5 provides the detailed missing-data methods, exploratory Gompertz analyses, design-stage specification, and conduct-stage construction.

## Citation

TBC
