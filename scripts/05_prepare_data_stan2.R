
# Preparing STAN data for trial sub-model
# Addition of participant characteristics also in IPD
# Looking for:
#' Age, sex, race
#' Smoking status
#' Duration of diabetes
#' BMI
#' Blood pressure
#' eGFR
#' Cholesterol
#' Comorbidities
#' Concomitant medications

source("scripts/00_functions_packages_etc.R")

a_imp_df <- readRDS("processed_data/dat_covars.rds")