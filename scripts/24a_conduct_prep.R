
# Predictors
# Trttype, smoke, sex, sbp
# Race (asian, black/aa, white, other)
# Number of con meds and comorbidities
# HbA1c, eGFR, duration of T2DM
# Total and HDL cholesterol
# BMI, age10

source("scripts/00_config.R")
library(tidyverse)
library(multinma)

# Import aggregated baseline characteristics
a_imp_agg <- read_csv("data/base_dsp.csv")
