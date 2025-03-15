
rm(list = ls())
gc()

library(tidyverse)
library(broom)
library(rstan)
library(bayesplot)
library(tidybayes)
library(markdown)

options(scipen = 999)

Sys.setenv(STANC3 = TRUE)