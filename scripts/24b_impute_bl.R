
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(multinma)

# Import bl chars
a_imp_bl <- readRDS("processed_data/bl_chars_noimp.rds")