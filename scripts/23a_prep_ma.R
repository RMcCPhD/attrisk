
source("scripts/00_config.R")
library(tidyverse)

# Import Vivli exports
a_imp_res <- read_csv("vivli/01_res.csv")
a_imp_vcov <- read_csv("vivli/03_cor.csv")

# Get unique model terms
unique_terms <- unique(a_imp_res$term)

# Complete terms set (i.e. were missing from IPD, convg issues)
b_terms_grid <- expand.grid(
  nct_id = unique(a_imp_res$nct_id),
  term = unique_terms
) %>% as_tibble()

# Join terms grid to add missing terms
# Add zero estimates and large standard errors where missing
b_compl_terms <- a_imp_res %>%
  full_join(b_terms_grid) %>% 
  mutate(
    estimate = if_else(is.na(estimate), 0, estimate),
    std.error = if_else(is.na(std.error), 99999, std.error)
  ) %>% 
  arrange(nct_id, term)

# Get VCOV matrices ------------------------------------------------------------

# Create terms grid for correlation matrices
c_terms_vcov <- unique_terms %>% sort()

# Function to complete correlation upper triangle per trial
compl_cor_tri <- function(df, term_set) {
  
  # Canonicalise current pairs
  canon <- df %>% 
    transmute(nct_id, i = pmin(row, col), j = pmax(row, col), value) %>% 
    distinct()
  
  # Get current terms
  obs_vars <- sort(unique(c(df$row, df$col)))
  
  # Include missing terms
  vars <- sort(union(obs_vars, term_set))
  
  # Build grid
  if (length(vars) < 2) {
    return(
      tibble(
        nct_id = unique(df$nct_id),
        row = character(),
        col = character(),
        value = numeric()
      )
    )
  }
  grid <- as_tibble(t(combn(vars, 2)))
  names(grid) <- c("i", "j")
  
  # Left-join observed terms, rename as row/col
  grid %>% 
    mutate(nct_id = unique(df$nct_id)) %>% 
    left_join(canon, by = c("nct_id", "i", "j")) %>% 
    mutate(value = if_else(is.na(value), 0, value)) %>% 
    transmute(nct_id, row = i, col = j, value)
  
}

# Apply function to each trial
c_compl_tri <- function(a_imp_vcov, c_terms_vcov) {
  a_imp_vcov %>% 
    group_by(nct_id) %>% 
    group_split() %>% 
    map_dfr(compl_cor_tri, c_terms_vcov)
}

c_compl_cor <- c_compl_tri(a_imp_vcov, c_terms_vcov)

# List unique covariate names per trial
c_cor_vars <- c_compl_cor %>% 
  pivot_longer(row:col, values_to = "variable") %>% 
  distinct(nct_id, variable) %>% 
  group_by(nct_id) %>% 
  reframe(vars = list(unique(variable)))

# Bind correlation and rename components for upper and lower triangles
# Create diagonals (i.e. cor = 1)
c_cor_parts <- c_compl_cor %>% 
  bind_rows(c_compl_cor %>% rename(col = row, row = col)) %>% 
  group_by(nct_id) %>% 
  group_split() %>% 
  map_dfr(function(data_part) {
    key <- data_part %>% select(nct_id) %>% slice(1)
    vars <- unique(c(data_part$row, data_part$col))
    
    # Get diagonals
    diags <- expand.grid(row = vars, col = vars) %>% 
      filter(row == col) %>% 
      mutate(value = 1) %>% 
      bind_cols(key[rep(1, nrow(.)), ])
    
    # Bind diagonals and triangles
    bind_rows(data_part, diags)
  })

# Construct correlation matrices
c_cor_construct <- c_cor_parts %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  mutate(
    cor_matrix = map(data, ~ {
      vars <- union(.x$row, .x$col)
      cor_tbl <- xtabs(value ~ row + col, data = .x)[vars, vars]
      matrix(
        as.numeric(cor_tbl),
        nrow = length(vars),
        ncol = length(vars),
        dimnames = list(vars, vars)
      )
    })
  ) %>% 
  select(-data)

# Get VCOV matrices
c_vcov <- c_cor_construct %>% 
  mutate(
    cov_matrix = pmap(list(nct_id, cor_matrix), function(id, cor_mat) {
      se_df <- b_compl_terms %>% filter(nct_id == id)
      se_vec <- setNames(se_df$std.error, se_df$term)
      se_order <- se_vec[rownames(cor_mat)]
      cov_mat <- outer(se_order, se_order) * cor_mat
      cov_mat
    })
  )

# Quick check that dimensions are correct for each (16x16)
for (i in seq_along(c_vcov$cov_matrix)) {print(dim(c_vcov$cov_matrix[[i]]))}

# Save outputs
saveRDS(b_compl_terms, "processed_data/conduct_cox_res.rds")
saveRDS(c_vcov, "processed_data/conduct_cox_vcov.rds")