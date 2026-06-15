
# Cox checks - proportional hazards, non-linearity
source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(gridExtra)
library(qpdf)
library(survival)
library(survminer)
library(splines)

# Import data, models and ph checks
a_imp_mdls <- readRDS("processed_data/cox_outputs/cox_mdls.rds")
a_imp_ph <- readRDS("processed_data/cox_outputs/cox_ph.rds")

# Proportional hazards ---------------------------------------------------------

# Get trials
b_trials <- a_imp_ph %>% 
  filter(term != "GLOBAL") %>% 
  group_by(nct_id) %>% 
  mutate(term = paste0(unique(term), collapse = ", ")) %>% 
  ungroup() %>% 
  select(nct_id, term) %>% 
  distinct() %>% 
  separate_rows(term, sep = ", ") %>% 
  left_join(a_imp_mdls)

# Get ph and plots
b_scho_all <- b_trials %>% 
  mutate(
    ph = map(cox, ~ cox.zph(.x, terms = FALSE)),
    plot_id = paste0(nct_id, "_", term),
    ph_plots = pmap(
      list(ph, term, nct_id),
      function(ph, var, trial) {
        plot <- ggcoxzph(
          ph, 
          var = var, 
          point.alpha = 0, 
          ggtheme = theme_bw()
        )
        plot_clean <- lapply(plot, function(p) {
          p + labs(title = NULL) + ylab(NULL)
        })
      }
    ),
    grob = map(ph_plots, ~ grid::grid.grabExpr(print(.x)))
  )

b_scho_split_all <- split(b_scho_all, b_scho_all$term)

walk2(
  b_scho_split_all,
  names(b_scho_split_all),
  ~ {
    
    # Get grobs and dimensions
    grob_lst <- .x$grob
    n_plots <- length(grob_lst)
    n_cols <- min(10, n_plots)
    n_rows <- ceiling(n_plots / n_cols)
    
    # Save as PNG
    file_name <- paste0("Outputs/cox_trial_level/ph/schoenfeld_plots_all/", .y, ".pdf")
    pdf(file_name, width = 20, height = 10)
    
    # Grid plot
    grid.arrange(
      grobs = grob_lst,
      ncol = n_cols,
      nrow = n_rows,
      top = .y
    )
    
    dev.off()
    
  }
)

b_pdf_paths <- list.files(
  "Outputs/cox_trial_level/ph/schoenfeld_plots_all/",
  ".pdf$",
  full.names = TRUE
)

qpdf::pdf_combine(
  b_pdf_paths, 
  "Outputs/cox_trial_level/ph/schoenfeld_combined_all.pdf"
)

# Non-linearity ----------------------------------------------------------------

# Target covariates
z_names <- c("age10", "hba1c", "durn", "bmi", "egfr", "sbp", "chol_hdl", "chol_total")

# Select covariates of interest
# Compute martingale residuals
# Join data
c_df <- a_imp_mdls %>% 
  mutate(
    covs = map(data_race, ~ .x %>% select(any_of(z_names))),
    resid = map(cox, ~ residuals(.x, type = "martingale")),
    cmbn = map2(covs, resid, ~ .x %>% cbind(.y) %>% rename(resid = ".y"))
  )

# Pivot nested datasets longer
# Unnest and renest by trial and covariate
c_transf <- c_df %>% 
  select(nct_id, cmbn) %>% 
  unnest(cmbn) %>% 
  pivot_longer(-c(nct_id, resid)) %>% 
  rename(cov_name = name) %>% 
  group_by(cov_name) %>% 
  nest(.key = "cov_df")

# Get plots
c_plot <- c_transf %>% 
  mutate(
    plots = map2(
      cov_df, cov_name,
      function(df, cov) {
        df %>% 
          ggplot(aes(x = value, y = resid)) +
          geom_point(size = 0.5, alpha = 0) +
          geom_smooth(method = "loess", se = TRUE, linewidth = 0.5) +
          facet_wrap(~nct_id, scales = "free") +
          scale_x_continuous(n.breaks = 5) +
          theme_bw() +
          theme(
            strip.background = element_rect(size = 0.05),
            strip.text = element_text(size = 5, margin = margin(2,2,2,2)),
            panel.spacing = unit(0.5, "lines"),
            axis.text.x = element_text(size = 5),
            axis.text.y = element_text(size = 5)
          ) +
          labs(x = NULL, y = NULL, title = cov)
      }
    )
  )

walk2(
  .x = c_plot$plots,
  .y = c_plot$cov_name,
  ~ ggsave(
    file.path("Outputs/cox_trial_level/non_linearity/lowess/", paste0(.y, ".pdf")),
    plot = .x,
    device = "pdf",
    width = 10, height = 6, units = "in"
  )
)

e_pdf_paths <- list.files(
  "Outputs/cox_trial_level/non_linearity/lowess/",
  ".pdf$",
  full.names = TRUE
)

qpdf::pdf_combine(
  e_pdf_paths, 
  "Outputs/cox_trial_level/lowess_all.pdf"
)

# Non-linearity (covariate vs martingale residuals) ----------------------------

# Basic working example with age10 in first trial
# Get residuals, get non-parametric relationship (LOWESS) with covariate
# Plot
tst <- a_imp_mdls$cox[[1]]
tst_df <- a_imp_mdls$data_race[[1]]
tst_resid <- residuals(tst, type = "martingale")
plot(tst_df$age10, tst_resid)
lines(lowess(tst_df$age10, tst_resid), col = "red")

# Get martingale residuals for each trial
# Join with numeric covariates per trial
# Unnest and then nest by trial and covariate
c_mtngale <- a_imp_mdls %>% 
  mutate(
    mtn_resid = map(cox, ~ residuals(.x, type = "martingale")),
    num_cov = map(data_race, ~ .x %>% select(any_of(z_names))),
    nl_df = map2(
      mtn_resid, num_cov, 
      ~ .x %>% 
        cbind(.y) %>% 
        rename(resid = ".") %>% 
        as_tibble()
    )
  ) %>% 
  select(nct_id, nl_df) %>% 
  unnest(nl_df) %>% 
  pivot_longer(-c(nct_id, resid)) %>% 
  group_by(nct_id, name) %>% 
  nest(.key = "cov_df") %>% 
  ungroup()

# Get lowess per nested covariate dataset
c_lowess <- c_mtngale %>% 
  mutate(
    lowess_cov = map(cov_df, ~ lowess(.x$value, .x$resid)),
    plot_df = map2(
      cov_df, lowess_cov, 
      ~ .x %>% 
        cbind(.y$x) %>% 
        cbind(.y$y) %>% 
        rename(
          nl_x = `.y$x`,
          nl_y = `.y$y`
        )
    )
  )

c_lowess$plot_df[[1]] %>% 
  ggplot(aes(x = value, y = resid)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) + 
  theme_classic() +
  labs(x = NULL, y = NULL)
