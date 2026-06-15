
# Check outputs of trial-level Cox models
# Checking for proportional hazard violations and non-linearity

source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(survival)
library(survminer)
library(gridExtra)

# Import model outputs
a_imp_out <- readRDS("Processed_data/06_cox_out.rds")

# Proportional hazards ---------------------------------------------------------

# Get PH outputs (to get term names)
b_ph_out <- a_imp_out %>% 
  select(sponsor, nct_id, ph_df) %>% 
  unnest(ph_df) %>% 
  ungroup() %>% 
  filter(term != "GLOBAL") %>% 
  group_by(nct_id) %>% 
  mutate(term = paste0(unique(term), collapse = ", ")) %>% 
  ungroup() %>% 
  select(nct_id, term) %>% 
  distinct() %>% 
  separate_rows(term, sep = ", ")

# Generate Schoenfeld plots
b_scho_all <- a_imp_out %>% 
  select(sponsor, nct_id, ph) %>% 
  left_join(b_ph_out) %>% 
  mutate(
    plot_id = paste0(nct_id, "_", term),
    ph_plots = pmap(
      list(ph, term, nct_id),
      function(ph, var, trial) {
        plot <- ggcoxzph(
          ph, 
          var = var, 
          point.alpha = 0, 
          ggtheme = theme_classic()
        )
        plot_clean <- lapply(plot, function(p) {
          p + labs(title = NULL) + ylab(NULL)
        })
      }
    ),
  )

# Capture plots using grid for combining together
b_scho_capture <- b_scho_all %>% 
  mutate(
    grob = map(ph_plots, ~ grid::grid.grabExpr(print(.x)))
  )

# Split by covariate and save as PDFs
b_scho_split_all <- split(b_scho_capture, b_scho_capture$term)

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
    file_name <- paste0("Outputs/Checks/", .y, ".pdf")
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

# Combine outputs into a single PDF
b_pdf_paths <- list.files(
  "Outputs/Checks/",
  ".pdf$",
  full.names = TRUE
)

qpdf::pdf_combine(
  b_pdf_paths, 
  "Outputs/ph_schoenfeld.pdf"
)


# Non-linearity ----------------------------------------------------------------

# Target covariates
z_names <- c(
  "age10", "hba1c", "durn", "bmi", "egfr", "sbp", "chol_hdl", "chol_total",
  "comorbs", "conmeds"
)

# Select covariates of interest
# Compute martingale residuals
# Join data
c_df <- a_imp_out %>% 
  mutate(
    covs = map(drp3, ~ .x %>% select(any_of(z_names))),
    resid = map(fit3, ~ residuals(.x, type = "martingale")),
    cmbn = map2(covs, resid, ~ .x %>% cbind(.y) %>% rename(resid = ".y"))
  )

# Pivot nested datasets longer
# Unnest and renest by trial and covariate
c_transf <- c_df %>% 
  select(sponsor, nct_id, cmbn) %>% 
  unnest(cmbn) %>% 
  ungroup() %>% 
  pivot_longer(-c(sponsor, nct_id, resid)) %>% 
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
          theme_classic() +
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
    file.path("Outputs/Non_linearity/", paste0(.y, ".pdf")),
    plot = .x,
    device = "pdf",
    width = 10, height = 6, units = "in"
  )
)

e_pdf_paths <- list.files(
  "Outputs/Non_linearity",
  ".pdf$",
  full.names = TRUE
)

qpdf::pdf_combine(
  e_pdf_paths, 
  "Outputs/non_linearity_loess.pdf"
)
