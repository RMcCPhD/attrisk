
source("scripts/00_functions_packages_etc.R")

imp_data <- read_rds("data/ipd_models.rds")
glimpse(imp_data)

# Get gompertz and log-normal model parameters
model_data <- imp_data$parameters_fit %>% 
  filter(dist == "gompertz" | dist == "lnorm") %>% 
  filter(!ctgov %in% c("NCT01131676", "NCT00274573"))

gompertz_data <- model_data %>% filter(dist == "gompertz")

gompertz_data %>% 
  ggplot(aes(x = est, colour = parameter)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free")

gompertz_data %>% 
  group_by(parameter) %>% 
  reframe(
    mean = mean(est),
    sd = sd(est)^2
  )

