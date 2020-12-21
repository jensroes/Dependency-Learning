library(tidyverse)
library(magrittr)
library(brms)
library(emmeans)

# Load model
m <- readRDS(file = "Group comparison/group_model.rda")


# Get summary
m %>% fixef() %>% as.data.frame() %>% 
  rownames_to_column("Predictor") %>%
  filter(grepl("COND", Predictor)) %>%
  dplyr::select(-`Est.Error`) %>%

model_summary %>% mutate_if(is.numeric, round, 2); model_summary

# Save model results
write_csv(model_summary, "Group comparison/group_model_summary.csv")