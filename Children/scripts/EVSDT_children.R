library(tidyverse)
#library(rstanarm) # stan_lmer
#library(polspline) # For BF
require(purrr)   # for data manipulations
library(magrittr)
library(scales)  # for percentage scales
library(brms)
library(rethinking)
source("functions/functions.R")

# Load data
d <- read_csv("data/expl_learning.csv") %>%
  mutate(dependency = factor(dependency))

d <- d %>% 
  mutate(sayDep = given == "k",
         isDep = expected == "k") %>% 
  dplyr::select(-expected, -given)

d %>% group_by(dependency) %>%
  summarise(mean = mean(acc),
            sd = sd(acc))
d %>% 
  group_by(isDep, dependency) %>%
  summarise(N = n(),
            M = mean(sayDep),
            SE = sqrt((mean(sayDep)*(1 - mean(sayDep)))/length(sayDep))) %>%
  ungroup() %>%
  mutate(dependency = ifelse(dependency == "adjacent", "Adjacency", "Non-adjacency"),
         isDep = ifelse(isDep == TRUE, "Hit", "False alarm")) %>%
  ggplot(aes(y = M, x = isDep, color = isDep)) +
  geom_line() +
  geom_point(shape = 1, size = 3.5) +
  labs(y = bquote("Percentage of"~italic("yes")~"responses with 95% CIs")) +
  facet_grid(~dependency) +
  scale_y_continuous(labels=scales::percent,breaks = c(.2,.3,.4,.5,.6,.7), limits = c(.4,.7)) +
  geom_hline(yintercept = .5, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin=M-1.96*SE, ymax=M+1.96*SE), width=.1, size = .75) +
  theme_minimal() +
  scale_color_manual("", values = c("forestgreen", "darkred")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16))

ggsave("plots/expljudge2.png",  width = 7, height =6)

# Fit model
nchains = ncores = 3
iter = 8000

# We can estimate the EVSDT model’s parameters for every subject and the population average 
# in one step using a Generalized Linear Mixed Model (GLMM).  Rouder and Lu (2005) and Rouder et al. (2007) 
# discuss hierarchical modeling in the context of signal detection theory.
# Dep describes the difference in c between dependencies, and the interaction term isdep:dep
# would describe the difference in d’ between groups. 
#https://vuorre.netlify.com/post/2017/10/12/bayesian-estimation-of-signal-detection-theory-models-part-2/

m <- bf(sayDep ~ Phi(dprime*isDep - c), 
         dprime ~ dependency + (dependency |s| subj) + (dependency |i| item), 
         c ~ dependency + (dependency |s| subj) + (dependency |i| item), 
         nl = TRUE)

priors <- c(prior(normal(0, 3), nlpar = "dprime"), # lb = 0
            prior(normal(0, 3), nlpar = "c"),
            prior(student_t(10, 0, 1), class = "sd", nlpar = "dprime"),
            prior(student_t(10, 0, 1), class = "sd", nlpar = "c"),
            prior(lkj(2), class = "cor"))

fit <- brm(m, 
            family = bernoulli(link="identity"), 
            data = d,
            prior = priors,
            control = list(adapt_delta = .99),
            cores = ncores,
            chains = nchains,
            iter = iter)

summary(fit)
pp_check(fit)
plot(fit, "^b_")
marginal_effects(fit)

# Save model
save(fit,
     file="stanout/EVSDT_children.rda",
     compress="xz")

#load(file="stanout/EVSDT_adults.rda")

# Model checks
samps <- brms::posterior_samples(fit, "b_") # It saves all the samples from the model.

write_csv(samps, "stanout/EVSDT_children_posterior.csv")

samps %>% 
  gather(key = param, value = value) %>%
  group_by(param) %>%
  summarise(Med = median(value),
            Lo = PI(value, prob = .95)[1],
            Up = PI(value, prob = .95)[2]) %>%
  mutate(param = c("-c (adjacent - non-adjacent)", "-c", "d' (adjacent - non-adjacent)", "d'")) %>%
  mutate(id = nchar(param)) %>%
  arrange(id, param) %>%
  dplyr::select(-id) %>%
  round_df(2) -> M;M


rownames(M) <- c("-\\textit{c}", "\\textit{d}'", "\\textit{c} (adjacent-non adjacent)","\\textit{d}' (adjacent-non adjacent)" )
colnames(M) <- c("Predictor", "\\hat{\\mu}", "PI [2.5\\%", "97.5\\%]")

library(xtable)
print(xtable(M, caption = "", label = "tab:expljudg"),
      table.placement = "!h",
      caption.placement = "top",  include.rownames = TRUE,  
      sanitize.rownames.function = identity, sanitize.colnames.function = identity, sanitize.text.function =identity)
             
