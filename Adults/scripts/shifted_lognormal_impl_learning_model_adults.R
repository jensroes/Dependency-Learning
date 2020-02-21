library(MASS)
library(tidyverse)
library(rethinking)
library(brms)
source("functions/functions.R")

# Load data
d <- read_csv("data/seqlearn.csv") %>%
  mutate(dependency = factor(dependency))

table(d$dependency2, d$block)

#d <- d[d$trial %in% filter(d, antecedent != "none" & acc == 1)$trial,]
d <- d %>% filter(acc == 1
                  , rt > 100
                  , rt < 10000
                  , !is.na(dependency2)
                  , dependency != "both"
)
table(d$dependency, d$block)

d$dependency <- factor(d$dependency)
d$block <- as.factor(d$block)
with(d, table(dependency, block))

# Determine contrasts
d$COND <- paste(d$dependency, d$block, sep = "_")
d$COND <- as.factor(d$COND)
summary(d$COND)

cmat <- fractions(matrix(c(
  # Main effects
  # Sum coding on adjacency
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2,-2,-2,-2,  # Dependency
  -1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,  # Adjacency
  # Forward difference coding for block
  1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,  # Block 1-2
  0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,  # Block 2-3
  0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0,  # Block 3-4
  0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0,  # Block 4-5
  0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0,  # Block 5-6
  0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1,  # Block 6-7
  # Interactions
  1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,  # Dep * Block 1-2
  -1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # Adj * Block 1-2
  0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-2, 2, 0, 0, 0, 0, # Dep * Block 2-3
  0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # Adj * Block 2-3
  0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-2, 2, 0, 0, 0, # Dep * Block 3-4
  0, 0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # Adj * Block 3-4
  0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-2, 2, 0, 0,  # Dep * Block 4-5
  0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, # Adj * Block 4-5
  0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-2, 2, 0, # Dep * Block 5-6
  0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, # Adj * Block 5-6
  0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0,-2, 2, # Dep * Block 6-7
  0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0 # Adj * Block 6-7
), nrow=21, byrow=F))

rown <- levels(d$COND)
rownames(cmat) <- rown

dep <- c("Dependency", "Adjacency") 
block <- c("Block 1-2", "Block 2-3", "Block 3-4", "Block 4-5", "Block 5-6", "Block 6-7")
coln <- c(dep, block)
for(bl in block){
  coln <- c(coln, paste(dep, bl, sep = " * "))
}
colnames(cmat)<- coln
cmat; colSums(cmat)

# transpose and inverse matrix to get contrast between the expected levels
inv.cmat<-fractions(t(ginv(cmat)))
rownames(inv.cmat) <- rown
colnames(inv.cmat) <- coln
inv.cmat

# Assign contrasts
contrasts(d$COND) <- inv.cmat

# Fit model
nchains = ncores = 3
iter = 8000

m <- bf(formula = rt ~ COND + scale(trial_in_block) + 
          (COND + scale(trial_in_block) | subj) + 
          (1 | stim), 
        ndt ~ 1 )

get_prior(m, data = d, family = shifted_lognormal())

priors <- c(prior(normal(6, 2), class = "Intercept"),
            prior(normal(0, 3), class = "b"),
            prior(normal(0, 5), class = "sd"),
            prior(lkj(3), class = "cor"),
            prior(cauchy(0,2.5), class = "Intercept", dpar = "ndt")
)

fit <- brm(m, 
           data = d,
           prior = priors,
           family = shifted_lognormal(),
           control = list(adapt_delta = .99),
           cores = ncores,
           chains = nchains,
           iter = iter)


# Save model
save(fit,
     file="stanout/shifted_lognormal_impl_learning_model_adults.rda",
     compress="xz")


summary(fit)
plot(fit, pars = "^b_") 
pp_check(fit)

marginal_effects(fit)

# Model checks
samps <- brms::posterior_samples(fit, "b_") # It saves all the samples from the model.

write_csv(samps, "stanout/shifted_lognormal_impl_learning_adults_posterior.csv")

samps <- read_csv("stanout/shifted_lognormal_impl_learning_adults_posterior.csv")


names(samps) <- c("Int", "Int2", coln, "trial_in_block")

samps %>%
  as_tibble() %>%
  dplyr::select(-trial_in_block, -Int2) %>%
#  mutate_at(-1, list(~ (exp(Int + .) - exp(Int)))) %>%
  dplyr::select(-Int) %>%
  gather(key = param, value = value) %>%
  mutate(param = factor(param, levels = rev(coln), ordered = T)) %>%
  group_by(param) %>%
  summarise(Med = dmode(value),
            Lo = HPDI(value, prob = .95)[1],
            Up = HPDI(value, prob = .95)[2]) %>%
  ggplot(aes(y = Med, x = param, ymin = Lo, ymax = Up)) +
  geom_errorbar(width = .25, color = "darkred") +
  geom_point(size = 2) +
  coord_flip() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey10") +
 # scale_y_continuous(breaks = seq(-600, 600, by = 200)) +
  theme(axis.title.y = element_text(angle = 360)) +
  labs(y = "Most probable parameter value\nof predictor with 95% HPDI",
       x = "Predictor") 

#ggsave("resultsSeqLearn/gfx/Figure3.png", width = 6, height = 6)
ggsave("plots/FigureAdultsRTlogmodel.png", width = 6, height = 6)

