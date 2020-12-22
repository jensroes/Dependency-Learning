library(MASS)
library(tidyverse)
library(brms)

# Load data
children <- read_csv("Children/data/seqlearn.csv") %>%
  mutate(group = "children")

# Load data
adults <- read_csv("Adults/data/seqlearn.csv") %>%
  mutate(group = "adults")


d <- bind_rows(children, adults) %>%
  mutate(dependency = factor(dependency),
         subj = as.numeric(factor(paste(group, subj)))) %>%
  filter(acc == 1
        , rt > 100
        , rt < 10000
        , !is.na(dependency2)
        , dependency != "both"
        , block %in% c(6,7)) 

# Determine contrasts
d$COND <- paste(d$group, d$dependency, d$block, sep = "_")
d$COND <- as.factor(d$COND)
summary(d$COND)

cmat <- fractions(matrix(c(
   # Main effects
   1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,  # Group
   1, 1, 1, 1,-2,-2, 1, 1, 1, 1,-2,-2,  # Dependency
  -1,-1, 1, 1, 0, 0,-1,-1, 1, 1, 0, 0,  # Adjacency
   1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, # Block
   # Interactions
   1,-1, 1,-1,-2, 2, 1,-1, 1,-1,-2, 2, # Dep * Block
  -1, 1, 1,-1, 0, 0,-1, 1, 1,-1, 0, 0, # Adj * Block 
   1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1,-1,  # Group * Block
   1, 1, 1, 1,-2,-2,-1,-1,-1,-1, 2, 2,   # Group * Dep
  -1,-1, 1, 1, 0, 0, 1, 1,-1,-1, 0, 0,  # Group * Adj
   1,-1, 1,-1,-2, 2,-1, 1,-1, 1, 2,-2,# Dep * Block * Group
  -1,-1, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0  # Adj * Block * Group
), nrow=length(levels(d$COND)), byrow=F));cmat

rown <- levels(d$COND)
rownames(cmat) <- rown

coln <- c("Age Group", "Dependency", "Adjacency", "Block", "Dep:Block", "Adj:Block",
          "Age:Block", "Age:Dep", "Age:Adj", "Dep:Block:Age", "Adj:Block:Age")

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
          (COND + scale(trial_in_block) |p| subj) + 
          (1 | stim), 
        ndt ~ 0 + group, sigma ~ 0 + group)

get_prior(m, data = d, family = shifted_lognormal())

priors <- c(prior(normal(7, 5), class = "Intercept"),
            prior(normal(0, 2), class = "b"),
            prior(cauchy(0, 2.5), class = "sd"),
            prior(normal(0, 2), class = "sd", group = "subj"),
            prior(normal(0, 2), class = "sd", group = "stim"),
            prior(lkj(2), class = "cor"),
            prior(lkj(2), class = "cor", group = "subj"),
            prior(normal(4, 2), class = "b", dpar = "ndt"),
            prior(cauchy(0, 2.5), class = "b", dpar = "sigma"))

fit <- brm(m, 
           data = d,
           prior = priors,
           family = shifted_lognormal(),
           cores = ncores,
           chains = nchains,
           iter = iter,
           sample_prior = TRUE,
           seed = 265,
           control = list(adapt_delta = .99,
                          max_treedepth = 16))


# Save model
saveRDS(fit,
        file="Group comparison/posterior/group_model.rda",
        compress="xz")

#summary(fit)
#plot(fit, pars = "^b_") 
#pp_check(fit)
conditional_effects(fit)

# Need betas, sigmas and all re sds
names(fit$fit)

samps <- brms::posterior_samples(fit, "b_|sd") # It saves all the samples from the model.
names(samps)

write_csv(samps, "Group comparison/posterior/group_model.csv")
