library(MASS)
library(tidyverse)
library(brms)

# Load data
children <- read_csv("Children/data/seqlearn.csv") %>%
  mutate(group = "children")

# Load data
adults <- read_csv("Adults/data/seqlearn.csv") %>%
  mutate(group = "adults")

# Combine data
d <- bind_rows(children, adults) %>%
  mutate(dependency = factor(dependency),
         subj = as.numeric(factor(paste(group, subj)))) %>%
  filter(acc == 1
        , rt > 100
        , rt < 10000
        , !is.na(dependency2)
        , dependency != "both"
        , block %in% c(1,6,7)) 

table(d$block)

# Determine contrasts
d$COND <- paste(d$group, d$dependency, d$block, sep = "_")
d$COND <- factor(d$COND)
summary(d$COND)
(nrow <- length(levels(d$COND)))

dim(contrasts(d$COND))

# Main effects matrix
cmat <- fractions(matrix(c(
   # Main effects
   1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  # Group
   1, 1, 1, 1, 1, 1,-2,-2,-2, 1, 1, 1, 1, 1, 1,-2,-2,-2,  # Dependency
  -1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,  # Adjacency
   1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0,  # 1 Block 1 vs 6
   0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1  # 2 Block 6 vs 7
), nrow = nrow, byrow = F))

# Interactions
groupXdep <- cmat[,1] * cmat[,2] # group * Dep
groupXadj <- cmat[,1] * cmat[,3] # group * Adj
groupXblock1 <- cmat[,1] * cmat[,4] # group * block 1
groupXblock2 <- cmat[,1] * cmat[,5] # group * block 2
depXblock1 <- cmat[,2] * cmat[,4] # Dep * block 1
depXblock2 <- cmat[,2] * cmat[,5] # Dep * block 2
adjXblock1 <- cmat[,3] * cmat[,4] # Adj * block 1
adjXblock2 <- cmat[,3] * cmat[,5] # Adj * block 2
depXblock1Xgroup <- cmat[,1] * cmat[,2] * cmat[,4]
depXblock2Xgroup <- cmat[,1] * cmat[,2] * cmat[,5]
adjXblock1Xgroup <- cmat[,1] * cmat[,3] * cmat[,4]
adjXblock2Xgroup <- cmat[,1] * cmat[,3] * cmat[,5]

# Combine MEs and interactions
cmat <- cbind(cmat, groupXdep, groupXadj, 
              groupXblock1, groupXblock2, 
              depXblock1, depXblock2,
              adjXblock1, adjXblock2, 
              depXblock1Xgroup, depXblock2Xgroup, 
              adjXblock1Xgroup, adjXblock2Xgroup)
  
rown <- levels(d$COND)
rownames(cmat) <- rown

coln <- c("Age Group", "Dependency", "Adjacency", 
          "Block 1 vs 6", "Block 6 vs 7", 
          "Age:Dep", "Age:Adj", 
          "Age:Block 1 vs 6", "Age:Block 6 vs 7",  
          "Dep:Block 1 vs 6", "Dep:Block 6 vs 7",  
          "Adj:Block 1 vs 6", "Adj:Block 6 vs 7",  
          "Age:Dep:Block 1 vs 6", "Age:Dep:Block 6 vs 7", 
          "Age:Adj:Block 1 vs 6", "Age:Adj:Block 6 vs 7")

colnames(cmat)<- coln
cmat; colSums(cmat)

# transpose and inverse matrix to get contrast between the expected levels
inv.cmat<-fractions(t(ginv(cmat)))
rownames(inv.cmat) <- rown
colnames(inv.cmat) <- coln
inv.cmat

dim(cmat)
# Assign contrasts
dim(contrasts(d$COND))
contrasts(d$COND) <- inv.cmat
dim(contrasts(d$COND))
dim(inv.cmat)

contrasts(d$COND)

# Fit model
nchains = ncores = 3
iter = 8000

m <- bf(formula = rt ~ COND + scale(trial_in_block) + 
          (dependency * block + scale(trial_in_block) | subj) + 
          (1 | stim), ndt ~ 0 + group, sigma ~ 0 + group)

get_prior(m, data = d, family = shifted_lognormal()) %>%
  as_tibble() %>%
  filter(dpar == "sigma")

priors <- c(prior(normal(7, 5), class = "Intercept"),
            prior(normal(0, 5), class = "b"),
            prior(cauchy(0, 2.5), class = "sd"),
            prior(normal(0, 2), class = "sd", group = "subj"),
            prior(normal(0, 2), class = "sd", group = "stim"),
            prior(lkj(2), class = "cor", group = "subj"),
            prior(normal(4, 2), class = "b", dpar = "ndt"),
            prior(cauchy(0, 1), class = "b", dpar = "sigma"))

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
        file="Group comparison/posterior/group_model_v3.rda",
        compress="xz")

#summary(fit)
#plot(fit, pars = "^b_") 
#plot(fit, pars = "prior") 
#pp_check(fit)
conditional_effects(fit)

# Need betas, sigmas and all re sds
names(fit$fit)

samps <- brms::posterior_samples(fit, "b_|sd") # It saves all the samples from the model.
names(samps)

write_csv(samps, "Group comparison/posterior/group_model_v3.csv")
