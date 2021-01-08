rm(list =ls())
library(MASS)
library(rethinking)
library(tidyverse)
library(brms)
source("./functions/functions.R")

# Load data
d <- read_csv("Children/data/seqlearn.csv")
table(d$dependency2, d$block)

#d <- d[d$trial %in% filter(d, antecedent != "none" & acc == 1)$trial,]
d <- d %>% filter(acc == 1
                  , !is.na(dependency2)
                  , dependency != "both"
)
table(d$dependency, d$block)

d$dependency <- factor(d$dependency)
d$block <- as.factor(d$block)
with(d, table(dependency, block))

# Determine contrasts
d$dependency <- gsub(pattern = "-", replacement = "", d$dependency)

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
inv.cmat <- t(ginv(cmat))
rownames(inv.cmat) <- rown
colnames(inv.cmat) <- coln
#inv.cmat

# Assign contrasts
contrasts(d$COND) <- inv.cmat

# Load model
#load( file="stanout/lognormal_impl_learning_model_adults.rda")

samps <- read_csv("Children/stanout/shifted_lognormal_impl_learning_children_posterior.csv")
samps[,c(2,ncol(samps))] <- NULL
#names(samps) <- c("Int", coln, "trial_in_block")
# Output summary
#samps <- brms::posterior_samples(fit, "b_") # It saves all the samples from the model.

names(samps) <- c("alpha[1]", paste0("beta[",1:20,"]"))

samps %<>%
  gather(Parameter, value)


# ----------------------
# Extract conditional posteriors
for(i in 1:nrow(inv.cmat)){
  assign(rownames(inv.cmat)[i], 
         with(
           samps
           , value[Parameter == 'alpha[1]'] +
             value[Parameter == 'beta[1]']*inv.cmat[i, 1] +
             value[Parameter == 'beta[2]']*inv.cmat[i, 2] +
             value[Parameter == 'beta[3]']*inv.cmat[i, 3] +
             value[Parameter == 'beta[4]']*inv.cmat[i, 4] +
             value[Parameter == 'beta[5]']*inv.cmat[i, 5] +
             value[Parameter == 'beta[6]']*inv.cmat[i, 6] +
             value[Parameter == 'beta[7]']*inv.cmat[i, 7] +
             value[Parameter == 'beta[8]']*inv.cmat[i, 8] +
             value[Parameter == 'beta[9]']*inv.cmat[i, 9] +
             value[Parameter == 'beta[10]']*inv.cmat[i, 10] +
             value[Parameter == 'beta[11]']*inv.cmat[i, 11] +
             value[Parameter == 'beta[12]']*inv.cmat[i, 12] +
             value[Parameter == 'beta[13]']*inv.cmat[i, 13] +
             value[Parameter == 'beta[14]']*inv.cmat[i, 14] +
             value[Parameter == 'beta[15]']*inv.cmat[i, 15] +
             value[Parameter == 'beta[16]']*inv.cmat[i, 16] +
             value[Parameter == 'beta[17]']*inv.cmat[i, 17] +
             value[Parameter == 'beta[18]']*inv.cmat[i, 18] +
             value[Parameter == 'beta[19]']*inv.cmat[i, 19] +
             value[Parameter == 'beta[20]']*inv.cmat[i, 20] 
         )
  )
}

conds <- rownames(inv.cmat); lenconds <- length(conds); lenconds == dim(inv.cmat)[1]
pred <- eval(parse(text = conds[1]))
for(i in 2:lenconds){
  pred <- cbind(pred, eval(parse(text = conds[i])))
}
colnames(pred) <- conds
pred %<>% as_tibble() %>%
  mutate(row = 1:n()) %>%
  gather(key = COND, value = RT, -row) %>%
  separate(COND, sep = "_", into = c("Dependency", "Block")) %>%
  mutate(RT = exp(RT)) %>%
  dplyr::select(-row)

summary(pred$RT)
write_csv(pred, "Children/stanout/post_by_cond.csv")

pred %>%
  mutate(Dependency = ifelse(Dependency == "none", "baseline", Dependency)) %>%
  mutate(Dependency = factor(Dependency, levels = c("adjacent", "nonadjacent", "baseline"), ordered = T)) %>%
  group_by(Block, Dependency) %>%
  dplyr::summarise(M = median(RT),
                   Lower = HPDI(RT, prob = .95)[1],
                   Upper = HPDI(RT, prob = .95)[2]) %>%
  ggplot(aes(x = Block, y = M, 
             color = Dependency, 
             linetype = Dependency, 
             shape = Dependency,
             ymin = Lower, ymax = Upper,
             group = interaction(Dependency))) +
  theme_minimal() +
  geom_errorbar(width = .0, position = position_dodge(.5)) +
  geom_point(size = 4, position = position_dodge(.5)) +
  geom_line(position = position_dodge(.5)) +
  ggtitle("Modelled serial reaction-time data (children)") +
  scale_colour_manual(name  = "Dependency:", values =  c("black", "black", "grey50")) +
  scale_linetype_manual(name  = "Dependency:",values = c("solid", "dashed", "dotted")) +
  scale_shape_manual("Dependency:", values = c(21,23,24)) +
  theme(legend.justification = "right",
        legend.position = "bottom",
        legend.key.width = unit(1.25, "cm")) +
  labs(y = bquote(atop("Most probable parameter value "~hat(mu),"with 95% HPDIs (in msecs)"))) +
  scale_y_continuous(breaks = seq(800, 2100, 100))

ggsave("plots/children_modelled_data_shifted_lognormal.png", width = 8, height = 5)

