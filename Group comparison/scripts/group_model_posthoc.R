library(MASS)
library(tidyverse)
library(brms)
library(rethinking)
source("./functions/functions.R")

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

# Load model
samps <- read_csv("Group comparison/posterior/group_model_posterior.csv")
names(samps)

sigma_adults <- exp(samps$b_sigma_groupadults)
sigma_children <- exp(samps$b_sigma_groupchildren)

samps$b_ndt_Intercept <- NULL
samps$b_scaletrial_in_block <- NULL

names(samps)[grepl("b_", names(samps))] <- c("alpha[1]", paste0("beta[",1:11,"]"))

samps %<>% gather(Parameter, value)

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
             value[Parameter == 'beta[11]']*inv.cmat[i, 11] 
         )
  )
}

conds <- rownames(inv.cmat); lenconds <- length(conds); lenconds == dim(inv.cmat)[1]
pred <- eval(parse(text = conds[1]))
for(i in 2:lenconds){
  pred <- cbind(pred, eval(parse(text = conds[i])))
}
colnames(pred) <- conds

pred_rt <- pred %>%
  as_tibble() %>%
  gather(key = COND, value = RT) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  mutate(RT = exp(RT))

summary(pred_rt$RT)
write_csv(pred_rt, "Group comparison/posterior/by_group.csv")


pred_std <- pred %>% as_tibble() %>%
  mutate(row = 1:n(),
         sigma_children = sigma_children,
         sigma_adults = sigma_adults) %>%
  gather(key = COND, value = RT, -row,-sigma_children,-sigma_adults) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  dplyr::select(-row) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(value = `7`- `6`,
         std_value = ifelse(Group == "adults", value / sigma_adults,
                    ifelse(Group == "children", value / sigma_children, NA))) %>%
  dplyr::select(-starts_with("sigma"), -`6`, -`7`, -value)
summary(pred_std)

pred_std_ests <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_value),
            Lo = rethinking::HPDI(std_value, prob = .95)[1],
            Up = rethinking::HPDI(std_value, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_")

pred_std_es <- pred_std %>%
  group_by(Group, Dependency) %>%
  mutate(id =  1:n()) %>%
  pivot_wider(names_from = c("Group", "Dependency"), values_from = std_value) %>%
  dplyr::select(adults_adjacent:children_none) %>%
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% left_join(pred_std_ests) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE)

pred_rt_es <- pred_rt %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(value = `7`- `6`) %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  left_join(pred_std_es, "Parameter") %>%
  separate(Parameter, into = c("Group", "Dep"), sep = "_")


write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize.csv")


pred_std <- pred %>% as_tibble() %>%
  mutate(adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_6 = (children_adjacent_6 + children_nonadjacent_6)/2,
         children_dependency_7 = (children_adjacent_7 + children_nonadjacent_7)/2) %>%
  dplyr::select(-contains("adjacent"), -contains("nonadjacent")) %>%
  mutate(row = 1:n(),
         sigma_children = sigma_children,
         sigma_adults = sigma_adults) %>%
  gather(key = COND, value = RT, -row,-sigma_children,-sigma_adults) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  dplyr::select(-row) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  group_by(Group, Dependency) %>%
  mutate(value = `7`- `6`,
         std_value = ifelse(Group == "adults", value / sigma_adults,
                            ifelse(Group == "children", value / sigma_children, NA))) %>%
  dplyr::select(-starts_with("sigma"), -`6`, -`7`, -value)


pred_std_ests <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_value),
            Lo = rethinking::HPDI(std_value, prob = .95)[1],
            Up = rethinking::HPDI(std_value, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_")

pred_std_es <- pred_std %>%
  group_by(Group, Dependency) %>%
  mutate(id =  1:n()) %>%
  pivot_wider(names_from = c("Group", "Dependency"), values_from = std_value) %>%
  dplyr::select(adults_none:children_dependency) %>% 
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% left_join(pred_std_ests) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE)

pred_rt_es <- pred %>%
  as_tibble() %>%
  mutate(adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_6 = (children_adjacent_6 + children_nonadjacent_6)/2,
         children_dependency_7 = (children_adjacent_7 + children_nonadjacent_7)/2) %>%
  dplyr::select(-contains("adjacent"), -contains("nonadjacent")) %>%
  gather(key = COND, value = RT) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  mutate(RT = exp(RT)) %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(value = `7`- `6`) %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  left_join(pred_std_es, "Parameter") %>%
  separate(Parameter, into = c("Group", "Dependency"), sep = "_")

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_dependency.csv")
