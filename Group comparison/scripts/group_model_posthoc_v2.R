library(MASS)
library(tidyverse)
library(brms)
library(magrittr)
library(rethinking)
library(bayestestR)
source("./functions/functions.R")

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
         , block %in% c(1,2,6,7)) 

table(d$block)

# Determine contrasts
d$COND <- paste(d$group, d$dependency, d$block, sep = "_")
d$COND <- factor(d$COND)
summary(d$COND)
(nrow <- length(levels(d$COND)))

dim(contrasts(d$COND))

cmat <- fractions(matrix(c(
  # Main effects
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  # Group
  1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,-1,-1,-1,-1, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,-1,-1,-1,-1,  # Dependency
  -1,-1,-1,-1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1, 1, 1, 1, 1, 0, 0, 0, 0,  # Adjacency
  1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0,  # 1 Block 1 vs 2
  0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,  # 2 Block 6 vs 7
  0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0   # 3 Block 2 vs 6
), nrow = nrow, byrow = F))


groupXdep <- cmat[,1] * cmat[,2] # group * Dep
groupXadj <- cmat[,1] * cmat[,3] # group * Adj
groupXblock1 <- cmat[,1] * cmat[,4] # group * block 1
groupXblock2 <- cmat[,1] * cmat[,5] # group * block 2
groupXblock3 <- cmat[,1] * cmat[,6] # group * block 3
depXblock1 <- cmat[,2] * cmat[,4] # Dep * block 1
depXblock2 <- cmat[,2] * cmat[,5] # Dep * block 2
depXblock3 <- cmat[,2] * cmat[,6] # Dep * block 3
adjXblock1 <- cmat[,3] * cmat[,4] # Adj * block 1
adjXblock2 <- cmat[,3] * cmat[,5] # Adj * block 2
adjXblock3 <- cmat[,3] * cmat[,6] # Adj * block 3
depXblock1Xgroup <- cmat[,1] * cmat[,2] * cmat[,4]
depXblock2Xgroup <- cmat[,1] * cmat[,2] * cmat[,5]
depXblock3Xgroup <- cmat[,1] * cmat[,2] * cmat[,6]
adjXblock1Xgroup <- cmat[,1] * cmat[,3] * cmat[,4]
adjXblock2Xgroup <- cmat[,1] * cmat[,3] * cmat[,5]
adjXblock3Xgroup <- cmat[,1] * cmat[,3] * cmat[,6]


cmat <- cbind(cmat, groupXdep, groupXadj, groupXblock1, groupXblock2, groupXblock3,
              depXblock1, depXblock2, depXblock3, adjXblock1, adjXblock2, adjXblock3,
              depXblock1Xgroup, depXblock2Xgroup, depXblock3Xgroup, adjXblock1Xgroup, adjXblock2Xgroup, adjXblock3Xgroup)

rown <- levels(d$COND)
rownames(cmat) <- rown

coln <- c("Age Group", "Dependency", "Adjacency", 
          "Block 1 vs 2", "Block 6 vs 7", "Block 2 vs 6", 
          "Age:Dep", "Age:Adj", 
          "Age:Block 1 vs 2", "Age:Block 6 vs 7", "Age:Block 2 vs 6", 
          "Dep:Block 1 vs 2", "Dep:Block 6 vs 7", "Dep:Block 2 vs 6", 
          "Adj:Block 1 vs 2", "Adj:Block 6 vs 7", "Adj:Block 2 vs 6", 
          "Age:Dep:Block 1 vs 2", "Age:Dep:Block 6 vs 7", "Age:Dep:Block 2 vs 6", 
          "Age:Adj:Block 1 vs 2", "Age:Adj:Block 6 vs 7", "Age:Adj:Block 2 vs 6")

colnames(cmat)<- coln
cmat; colSums(cmat)

# transpose and inverse matrix to get contrast between the expected levels
inv.cmat<-fractions(t(ginv(cmat)))
rownames(inv.cmat) <- rown
colnames(inv.cmat) <- coln
inv.cmat

# Load model
samps <- read_csv("Group comparison/posterior/group_model_v2.csv")
names(samps)

sigma_adults <- exp(samps$b_sigma_groupadults)
sigma_children <- exp(samps$b_sigma_groupchildren)

samps$b_ndt_groupadults <- NULL
samps$b_ndt_groupchildren <- NULL
samps$b_scaletrial_in_block <- NULL
samps$b_sigma_groupadults <- NULL
samps$b_sigma_groupchildren <- NULL
samps$prior_b_ndt <- NULL
samps$prior_b_sigma <- NULL

names(samps)[grepl("b_", names(samps))] <- c("alpha[1]", paste0("beta[",1:23,"]"))

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
             value[Parameter == 'beta[11]']*inv.cmat[i, 11] +
             value[Parameter == 'beta[12]']*inv.cmat[i, 12] +
             value[Parameter == 'beta[13]']*inv.cmat[i, 13] +
             value[Parameter == 'beta[14]']*inv.cmat[i, 14] +
             value[Parameter == 'beta[15]']*inv.cmat[i, 15] +
             value[Parameter == 'beta[16]']*inv.cmat[i, 16] +
             value[Parameter == 'beta[17]']*inv.cmat[i, 17] +
             value[Parameter == 'beta[18]']*inv.cmat[i, 18] +
             value[Parameter == 'beta[19]']*inv.cmat[i, 19] +
             value[Parameter == 'beta[20]']*inv.cmat[i, 20] +
             value[Parameter == 'beta[21]']*inv.cmat[i, 21] + 
             value[Parameter == 'beta[22]']*inv.cmat[i, 22] +
             value[Parameter == 'beta[23]']*inv.cmat[i, 23] 
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
  mutate(RT = exp(RT),
         Block = paste0("Block: ", Block))

summary(pred_rt$RT)
write_csv(pred_rt, "Group comparison/posterior/by_group_v2.csv")

pred_rt %>% count(Block)

pred_std <- pred %>% as_tibble() %>%
  mutate(row = 1:n(),
         sigma_children = sigma_children,
         sigma_adults = sigma_adults) %>%
  gather(key = COND, value = RT, -row,-sigma_children,-sigma_adults) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  dplyr::select(-row) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(diff_transfer = `7`-`6`,
         diff_learning = `2`-`1`,
         std_transfer = ifelse(Group == "adults", diff_transfer / sigma_adults,
                     ifelse(Group == "children", diff_transfer / sigma_children, NA)),
         std_learning = ifelse(Group == "adults", diff_learning / sigma_adults,
                               ifelse(Group == "children", diff_learning / sigma_children, NA))) %>%
  dplyr::select(-starts_with("sigma"), -`7`, -`6`, -`1`, -`2`, -diff_transfer, -diff_learning)
summary(pred_std)

pred_std

pred_std_ests_learn <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_learning),
            Lo = rethinking::HPDI(std_learning, prob = .95)[1],
            Up = rethinking::HPDI(std_learning, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  mutate(Parameter = paste0("std_learning_", Parameter))

pred_std_ests_transf <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_transfer),
            Lo = rethinking::HPDI(std_transfer, prob = .95)[1],
            Up = rethinking::HPDI(std_transfer, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  mutate(Parameter = paste0("std_transfer_", Parameter))

pred_std_ests <- bind_rows(pred_std_ests_learn, pred_std_ests_transf)

pred_std_es <- pred_std %>%
  group_by(Group, Dependency) %>%
  mutate(id =  1:n()) %>%
  pivot_wider(names_from = c("Group", "Dependency"), values_from = c(std_learning, std_transfer)) %>%
  dplyr::select(-id) %>%
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% left_join(pred_std_ests) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE) %>%
  mutate(Parameter = gsub("std_", "", Parameter))

pred_rt_es <- pred_rt %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(diff_transfer = `Block: 7`- `Block: 6`,
         diff_learning = `Block: 2` - `Block: 1`) %>%
  select(-`Block: 1`:-`Block: 7`) %>%
  pivot_longer(cols = starts_with("diff_")) %>%
  group_by(Group, Dependency, name) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  mutate(name = gsub("diff_", "", name)) %>%
  unite("Parameter", name:Parameter, sep = "_") %>%
  left_join(pred_std_es, "Parameter") %>%
  separate(Parameter, into = c("Parameter", "Group", "Dep"), sep = "_")

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_v2.csv")


pred_std <- pred %>% as_tibble() %>%
  mutate(adults_dependency_1 = (adults_adjacent_1 + adults_nonadjacent_1)/2,
         adults_dependency_2 = (adults_adjacent_2 + adults_nonadjacent_2)/2,
         adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_1 = (children_adjacent_1 + children_nonadjacent_1)/2,
         children_dependency_2 = (children_adjacent_2 + children_nonadjacent_2)/2,
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
  mutate(transfer = `7`-`6`,
         learning = `2`-`1`,
         std_transfer = ifelse(Group == "adults", transfer / sigma_adults,
                            ifelse(Group == "children", transfer / sigma_children, NA)),
         std_learning = ifelse(Group == "adults", learning / sigma_adults,
                               ifelse(Group == "children", learning / sigma_children, NA))) %>%
  dplyr::select(-starts_with("sigma"), -`7`, -`2`, -`6`, -`1`, -transfer, -learning)


pred_std_ests_trans <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_transfer),
            Lo = rethinking::HPDI(std_transfer, prob = .95)[1],
            Up = rethinking::HPDI(std_transfer, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  mutate(effect = "transfer")

pred_std_ests_learn <- pred_std %>%
  group_by(Group, Dependency) %>%
  summarise(Med = dmode(std_learning),
            Lo = rethinking::HPDI(std_learning, prob = .95)[1],
            Up = rethinking::HPDI(std_learning, prob = .95)[2]) %>%
  unite("Parameter", Group:Dependency, sep = "_") %>%
  mutate(effect = "learning")

pred_std_ests <- bind_rows(pred_std_ests_trans, pred_std_ests_learn) %>%
  mutate(Parameter = paste(effect, Parameter, sep = "_"))

pred_std_es <- pred_std %>%
  group_by(Group, Dependency) %>%
  mutate(id =  1:n()) %>%
  pivot_wider(names_from = c("Group", "Dependency"), values_from = c(std_transfer, std_learning)) %>%
  dplyr::select(-id) %>% 
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% mutate(Parameter = gsub("std_", "", Parameter)) %>%
  left_join(pred_std_ests) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE)

pred_rt_es <- pred %>% as_tibble() %>%
  mutate(adults_dependency_1 = (adults_adjacent_1 + adults_nonadjacent_1)/2,
         adults_dependency_2 = (adults_adjacent_2 + adults_nonadjacent_2)/2,
         adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_1 = (children_adjacent_1 + children_nonadjacent_1)/2,
         children_dependency_2 = (children_adjacent_2 + children_nonadjacent_2)/2,
         children_dependency_6 = (children_adjacent_6 + children_nonadjacent_6)/2,
         children_dependency_7 = (children_adjacent_7 + children_nonadjacent_7)/2) %>%
  dplyr::select(-contains("adjacent"), -contains("nonadjacent")) %>%
  gather(key = COND, value = RT) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  mutate(RT = exp(RT)) %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(transfer = `7`- `6`,
         learning = `2` - `1`,
         std_transfer = ifelse(Group == "adults", transfer / sigma_adults,
                               ifelse(Group == "children", transfer / sigma_children, NA)),
         std_learning = ifelse(Group == "adults", learning / sigma_adults,
                               ifelse(Group == "children", learning / sigma_children, NA))) %>%
  pivot_longer(std_transfer:std_learning) %>%
  group_by(Group, Dependency, name) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  mutate(name = gsub("std_", "", name)) %>%
  unite("Parameter", c(name, Group, Dependency), sep = "_") %>%
  left_join(pred_std_es, "Parameter") %>%
  separate(Parameter, into = c("Effect", "Group", "Dependency"), sep = "_")

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_dependency_v2.csv")
