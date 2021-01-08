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
          "Dep:Block 1 vs 2", "Dep:Block 6 vs 7",  
          "Adj:Block 1 vs 2", "Adj:Block 6 vs 7",  
          "Age:Dep:Block 1 vs 6", "Age:Dep:Block 6 vs 7", 
          "Age:Adj:Block 1 vs 6", "Age:Adj:Block 6 vs 7")

colnames(cmat)<- coln
cmat; colSums(cmat)

# transpose and inverse matrix to get contrast between the expected levels
inv.cmat<-fractions(t(ginv(cmat)))
rownames(inv.cmat) <- rown
colnames(inv.cmat) <- coln
inv.cmat

contrasts(d$COND)

# Load model
samps <- read_csv("Group comparison/posterior/group_model_v3.csv")
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

names(samps)[grepl("b_", names(samps))] <- c("alpha[1]", paste0("beta[",1:17,"]"))

samps %<>% gather(Parameter, value) %>% filter(grepl("alpha|beta", Parameter)) 

table(samps$Parameter)

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
             value[Parameter == 'beta[17]']*inv.cmat[i, 17] 
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
write_csv(pred_rt, "Group comparison/posterior/by_group_v3.csv")

pred_rt %>% count(Block)

sigmas <- tibble(adults = sigma_adults, children = sigma_children) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "sigma") %>%
  group_by(Group) %>%
  summarise(sigma = mean(sigma))

pred_std <- pred %>% as_tibble() %>%
  mutate(row = 1:n()) %>%
  gather(key = COND, value = RT, -row) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  group_by(Block) %>%
  mutate(row = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(diff_transfer = `7`-`6`,
         diff_learning = `6`-`1`) %>%
  dplyr::select(-`7`, -`6`, -`1`) %>%
  pivot_longer(diff_transfer:diff_learning) %>%
  left_join(sigmas, by = "Group") %>%
  mutate(value = value / sigma) %>%
  mutate(name = gsub("diff", "sdt", name)) 
  
pred_std_sum <- pred_std %>% group_by(Group, Dependency, name, sigma) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2]) %>%
  unite("Parameter", c(name,Group,Dependency), sep = "_") %>% dplyr::select(-sigma) %>%
  arrange(Parameter)

pred_std_es <- pred_std %>%
  group_by(Group, Dependency, name) %>%
  mutate(row =  1:n()) %>%
  pivot_wider(-sigma, names_from = c("name", "Group", "Dependency"), values_from = value) %>%
  dplyr::select(-row) %>%
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% left_join(pred_std_sum) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE) %>%
  mutate(Parameter = gsub("sdt_", "", Parameter))

pred_rt_es <- pred_rt %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(diff_transfer = `Block: 7`- `Block: 6`,
         diff_learning = `Block: 6` - `Block: 1`) %>%
  dplyr::select(-`Block: 1`:-`Block: 7`) %>%
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

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_v3.csv")


pred_std <- pred %>% as_tibble() %>%
  mutate(adults_dependency_1 = (adults_adjacent_1 + adults_nonadjacent_1)/2,
         adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_1 = (children_adjacent_1 + children_nonadjacent_1)/2,
         children_dependency_6 = (children_adjacent_6 + children_nonadjacent_6)/2,
         children_dependency_7 = (children_adjacent_7 + children_nonadjacent_7)/2) %>%
  dplyr::select(-contains("adjacent"), -contains("nonadjacent")) %>%
  mutate(row = 1:n()) %>%
  gather(key = COND, value = RT, -row) %>%
  separate(COND, sep = "_", into = c("Group", "Dependency", "Block")) %>%
  dplyr::select(-row) %>%
  group_by(Block) %>%
  mutate(row = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(diff_transfer = `7`-`6`,
         diff_learning = `6`-`1`) %>%
  dplyr::select(-`7`, -`6`, -`1`) %>%
  pivot_longer(diff_transfer:diff_learning) %>%
  left_join(sigmas, by = "Group") %>%
  mutate(value = value / sigma) %>%
  mutate(name = gsub("diff_", "", name)) 

pred_std_sum <- pred_std %>% group_by(Group, Dependency, name, sigma) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2]) %>%
  unite("Parameter", c(name,Group,Dependency), sep = "_") %>% dplyr::select(-sigma) %>%
  arrange(Parameter) 


pred_std_es <- pred_std %>%
  group_by(Group, Dependency, name) %>%
  mutate(row =  1:n()) %>%
  pivot_wider(-sigma, names_from = c("name", "Group", "Dependency"), values_from = value) %>%
  dplyr::select(-row) %>% 
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% left_join(pred_std_sum) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Parameter, Med, Lo, Up, ROPE)

pred_rt_es <- pred %>% as_tibble() %>%
  mutate(adults_dependency_1 = (adults_adjacent_1 + adults_nonadjacent_1)/2,
         adults_dependency_6 = (adults_adjacent_6 + adults_nonadjacent_6)/2,
         adults_dependency_7 = (adults_adjacent_7 + adults_nonadjacent_7)/2,
         children_dependency_1 = (children_adjacent_1 + children_nonadjacent_1)/2,
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
         learning = `6` - `1`) %>%
  pivot_longer(transfer:learning) %>%
  group_by(Group, Dependency, name) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  unite("Parameter", c(name, Group, Dependency), sep = "_") %>%
  left_join(pred_std_es, "Parameter") %>%
  separate(Parameter, into = c("Effect", "Group", "Dependency"), sep = "_")

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_dependency_v3.csv")


# Age X dependency
pred_std <- pred %>% as_tibble() %>%
  transmute(adults_dependency = (adults_adjacent_1 + adults_nonadjacent_1 + adults_adjacent_6 + adults_nonadjacent_6 + adults_adjacent_7 + adults_nonadjacent_7)/6,
            adults_none = (adults_none_1 + adults_none_6 + adults_adjacent_7)/3,
            children_dependency = (children_adjacent_1 + children_nonadjacent_1 + children_adjacent_6 + children_nonadjacent_6 + children_adjacent_7 + children_nonadjacent_7)/6,
            children_none = (children_none_1 + children_none_6 + children_adjacent_7)/3) %>%
  pivot_longer(everything()) %>%
  separate(name, sep = "_", into = c("Group", "Dependency")) %>%
  group_by(Dependency) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Dependency, values_from = value) %>%
  mutate(dependency = `dependency`-`none`) %>% dplyr::select(Group, dependency) %>%
  left_join(sigmas, by = "Group") %>% mutate(value = dependency / sigma) 

pred_std_sum <- pred_std %>% group_by(Group, sigma) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2]) %>% dplyr::select(-sigma) 

pred_std_es <- pred_std %>%
  group_by(Group) %>%
  mutate(row =  1:n()) %>%
  pivot_wider(-sigma:-dependency, names_from = c("Group"), values_from = value) %>%
  dplyr::select(-row) %>% 
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% rename(Group = Parameter) %>%
  left_join(pred_std_sum) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  dplyr::select(Group, Med, Lo, Up, ROPE)

pred_rt_es <- pred %>% as_tibble() %>%
  transmute(adults_dependency = (adults_adjacent_1 + adults_nonadjacent_1 + adults_adjacent_6 + adults_nonadjacent_6 + adults_adjacent_7 + adults_nonadjacent_7)/6,
            adults_none = (adults_none_1 + adults_none_6 + adults_adjacent_7)/3,
            children_dependency = (children_adjacent_1 + children_nonadjacent_1 + children_adjacent_6 + children_nonadjacent_6 + children_adjacent_7 + children_nonadjacent_7)/6,
            children_none = (children_none_1 + children_none_6 + children_adjacent_7)/3) %>%
  pivot_longer(everything()) %>%
  separate(name, sep = "_", into = c("Group", "Dependency")) %>%
  mutate(RT = exp(value)) %>%
  group_by(Dependency) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(-value, names_from = Dependency, values_from = RT) %>%
  mutate(dependency = `dependency`-`none`) %>% dplyr::select(Group, dependency) %>%
  rename(value = dependency) %>%
  group_by(Group) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  left_join(pred_std_es, "Group") 

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_dependency_agegroup_v3.csv")



# Age X learning
pred_std <- pred %>% as_tibble() %>%
  transmute(adults_block1 = (adults_adjacent_1 + adults_nonadjacent_1 + adults_none_1)/3,
            adults_block6 = (adults_adjacent_6 + adults_nonadjacent_6 + adults_none_6)/3,
            adults_block7 = (adults_adjacent_7 + adults_nonadjacent_7 + adults_none_7)/3,
            children_block1 = (children_adjacent_1 + children_nonadjacent_1 + children_none_1)/3,            
            children_block6 = (children_adjacent_6 + children_nonadjacent_6 + children_none_6)/3,
            children_block7 = (children_adjacent_7 + children_nonadjacent_7 + children_none_7)/3) %>%
  pivot_longer(everything()) %>%
  separate(name, sep = "_", into = c("Group", "Block")) %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = value) %>%
  mutate(learning = block6 - block1,
         transfer = block7 - block6) %>% dplyr::select(Group, learning, transfer) %>%
  left_join(sigmas, by = "Group") %>% mutate(learning = learning / sigma,
                                             transfer = transfer / sigma) %>% dplyr::select(-sigma)

pred_std_sum <- pred_std %>% pivot_longer(learning:transfer) %>%
  group_by(Group, name) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2])

pred_std_es <- pred_std %>% pivot_longer(learning:transfer) %>%
  group_by(Group, name) %>%
  mutate(row =  1:n()) %>%
  pivot_wider(names_from = c("Group", "name"), values_from = value) %>%
  dplyr::select(-row) %>% 
  rope(range = c(-.1,.1), ci = .89) %>%
  as_tibble() %>% separate(Parameter, into = c("Group", "name")) %>%
  left_join(pred_std_sum) %>%
  mutate(ROPE = ROPE_Percentage*100) %>%
  rename(Effect = name) %>%
  dplyr::select(Group, Effect,  Med, Lo, Up, ROPE)

pred_rt_es <- pred %>% as_tibble() %>%
  transmute(adults_block1 = (adults_adjacent_1 + adults_nonadjacent_1 + adults_none_1)/3,
            adults_block6 = (adults_adjacent_6 + adults_nonadjacent_6 + adults_none_6)/3,
            adults_block7 = (adults_adjacent_7 + adults_nonadjacent_7 + adults_none_7)/3,
            children_block1 = (children_adjacent_1 + children_nonadjacent_1 + children_none_1)/3,            
            children_block6 = (children_adjacent_6 + children_nonadjacent_6 + children_none_6)/3,
            children_block7 = (children_adjacent_7 + children_nonadjacent_7 + children_none_7)/3) %>%
  pivot_longer(everything()) %>%
  separate(name, sep = "_", into = c("Group", "Block")) %>%
  mutate(RT = exp(value)) %>% dplyr::select(-value) %>%
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Block, values_from = RT) %>%
  mutate(learning = block6 - block1,
         transfer = block7 - block6) %>% dplyr::select(Group, learning, transfer) %>%
  pivot_longer(learning:transfer) %>%
  group_by(Group, name) %>%
  summarise(Med = dmode(value),
            Lo = rethinking::HPDI(value, prob = .95)[1],
            Up = rethinking::HPDI(value, prob = .95)[2],
            P = mean(value < 0)) %>%
  rename(Effect = name) %>%
  left_join(pred_std_es, c("Group", "Effect")) 

write_csv(pred_rt_es, "Group comparison/posterior/by_group_effectsize_learning_agegroup_v3.csv")

