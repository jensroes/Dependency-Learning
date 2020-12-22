library(tidyverse)
library(brms)
library(rethinking)
source("./functions/functions.R")

ps <- read_csv("Group comparison/posterior/by_group.csv")

# Difference block 6 to block 7

ps %>% 
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = "Block", values_from = "RT") %>%
  mutate(RT = `7` - `6`) %>% dplyr::select(RT, Group, Dependency) %>%
  group_by(Group, Dependency) %>% mutate(id = 1:n()) %>%
  pivot_wider(names_from = Group:Dependency, values_from = RT) %>% dplyr::select(-id) %>%
  rope(c(-10, 10))


ps %>% 
  group_by(Block) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = "Block", values_from = "RT") %>%
  mutate(RT = `7` - `6`) %>%
  group_by(Group, Dependency) %>%
  dplyr::summarise(M = median(RT),
                   Lower = HPDI(RT, prob = .95)[1],
                   Upper = HPDI(RT, prob = .95)[2]) %>% ungroup() 
  
