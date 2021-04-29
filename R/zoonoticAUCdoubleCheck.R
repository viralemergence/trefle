
setwd("~/Github/trefle")
library(tidyverse)
library(PresenceAbsence)

tre <- read_csv("hpc/outputs/predictions.csv")
clo <- read_csv("data/clover.csv")

clo %>%
  select(Virus, Host) %>% 
  group_by(Virus) %>%
  summarize(zoonotic = max(Host %in% 'Homo sapiens')) -> zoon

tre %>% filter(host == 'Homo sapiens') %>%
  select(virus, updated) %>%
  rename(Virus = 'virus') %>%
  mutate(updated = ((updated - min(updated))/(max(updated)-min(updated)))) -> P

df <- left_join(zoon, P)

df %>% ggplot(aes(x=factor(zoonotic), y=updated)) +
  geom_violin()
PresenceAbsence::auc.roc.plot(df)

tre %>% filter(host == 'Homo sapiens') %>%
  select(virus, initial, updated) %>%
  rename(Virus = 'virus') -> P

df <- left_join(zoon, P)
