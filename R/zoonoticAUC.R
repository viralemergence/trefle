
library(tidyverse)
library(PresenceAbsence)

pred <- read_csv("~/Github/trefle/hpc/outputs/predictions.csv")
pred <- read_csv("~/Github/trefle/hpc/outputs/predictions.csv")

pred %>% filter(host == 'Homo sapiens') %>%
  select(virus, value, updated) %>% filter(value == TRUE)

# Oh no

clo <- read_csv("~/Github/clover/clover/Clover_v1.0_NBCIreconciled_20201218.csv")

clo %>% filter(Host == 'Homo sapiens') %>% pull(Virus) -> zoonoses

pred %>% filter(host == 'Homo sapiens') %>%
  select(virus, value, updated) %>%
  mutate(value = (virus %in% zoonoses)) %>%
  mutate(updated = scales::rescale(updated, to = c(0, 1))) -> pred

pred 

pred %>% mutate(value = as.numeric(value)) %>% 
  mutate(virus = row_number()) -> pred

auc.roc.plot(pred)
