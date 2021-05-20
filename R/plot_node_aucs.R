## Plot AUCs for individual nodes in the network

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ModelMetrics)

clover <- read_csv("data/clover.csv", col_types = cols(.default = "c"))

svd_preds <- read_csv("hpc/outputs/predictions.csv", 
                      col_types = cols(.default = "d", value = "l", host = "c", virus = "c"))


# ---- Add missing predictions: --------------------------------------------------------------------
# Some host-virus pairs could not be predicted at all - assume these are predicted as very unlikely
# and rank them last (using a value much lower than any predicted score)
missing_preds <- table(host = svd_preds$host, virus = svd_preds$virus) %>% 
  as_tibble() %>% 
  filter(.data$n == 0) %>% 
  select(-.data$n) %>% 
  mutate(updated = -999)

svd_preds <- bind_rows(svd_preds, missing_preds)

# ---- Host AUC: -----------------------------------------------------------------------------------
get_host_auc <- function(host_spp, obs_data = clover, pred_data = svd_preds) {
  known_positives <- obs_data %>% 
    filter(.data$Host == host_spp) %>% 
    pull(.data$Virus)
  
  label_data <- obs_data %>% 
    distinct(.data$Virus) %>% 
    mutate(label = .data$Virus %in% known_positives)
  
  preds <- pred_data %>% 
    filter(.data$host == host_spp) %>% 
    left_join(label_data, by = c("virus" = "Virus"))
  
  auc(actual = preds$label, predicted = preds$updated)
}

host_aucs <- clover %>% 
  group_by(.data$Host) %>% 
  summarise(n_viruses = n_distinct(.data$Virus)) %>% 
  group_by(.data$Host) %>% 
  mutate(auc = get_host_auc(.data$Host))

outlier_auc <- host_aucs %>% 
  filter(.data$n_viruses >= 50)

human_auc <- host_aucs %>% 
  filter(.data$Host == "Homo sapiens")

other_aucs <- host_aucs %>% 
  filter(!.data$Host %in% human_auc$Host)

p_host <- ggplot(other_aucs, aes(x = log10(n_viruses), y = auc)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey20") +
  geom_point(colour = "#4477AA") +
  geom_point(colour = "#EE6677", data = human_auc) +
  geom_text(aes(label = Host), nudge_x = 0.025, hjust = 0, size = 2, colour = "grey20", 
            fontface = "italic", data = outlier_auc) +
  scale_x_continuous(breaks = seq(0, 5, by = 0.5),
                     expand = expansion(mult = c(0.015, 0.09))) +
  labs(x = expression(log[10]*bgroup("(", textstyle("Number of known viruses"), ")")),
       y = "AUC",
       colour = NULL) +
  theme_bw()

# ---- Virus AUC: ----------------------------------------------------------------------------------
get_virus_auc <- function(virus_spp, obs_data = clover, pred_data = svd_preds) {
  known_positives <- obs_data %>% 
    filter(.data$Virus == virus_spp) %>% 
    pull(.data$Host)
  
  label_data <- obs_data %>% 
    distinct(.data$Host) %>% 
    mutate(label = .data$Host %in% known_positives)
  
  preds <- pred_data %>% 
    filter(.data$virus == virus_spp) %>% 
    left_join(label_data, by = c("host" = "Host"))
  
  auc(actual = preds$label, predicted = preds$updated)
}

virus_aucs <- clover %>% 
  group_by(.data$Virus) %>% 
  summarise(n_hosts = n_distinct(.data$Host)) %>% 
  group_by(.data$Virus) %>% 
  mutate(auc = get_virus_auc(.data$Virus))

outlier_auc <- virus_aucs %>% 
  filter(.data$n_hosts >= 90)

p_virus <- ggplot(virus_aucs, aes(x = log10(n_hosts), y = auc)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey20") +
  geom_point(colour = "#4477AA") +
  geom_text(aes(label = Virus), nudge_x = 0.025, hjust = 0, size = 2, colour = "grey20", 
            fontface = "italic", data = outlier_auc) + 
  scale_x_continuous(breaks = seq(0, 5, by = 0.5),
                     expand = expansion(mult = c(0.015, 0.12))) +
  labs(x = expression(log[10]*bgroup("(", textstyle("Number of known hosts"), ")")),
       y = "AUC",
       colour = NULL) +
  theme_bw()


# ---- Combine / save ------------------------------------------------------------------------------
p_combined <- plot_grid(p_host, p_virus, ncol = 1, labels = c("A", "B"))

ggsave2("figures/per-node-auc.pdf", p_combined, width = 7.5, height = 7.5)
