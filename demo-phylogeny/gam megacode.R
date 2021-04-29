
# 3x2 panels
# columns: sharing across all species, sharing with humans
# rows: binary sharing, jaccard, number of viruses
# on each plot, same axis, two lines - one pre one post gam fit

library(ape)
library(ggplot2)
library(magrittr)
library(patchwork)  
library(tidyverse)

# 1. Turn CLOVER and TREFLE into three (each) pairwise sharing matrices by-host

setwd("~/Github/trefle")

tre <- read_csv("artifacts/trefle.csv")

library(igraph)

### Some hijacked code from Albersnet
m <- table(tre[,1:2])
M <- as.matrix(m)
bipgraph <- graph.incidence(M, weighted = NULL)
Hostgraph <- bipartite.projection(bipgraph)$proj2
HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(tre$host)
Remove <- which(rowSums(HostAdj)==diag(HostAdj))
HostAdj <- HostAdj[-Remove,-Remove]
###

tre.df <- (data.frame(Host = row.names(HostAdj), HostAdj) %>% reshape2::melt())
tre.df %>% rename(Counts.Post = value) %>% mutate(Sharing.Post = as.numeric(Counts.Post > 0)) -> tre.df

clo <- read_csv("data/clover.csv")
clo %<>% select(Virus, Host)

### Some hijacked code from Albersnet
m <- table(clo[,1:2])
M <- as.matrix(m)
bipgraph <- graph.incidence(M, weighted = NULL)
Hostgraph <- bipartite.projection(bipgraph)$proj2
HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(tre$host)
Remove <- which(rowSums(HostAdj)==diag(HostAdj))
HostAdj <- HostAdj[-Remove,-Remove]
###

clo.df <- (data.frame(Host = row.names(HostAdj), HostAdj) %>% reshape2::melt())
clo.df %>% rename(Counts.Pre = value) %>% mutate(Sharing.Pre = as.numeric(Counts.Pre > 0)) -> clo.df

sharing.df <- left_join(tre.df, clo.df)




# 2. Generate GAMS against total pairwise matrix in Upham 

unames <- read_csv("~/Github/clover/clover/phylogenies/mammal_phylo_translations.csv")
sharing.df %>% select(Host) %>% unique() -> hosts
hosts <- merge(hosts, unames, by="Host")

## load phylogeny
tree <- read.nexus("data/upham_tree.nex")

## fix tips
tree$tip.label <- gsub("_"," ",tree$tip.label)

## check difference
setdiff(hosts$Host_Upham,tree$tip.label)

## trim tree
tree <- keep.tip(tree,hosts$Host_Upham)

## get distance from Homo sapiens
cdist <- cophenetic.phylo(tree)

pairwise.dist <- (data.frame(Host = row.names(cdist), cdist) %>% reshape2::melt())

all.data <- left_join(sharing.df, pairwise.dist)
all.data %<>% rename(PhyloDist = value)

# 3. Generate GAMs against log distance from humans in Upham

all.data %<>% as_tibble()

all.data %>%
  pivot_longer(
    cols = ends_with(c("Pre","post")),
    names_to = "Variable",
    values_to = "Value",
    values_drop_na = FALSE
  ) %>%
  separate(Variable, into = c("Variable","Source"), sep = "\\.") %>%
  pivot_wider(
    names_from = "Variable",
    values_from = Value
  ) -> all.data
  
all.data %>% 
  filter(PhyloDist > 0) %>% # kill self-similarity 
  filter(log(PhyloDist) < 5.5) -> all.data #kill monotremes

all.data %<>% rename(Dataset = "Source")
all.data %<>% mutate(Dataset = recode(Dataset, !!!c("Post" = "Imputed network", "Pre" = "Raw network")))

all.data %>% 
  ggplot(aes(x = PhyloDist, y = Counts, group = Dataset, color = Dataset)) +  
  geom_smooth(method = 'glm',  
              method.args = list(family = poisson)) + 
  xlab("Phylogenetic distance (pairwise)") +
  ylab("Number of viruses shared") + 
  theme_bw() + 
  scale_color_manual(values = c("#d1495b","#00798c"))  -> upper.right

all.data %>% 
  ggplot(aes(x = PhyloDist, y = Sharing, group = Dataset, color = Dataset)) + 
  geom_smooth(method = 'glm',  
              method.args = list(family = binomial(link = 'logit'))) + 
  xlab("Phylogenetic distance (pairwise)") +
  ylab("Viral sharing") + 
  theme_bw() + 
  scale_color_manual(values = c("#d1495b","#00798c")) -> upper.left

all.data %>% 
  filter(Host == "Homo sapiens") %>%
  ggplot(aes(x = PhyloDist, y = Counts, group = Dataset, color = Dataset)) + 
  geom_smooth(method = 'glm',  
              method.args = list(family = poisson)) + 
  xlab("Phylogenetic distance from humans") +
  ylab("Number of viruses shared") + 
  theme_bw() + 
  scale_color_manual(values = c("#d1495b","#00798c"))  -> bottom.right

all.data %>% 
  filter(Host == "Homo sapiens") %>%
  ggplot(aes(x = PhyloDist, y = Sharing, group = Dataset, color = Dataset)) + 
  #geom_point() + 
  geom_smooth(method = 'glm',  
              method.args = list(family = binomial(link = 'logit'))) + 
  xlab("Phylogenetic distance from humans") +
  ylab("Viral sharing") + 
  theme_bw() + 
  scale_color_manual(values = c("#d1495b","#00798c")) -> bottom.left

(upper.left + upper.right) / (bottom.left + bottom.right) + plot_layout(guides = 'collect') &
  theme(legend.position='bottom')

######

library(rsq)



m0 <- glm(Sharing ~ PhyloDist, family = binomial(link = 'logit'), data = (all.data %>% filter(Dataset == 'Raw network')))
s0 <- summary(m0); s0 
rsq(m0, adj = TRUE)

m1 <- glm(Sharing ~ PhyloDist, family = binomial(link = 'logit'), data = (all.data %>% filter(Dataset == 'Imputed network')))
s1 <- summary(m1); s1
rsq(m1, adj = TRUE)


m0 <- glm(Sharing ~ PhyloDist, family = binomial(link = 'logit'), data = (all.data %>% filter(Dataset == 'Raw network', Host == 'Homo sapiens')))
s0 <- summary(m0); s0 
rsq(m0, adj = TRUE)

m1 <- glm(Sharing ~ PhyloDist, family = binomial(link = 'logit'), data = (all.data %>% filter(Dataset == 'Imputed network', Host == 'Homo sapiens')))
s1 <- summary(m1); s1
rsq(m1, adj = TRUE)










m0 <- glm(Counts ~ PhyloDist, family = binomial(link = 'logit'), data = (all.data %>% filter(Dataset == 'Raw network')))
s0 <- summary(m0); s0 
rsq(m0, adj = TRUE)

m1 <- glm(Counts ~ PhyloDist, family = "poisson", data = (all.data %>% filter(Dataset == 'Imputed network')))
s1 <- summary(m1); s1
rsq(m1, adj = TRUE)


m0 <- glm(Counts ~ PhyloDist, family = "poisson", data = (all.data %>% filter(Dataset == 'Raw network', Host == 'Homo sapiens')))
s0 <- summary(m0); s0 
rsq(m0, adj = TRUE)

m1 <- glm(Counts ~ PhyloDist, family = "poisson", data = (all.data %>% filter(Dataset == 'Imputed network', Host == 'Homo sapiens')))
s1 <- summary(m1); s1
rsq(m1, adj = TRUE)








