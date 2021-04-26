
library(magrittr)
library(tidyverse)

setwd("~/Github/trefle")

tre <- read_csv("artifacts/trefle.csv")
tre %<>% rename(Virus = 'virus', Host = 'host')
tre %<>% mutate(Trefle = 1)

clo <- read_csv("data/clover.csv")
clo %<>% select(Virus, Host)
clo %<>% mutate(Clover = 1)

full_join(clo, tre) %>% 
  mutate(Clover = replace_na(Clover, 0)) %>%
  mutate(Trefle = replace_na(Trefle, 0)) -> df

df %<>% dplyr::mutate(Missing = as.numeric((Trefle-Clover)==1))

df %<>% group_by(Host) %>% summarize(Missing = sum(Missing), Clover = sum(Clover))

# Spatial stuff

library(sf)
library(fasterize)
library(raster)

iucn <- st_read(dsn = 'C:/Users/cjcar/Dropbox/HowManyHelminths2019',
                layer = 'TERRESTRIAL_MAMMALS')
r <- getData("worldclim",var="alt",res=5) # Make a blank raster

iucn %<>% filter(binomial %in% df$Host)
iucn %<>% left_join(df, by = c("binomial" = "Host"))

map.n <- fasterize(iucn, r, field = NULL, fun = 'count')
map <- fasterize(iucn, r, field = "Missing", fun = 'sum')
map.c <- fasterize(iucn, r, field = "Clover", fun = 'sum')
plot(map/map.n)

map.d <- ((map-cellStats(map, min))/(cellStats(map, max)-cellStats(map, min)))-((map.n-cellStats(map.n, min))/(cellStats(map.n, max)-cellStats(map.n, min)))

library(RColorBrewer)

cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50)
cols2 <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(50)

par(mar = c(0.3,2,0.3,2))
par(mfrow=c(2,2))

plot(map.n, axes = F, box = T, main = '', col = rev(cols))
maps::map('world', interior = F, add = T)
title("(A)", adj = 0.005, line = -1.1, outer = TRUE)

plot(map.c/1000, axes = F, box = T, main = '', col = rev(cols))
maps::map('world', interior = F, add = T)
title("(B)", adj = 0.513, line = -1.1, outer = TRUE)

plot(map/1000, axes = F, box = T, main = '', col = rev(cols))
maps::map('world', interior = F, add = T)
title("(C)", adj = 0.005, line = -16.3, outer = TRUE)

plot(map.d, axes = F, box = T, main = '', col = rev(cols2), zlim = c(-0.45, 0.45))
maps::map('world', interior = F, add = T)
title("(D)", adj = 0.513, line = -16.3, outer = TRUE)




#title("(B)", adj = 0.05, line = -1)