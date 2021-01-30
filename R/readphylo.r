library(ape)
upham <- read.nexus("data/upham_tree.nex")
D <- cophenetic(upham)

df <- data.frame(D["Homo_sapiens",])
colnames(df) <- c("phylodist")
write.csv(df, "artifacts/phylo_distance_to_human.csv")
