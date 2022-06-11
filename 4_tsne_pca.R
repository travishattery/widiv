# t-SNE script
setwd("C:/<WORKING DIRECTORY>")
library(Rtsne)
input <- read.csv("2_widiv_3env_plotavg_log10_wthr50.csv")
tsne_out <- Rtsne(as.matrix(input[,11:57]),pca=T, pca_center=T, pca_scale=T)
write.csv(tsne_out$Y,"tsne_results.csv", row.names = FALSE)

# PCA script
#install.packages("factoextra")
#library(factoextra)
setwd("C:/<WORKING DIRECTORY>")
input <- read.csv("2_widiv_3env_plotavg_wthr50.csv")
pca <- prcomp(input[,c(11:57)], center = TRUE,scale. = TRUE)

summary(pca)
loadings <- pca$rotation
