rm(list=ls())
load("../results/step1_imputing.RData")

svd_dat <- svd(dat_impute)
zz <- svd_dat$u[,1:2]%*%diag(sqrt(svd_dat$d[1:2]))
plot(zz[,1], zz[,2], asp = T, pch = 16)

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
table(cell_type_vec)
col_idx <- as.numeric(as.factor(cell_type_vec))

png("../figure/experiment/Week31_marques_imputed_svd.png", width = 1000, height = 1400,
    res = 300, units = "px")
plot(zz[,1], zz[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "SVD on imputed data")
legend("topleft", c("Pdgfra+", "Precursors",
                        "Differentiated\nPrecursors", "Newly formed",
                        "Myelin forming", "Mature OL"), cex=0.5, fill= c(6,5,1,4,2,3))
graphics.off()
# magenta - cyan - black - blue - red - green

set.seed(10)
tsne_res <- tsne::tsne(dat_impute)
png("../figure/experiment/Week31_marques_impute_tsne.png", width = 1000, height = 1400,
    res = 300, units = "px")
plot(tsne_res[,1], tsne_res[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "T-SNE on imputed data")
graphics.off()


###################

# let's try with just normal data
load("../results/step0_screening.RData")
dat <- log2(dat+1)
set.seed(10)
tsne_res <- tsne::tsne(dat, perplexity = 80)
png("../figure/experiment/Week31_marques_raw_tsne.png", width = 1000, height = 1400,
    res = 300, units = "px")
plot(tsne_res[,1], tsne_res[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "T-SNE on raw data")
graphics.off()

svd_dat <- svd(dat)
zz <- svd_dat$u[,1:2]%*%diag(sqrt(svd_dat$d[1:2]))
png("../figure/experiment/Week31_marques_raw_svd.png", width = 1000, height = 1400,
    res = 300, units = "px")
plot(zz[,1], zz[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "SVD on raw data")
legend("bottomright", c("Pdgfra+", "Precursors",
                        "Differentiated\nPrecursors", "Newly formed",
                        "Myelin forming", "Mature OL"), cex=0.5, fill= c(6,5,1,4,2,3))
graphics.off()

#####################

rm(list=ls())
load("../results/step4_factorization_tmp.RData")
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
col_idx <- as.numeric(as.factor(cell_type_vec))

png("../figure/experiment/Week31_marques_imputed_our.png", width = 1600, height = 1400,
    res = 300, units = "px")
plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Our method on imputed data")
graphics.off()

#####################

xx <- as.numeric(dat)
yy <- as.numeric(dat_impute)

png("../figure/experiment/Week31_marques_imputation.png", width = 1600, height = 1400,
    res = 300, units = "px")
plot(xx, yy, asp = T, pch = 16, col = rgb(0,0,0,0.1),
     xlab = "Raw data (transformed)", ylab = "After imputing")
lines(c(-10000,10000), c(-10000,10000), col = "red", lwd = 2, lty = 2)
graphics.off()
