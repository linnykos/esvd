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

png("../figure/experiment/Week30_marques_all.png", width = 1000, height = 1200,
    res = 300, units = "px")
plot(zz[,1], zz[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2")
legend("topleft", c("Precursors", "Newly formed",
                    "Myelin forming", "Mature OL"), cex=0.6,
       bty="n", fill= c("black", "blue", "red", "green"))
graphics.off()

## we can also do a t-sne
tsne_res <- tsne::tsne(dat_impute)
png("../figure/experiment/Week30_marques_all_tsne.png", width = 1000, height = 1200,
    res = 300, units = "px")
plot(tsne_res[,1], tsne_res[,2], asp = T, pch = 16, col = col_idx,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2")
legend("topleft", c("Precursors", "Newly formed",
                    "Myelin forming", "Mature OL"), cex=0.6,
       bty="n", fill= c("black", "blue", "red", "green"))
graphics.off()
