rm(list=ls())
load("../results/step4_factorization_100.RData")

png("../figure/experiment/Week34_marques_imputed.png", height = 1000, width = 1800, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(t(dat_impute))
graphics.off()

png("../figure/experiment/Week34_marques_raw.png", height = 1000, width = 1800, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(t(dat))
graphics.off()

length(which(dat != 0))
length(which(dat_impute != 0))

length(which(dat != 0))/prod(dim(dat))
length(which(dat_impute != 0))/prod(dim(dat))

png("../figure/experiment/Week34_marques_hist.png", height = 1000, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
max_val <- 40
hist(dat[intersect(which(dat != 0), which(dat < max_val))], breaks = 20, col = "gray",
     ylim = c(0, 1.25e5), xlab = "Values", ylab = "Number of instances", main = "Before imputation")
hist(dat_impute[intersect(which(dat_impute != 0), which(dat_impute < max_val))], breaks = 20, col = "gray",
     ylim = c(0, 1.25e5), xlab = "Values", ylab = "Number of instances", main = "After imputation")
graphics.off()

############################

rm(list=ls())
load("../results/step4_factorization.RData")
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.3
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))
i1 <- 1; i2 <- 2
xlim <- range(res_our$u_mat[,i1])
ylim <- range(res_our$u_mat[,i2])
order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
              "Myelin-forming (5)", "Mature (6)")

png("../figure/experiment/Week34_marques_latent_12.png", height = 1500, width = 2000, res = 300, units = "px")
par(mfrow = c(2,3), mar = c(4,4,4,0.5))
for(i in 1:6){
  idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
  plot(res_our$u_mat[-idx,i1], res_our$u_mat[-idx,i2], asp = T, pch = 16,
       col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
       main = name_vec[i], xlab = paste0("Latent dimension ", i1),
       ylab = paste0("Latent dimension ", i2))

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

  points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
         col = rgb(1,1,1), cex = 1.5)

  points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
         col = col_idx[idx], cex = 1.5)
}
graphics.off()

#####

num_cell <- length(unique(cell_type_vec))
i1 <- 2; i2 <- 3
xlim <- range(res_our$u_mat[,i1])
ylim <- range(res_our$u_mat[,i2])

png("../figure/experiment/Week34_marques_latent_23.png", height = 1500, width = 2000, res = 300, units = "px")
par(mfrow = c(2,3), mar = c(4,4,4,0.5))
for(i in 1:6){
  idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
  plot(res_our$u_mat[-idx,i1], res_our$u_mat[-idx,i2], asp = T, pch = 16,
       col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
       main = name_vec[i], xlab = paste0("Latent dimension ", i1),
       ylab = paste0("Latent dimension ", i2))

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

  points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
         col = rgb(1,1,1), cex = 1.5)

  points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
         col = col_idx[idx], cex = 1.5)
}
graphics.off()

