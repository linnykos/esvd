load(paste0("../results/step5_clustering", suffix, ".RData"))

idx <- which.min(quality_vec)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
tmp <- cbind(pred_mat[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)
plot(-pred_mat[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2), main = "Our embedding",
     xlab = "Predicted value", ylab = "Observed and masked value")
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)



#########
#under construction

d <- 3
dat <- res_our$u_mat[,1:d]

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.2
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))

combn_mat <- combn(3,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(dat[,i1])
  ylim <- range(dat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/main/marques_latent_our_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Our latent dimension ", i1),
         ylab = paste0("Our latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()

}
