rm(list=ls())
load("../results/step5_clustering.RData")

mean_vec <- log(apply(dat_impute, 2, mean))
sd_vec <- log(apply(dat_impute, 2, sd))

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)

p_vec <- sapply(1:ncol(dat_impute), function(i){
  tmp_df <- data.frame(val = dat_impute[,i], type = cluster_labels)
  fit <- stats::aov(val ~ type, data = tmp_df)
  sum_fit <- summary(fit)
  sum_fit[[1]][[5]][1]
})

len <- 5
col_palette <- grDevices::colorRampPalette(c("black", rgb(240/255, 228/255, 66/255)))(len)
idx_vec <- sapply(p_vec, function(x){
  which.min(abs(x - seq(0, 1, length.out = len)))
})

plot(NA, asp = T,
     xlim = range(mean_vec),
     ylim = range(sd_vec),
     xlab = "Mean of expression (Log)",
     ylab = "Standard deviation of expression (Log)", main = "Standard deviation verses mean per gene")
lines(c(-1e2,1e2), c(-1e2,1e2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e2,1e2), col = "red", lwd = 2)
lines(c(-1e2,1e2), rep(0,2), col = "red", lwd = 2)
for(i in 1:len){
  idx <- which(idx_vec == i)
  points(mean_vec[idx], sd_vec[idx], pch = 16, col = col_palette[i])
}


legend("topleft", c("More evenly expressed across cell types", "Less evenly expressed across cell types"),
       fill=c( rgb(240/255, 228/255, 66/255), "black"))

############

idx <- which.min(p_vec)
# idx <- which.max(p_vec)
tmp_df <- data.frame(val = dat_impute[,idx], type = cluster_labels)
boxplot(val ~ type, data = tmp_df)
