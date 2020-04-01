rm(list=ls())
load("../results/step5_clustering.RData")

# prepare things for sd-mean plot
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

# prepare things for vioplot plot
zz <- prcomp(dat_impute, center = T, scale. = F)
vec <- zz$rotation[,1]
vec[vec <= 0] <- 0
vec <- (vec - min(vec))/(max(vec) - min(vec))
vec <- vec/sum(vec)
xx <- dat_impute %*% vec

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(210/255, 198/255, 36/255, alpha)) #darker yellow
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)

custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

col_vec <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

tmp_df <- data.frame(val = xx, type = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))})[cluster_labels])


png("../../esvd_results/figure/experiment/Revision_writeup6_svd_preview.png", height = 1200, width = 2200, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(5,6,4,2))
plot(NA, asp = T,
     xlim = range(mean_vec),
     ylim = range(sd_vec),
     xlab = "Mean of expression (Log)",
     ylab = "Standard deviation\nof expression (Log)", main = "Standard deviation verses\nmean per gene",
     cex.lab = 1.25)
lines(c(-1e2,1e2), c(-1e2,1e2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e2,1e2), col = "red", lwd = 2)
lines(c(-1e2,1e2), rep(0,2), col = "red", lwd = 2)
for(i in 1:len){
  idx <- which(idx_vec == i)
  points(mean_vec[idx], sd_vec[idx], pch = 16, col = col_palette[i])
}


legend("topleft", c("Evenly expressed", "Unevenly expressed"),
       fill=c( rgb(240/255, 228/255, 66/255), "black"), cex = 0.9)

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6],
                 col = sapply(1:6, function(x){unique(col_info$col_code[which(col_info$level == x)])}),
                 pchMed = 21,
                 colMed = "black", colMed2 = "white",
                 xlab = "", names = rep("", 6))
title(ylab = "Weighted expression based\non 1st principal component",
      main = "Average gene expression\nper cell type", cex.lab = 1.25)
text(1:6, par("usr")[3]-10,
     srt = -45, xpd = TRUE,
     labels = c("Pdgfra+", "OPC", "COP", "NFOL", "MFOL", "MOL"), cex=1)
graphics.off()

#######################################


