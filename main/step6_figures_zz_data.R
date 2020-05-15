var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step6_figures", suffix, ".RData"))


########################

# make the plot displaying how evenly expressed cells are

mean_vec <- log(apply(dat_impute, 2, mean))
sd_vec <- log(apply(dat_impute, 2, sd))

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

#########################################

# prepare things for vioplot plot

zz <- stats::prcomp(dat_impute, center = T, scale. = F)
vec <- zz$rotation[,1]
vec[vec <= 0] <- 0
vec <- (vec - min(vec))/(max(vec) - min(vec))
vec <- vec/sum(vec)
xx <- dat_impute %*% vec

custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))
tmp_df <- data.frame(val = xx, type = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))})[cluster_labels])

##########################################

png("../../esvd_results/figure/main/data_overview.png", height = 1200, width = 2200, res = 300, units = "px")
graphics::par(mfrow = c(1,2), mar = c(5,6,4,2))
graphics::plot(NA, asp = T,
     xlim = range(mean_vec),
     ylim = range(sd_vec),
     xlab = "Mean of expression (Log)",
     ylab = "Standard deviation\nof expression (Log)", main = "Standard deviation verses\nmean per gene",
     cex.lab = 1.25)
graphics::lines(c(-1e2,1e2), c(-1e2,1e2), col = "red", lwd = 2, lty = 2)
graphics::lines(rep(0,2), c(-1e2,1e2), col = "red", lwd = 2)
graphics::lines(c(-1e2,1e2), rep(0,2), col = "red", lwd = 2)
for(i in 1:len){
  idx <- which(idx_vec == i)
  graphics::points(mean_vec[idx], sd_vec[idx], pch = 16, col = col_palette[i])
}


graphics::legend("topleft", c("Evenly expressed", "Unevenly expressed"),
       fill=c( rgb(240/255, 228/255, 66/255), "black"), cex = 0.9)

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6],
                 col = sapply(1:6, function(x){unique(col_info_svd$col_code[which(round(col_info_svd$order) == x)])}),
                 pchMed = 21,
                 colMed = "black", colMed2 = "white",
                 xlab = "", names = rep("", 6))
graphics::title(ylab = "Weighted expression based\non 1st principal component",
      main = "Average gene expression\nper cell type", cex.lab = 1.25)
graphics::text(1:6, par("usr")[3]-10,
     srt = -45, xpd = TRUE,
     labels = c("Pdgfra+", "OPC", "COP", "NFOL", "MFOL", "MOL"), cex=1)
grDevices::graphics.off()
