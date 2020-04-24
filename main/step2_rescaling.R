set.seed(10)
load(paste0("../results/step1_naive_svd", suffix, ".RData"))


dat_impute <- dat
dat_impute <- dat_impute * 1000/max(dat_impute)
dim(dat_impute)

###########

# png("../figure/main/data_impute.png", height = 2400, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# eSVD:::.plot_singlecell(dat_impute)
# graphics.off()

rm(list = c("dat"))
print(paste0(Sys.time(), ": Finished rescaling"))
save.image(paste0("../results/step1_rescaling", suffix, ".RData"))
