set.seed(10)
load(paste0("../results/step1_naive_svd", suffix, ".RData"))

# this step (rescaling values so the maximum value is 1000) is purely to prevent underflow downstream
rescaling_factor <- 1000/max(dat)
dat_impute <- dat
dat_impute <- dat_impute * rescaling_factor
dim(dat_impute)


rm(list = c("dat"))
print(paste0(Sys.time(), ": Finished rescaling"))
save.image(paste0("../results/step2_rescaling", suffix, ".RData"))
print(warnings())


###########

# png("../figure/main/data_impute.png", height = 2400, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# eSVD:::.plot_singlecell(dat_impute)
# graphics.off()
