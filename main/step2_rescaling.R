set.seed(10)
load(paste0("../results/step1_naive_svd", suffix, ".RData"))

# this step (rescaling values so the maximum value is 1000) is purely to prevent underflow downstream
rescaling_factor <- 1000/max(dat)
dat_impute <- dat
dat_impute <- dat_impute * rescaling_factor
dim(dat_impute)


rm(list = c("dat"))
print(paste0(Sys.time(), ": Finished rescaling"))
source_code_info <- c(source_code_info, readLines("../main/step2_rescaling.R"))
save.image(paste0("../results/step2_rescaling", suffix, ".RData"))
print(warnings())
