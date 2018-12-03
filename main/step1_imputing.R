set.seed(10)
load("../results/step0_screening.RData")

res_hvg <- descend::findHVG(res_descend, threshold = 12)
length(res_hvg$HVG.genes)

idx1 <- sort(unlist(apply(res_list[[5]]$v, 2, function(x){which(x != 0)})))
idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
idx <- sort(unique(c(idx1, idx2)))

dat <- dat[,idx]
dim(dat)

impute_res <- SAVER::saver(t(dat), ncores = 10)
dat_impute <- t(impute_res$estimate)

mean_val <- mean(dat_impute)
dat_impute <- t(apply(dat_impute, 1, function(x){x/sum(x)}))
dat_impute <- dat_impute * mean_val/mean(dat_impute)

rm(list = c("idx1", "idx2"))
print(paste0(Sys.time(), ": Finished imputing"))
save.image("../results/step1_imputing.RData")
