load("../results/step0_screening.RData")

idx <- which(colnames(dat) %in% res_hvg$HVG.genes)
