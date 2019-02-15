rm(list=ls())
load("../../../SOUPR/data/zeisel.rda")
dat <- zeisel$counts
dim(dat)
dat[1:5,1:5]

dat2 <- read.table("../../raw_data/GSE75330_Marques_et_al_mol_counts2.tab")
dim(dat2)
colnames(dat2) <- dat2[1,]
rownames(dat2) <- dat2[,1]
dat2 <- dat2[-1,-1]
dat2 <- t(dat2)
dat2 <- as.data.frame(dat2)
for(i in 1:nrow(dat2)){
  dat2[,i] <- as.numeric(dat2[,i])
}
dat2[1:5,1:5]

#############################
# first see if the gene names agree
length(which(colnames(dat) %in% colnames(dat2)))/ncol(dat) #94%...
length(which(colnames(dat2) %in% colnames(dat)))/ncol(dat2) #66%...

# see which cells overlap
name_vec <- as.vector(sapply(rownames(dat2), function(x){
  zz <- strsplit(x, split = "-")[[1]]
  paste0(paste0(zz[2], zz[3]), "_", zz[4])
}))
length(which(rownames(dat) %in% name_vec))/nrow(dat) #11%...
length(which(name_vec %in% rownames(dat)))/length(name_vec) #6%...

# does the intersection even overlap?
dat_b <- dat[which(rownames(dat) %in% name_vec), which(colnames(dat) %in% colnames(dat2))]
dat2_b <- dat2[which(name_vec %in% rownames(dat)), which(colnames(dat2) %in% colnames(dat))]
dat_b <- dat_b[order(rownames(dat_b)), order(colnames(dat_b))]
dat2_b <- as.matrix(dat2_b[order(rownames(dat2_b)), order(colnames(dat2_b))])
zz <- as.numeric(dat2_b)
zz <- matrix(zz, ncol = ncol(dat2_b), nrow = nrow(dat2_b))
colnames(zz) <- colnames(dat2_b)
rownames(zz) <- rownames(dat2_b)
dat2_b <- zz
sum(abs(dat2_b - dat_b)) #it's a little different hm
