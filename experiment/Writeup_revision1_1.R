rm(list=ls())
library(xtable)

labels_vec <- c("../../raw_data/Lingxue_attachment/Cleaned/Baron_human1/Baron_human1_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human2/Baron_human2_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human3/Baron_human3_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human4/Baron_human4_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse1/Baron_mouse1_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse2/Baron_mouse2_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Darmanis/Darmanis_label.csv")
dat_vec <- c("../../raw_data/Lingxue_attachment/Cleaned/Baron_human1/Baron_human1_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human2/Baron_human2_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human3/Baron_human3_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human4/Baron_human4_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse1/Baron_mouse1_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse2/Baron_mouse2_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Darmanis/Darmanis_expr.csv")

i <- 7
labels <- read.csv(labels_vec[i])
dat <- read.csv(dat_vec[i])
dat <- dat[,-1]
dim(dat)

zz <- apply(dat, 2, function(x){sum(x != 0)})
length(which(zz >= floor(nrow(dat)/100)))
length(which(zz >= floor(nrow(dat)/100)))/length(zz)

head(labels)
zz <- table(labels$cluster)
tmp <- data.frame(Cluster = names(zz), Number = as.integer(zz))
xtable::xtable(tmp)
