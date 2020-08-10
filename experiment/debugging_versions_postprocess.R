rm(list=ls())
load("../experiment/debugging_0025.RData")

dat_alt <- dat
dat2_alt <- dat2
nat_mat_alt <- nat_mat
u_mat_alt <- u_mat
u_mat2_alt <- u_mat2
u_mat3_alt <- u_mat3
v_mat_alt <- v_mat
v_mat2_alt <- v_mat2
v_mat3_alt <- v_mat3
init_alt <- init

var_list <- ls()
var_list <- var_list[-grep("*alt", var_list)]
rm(list = var_list)

load("../experiment/debugging_0037.RData")

dat[1:5,1:5]
dat_alt[1:5,1:5]

dat2[1:5,1:5]
dat2_alt[1:5,1:5]

nat_mat[1:5,1:5]
nat_mat_alt[1:5,1:5]

head(u_mat)
head(u_mat_alt)

head(u_mat2)
head(u_mat2_alt)

head(u_mat3)
head(u_mat3_alt)

dim(dat)
dim(dat_alt)
sum(abs(dat - dat_alt))
sum(abs(dat2 - dat2_alt))
sum(abs(nat_mat - nat_mat_alt))
sum(abs(u_mat - u_mat_alt))
sum(abs(u_mat2 - u_mat2_alt))
sum(abs(u_mat3 - u_mat3_alt))
sum(abs(v_mat - v_mat_alt))
sum(abs(v_mat2 - v_mat2_alt))
sum(abs(v_mat3 - v_mat3_alt))
sum(abs(init$u_mat - init_alt$u_mat))
sum(abs(init$v_mat - init_alt$v_mat))
