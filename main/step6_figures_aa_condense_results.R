load(paste0("../results/step5_trajectory", suffix, ".RData"))

# from https://www.r-bloggers.com/size-of-each-object-in-rs-workspace/
for (thing in ls()) {
  message(thing)
  print(object.size(get(thing)), units='auto')
}

rm(list=c("nat_mat_list_list", "dat_count"))

save.image(paste0("../results/step6_figures_condensed", suffix, ".RData"))

