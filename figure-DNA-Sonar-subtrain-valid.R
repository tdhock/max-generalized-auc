f.csv.gz <- "figure-DNA-Sonar-subtrain-valid-data.csv.gz"
if(!file.exists(f.csv.gz)){
  u <- paste0("http://ml.nau.edu/data/", f.csv.gz)
  download.file(u, f)
}
loss.dt <- data.table::fread("figure-DNA-Sonar-subtrain-valid-data.csv.gz")

(selected.dt <- loss.dt[
  set.name=="validation",
  .SD[which.max(auc), .(iteration, step.size)],
  by=.(data.name, validation.fold, loss.name)])
#=> Some step sizes selected include 100 and 0.01, which means seq(2,
#-2) is not broad enough, so we should try -3, 3 next time.
