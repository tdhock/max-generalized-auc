library(data.table)
library(ggplot2)

data.vec <- c(
  "CIFAR10_N=5623",
  "FashionMNIST_N=10000",
  "MNIST_N=18032",
  "STL10_N=1778")
data.table(pre=data.vec)[, data.table(
  seed_csv=Sys.glob(paste0("data_Classif_scaled/",pre,"_seed*.csv"))
)[, fread(seed_csv,nrow=1),by=seed_csv]
]#check if weights are different in different seed files.
param.dt <- data.table(pre=data.vec)[, CJ(
  seed_csv=Sys.glob(paste0("data_Classif_scaled/",pre,"_seed*.csv")),
  objective=c("min.aum","max.auc"),
  set.obj=c("subtrain","validation")
), by=pre]
if(FALSE){
  seed_csv <- "data_Classif_scaled/FashionMNIST_N=10000_seed=4.csv"
  objective <- "min.aum"
  set.obj <- "subtrain"
}
run_one <- function(seed_csv, objective, set.obj, ...){
  library(data.table)
  seed.dt <- fread(seed_csv)
  data_csv <- sub("_seed=[0-9]", "", seed_csv)
  data.dt <- fread(data_csv)
  feature.list <- list()
  diffs.list <- list()
  for(set.name in c("subtrain","validation")){
    set.data <- data.dt[set.name==set]
    feature.list[[set.name]] <- as.matrix(set.data[,-(1:2)])
    diffs.list[[set.name]] <- aum::aum_diffs_binary(set.data$y)
  }
  fit <- aum:::aum_linear_model(
    feature.list,
    diffs.list,
    improvement.thresh=1e-3,
    maxIterations=objective,
    initial.weight.fun=function(...)seed.dt$weight,
    line.search.set=set.obj)
  out_csv <- sub(
    ".csv",
    paste0(
      "_objective=",objective,
      "_set=",set.obj,
      ".csv"),
    sub("scaled","line_search",seed_csv))
  fwrite(fit$loss,out_csv)
}

reg.dir <- "data_Classif_line_search"
unlink(reg.dir, recursive=TRUE)
reg <- batchtools::makeRegistry(reg.dir)
batchtools::batchMap(run_one, args=param.dt)
if(FALSE){
  batchtools::testJob(1,reg=reg)
}
(job.table <- batchtools::getJobTable(reg=reg))
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 4*60*60,#seconds
  memory = 4000,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=TRUE), reg=reg)
reg <- batchtools::loadRegistry(reg.dir)
batchtools::getStatus(reg=reg)

result.old <- if(file.exists("data_Classif_batchtools.csv")){
  fread("data_Classif_batchtools.csv")
}else{
  data.table()
}

f <- function(vname,...)list("_",nc::field(vname,"=",".*?",...))
result.new <- nc::capture_first_glob(
  "data_Classif_constant/*.csv",
  #"data_Classif_constant/STL10_N=1778_seed=4_lr=1.0_loss=SquaredHinge.csv"
  "/",
  data.name=".*?",
  f("N",as.integer),
  f("seed",as.integer),
  f("lr",as.numeric),
  f("loss"),
  "[.]csv")

result.dt <- rbind(result.old, result.new)
fwrite(result.dt, "data_Classif_batchtools.csv")
unlink("data_Classif_constant/*.csv")

best.valid <- result.dt[
  set_name=="validation",
  .SD[which.max(auc)],
  keyby=.(data.name,N,loss,seed)]
result.dt[step_number==0 & lr==1 & set_name=="subtrain"][order(data.name,loss,seed), .(data.name,loss,seed,auc)]##these should be different.
select.dt <- best.valid[, .(data.name, N, seed, loss, lr)]
select.result <- result.dt[select.dt, on=names(select.dt)]

gg <- ggplot()+
  geom_line(aes(
    step_number, auc,
    color=loss),
    data=select.result[set_name=="validation"])+
  ## geom_line(aes(
  ##   step_number, auc,
  ##   linetype=set_name,
  ##   color=loss),
  ##   data=select.result)+
  geom_point(aes(
    step_number, auc,
    color=loss),
    fill="black",
    shape=21,
    data=best.valid)+
  scale_x_log10()+
  coord_cartesian(
    ##xlim=c(0,1000),
    ylim=c(0.5,1))+
  facet_grid(data.name+N~seed, labeller=label_both)
png("data_Classif_batchtools.png", width=14, height=8, units="in", res=200)
print(gg)
dev.off()
