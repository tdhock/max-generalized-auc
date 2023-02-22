param.dt <- data.table::CJ(
  loss_name=c("logistic","balanced","AUM","AUM_rate"),
  seed=1:4,
  lr=10^seq(-4, 2),
  data_name=c("MNIST", "FashionMNIST"),
  batch_size=1000)
MyFun <- function(loss_name, seed, lr, data_name, batch_size){
  cmd <- paste(
    "python figure-aum-neural-networks-data.py",
    loss_name, seed, lr, data_name, batch_size)
  status <- system(cmd)
  if(status != 0){
    stop("error code ", status)
  }
}
unlink("registry",recursive=TRUE)
reg <- batchtools::makeRegistry("registry")
batchtools::batchMap(
  MyFun,
  args=param.dt,
  reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 1*24*60*60,#seconds
  memory = 8000,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=TRUE), reg=reg)
reg <- batchtools::loadRegistry("registry")
batchtools::getStatus(reg=reg)
jt <- batchtools::getJobTable(reg=reg)
jt[!is.na(error)]

"figure-aum-neural-networks-data/AUM/1/1e-06/FashionMNIST/1000/steps.csv"
(steps.csv.vec <- Sys.glob(
  "figure-aum-neural-networks-data/*/*/*/*/*/steps.csv"))
system(paste("wc -l", paste(steps.csv.vec, collapse=" ")))
unlink(grep("_count", steps.csv.vec, value=TRUE))

library(data.table)
steps.dt <- data.table(steps.csv=steps.csv.vec)[, {
  steps <- fread(
    steps.csv, colClasses=list(
      integer=c("epoch","step"),
      numeric="out_value",
      character=c("set_name","out_name")))
  meta <- nc::capture_first_vec(
    steps.csv,
    "figure-aum-neural-networks-data/",
    loss=".*?",
    "/",
    seed=".*?",
    "/",
    lr=".*?",
    "/",
    data_set=".*?",
    "/",
    batch_size=".*?",
    "/steps.csv")
  data.table(meta, steps)
}, by=steps.csv]
library(ggplot2)
steps.dt[, iteration := epoch*(1+max(step))+step]
steps.dt[, step.size := as.numeric(lr)]
one <- steps.dt[
  data_set=="MNIST" & out_name=="AUC" & set_name=="validation"]
ggplot()+
  facet_grid(step.size ~ seed, labeller=label_both)+
  geom_line(aes(
    iteration/(1+max(step)), out_value, color=loss),
    data=one)+
  scale_x_continuous("epoch")

ggplot()+
  facet_grid(step.size ~ seed, labeller=label_both)+
  geom_line(aes(
    iteration, out_value, color=loss),
    data=one[epoch < 10])
