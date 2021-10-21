source("packages.R")

folds.dt <- fread("../feature-learning-benchmark/labeled_problems_folds.csv")
addMeta <- function(dt){
  dt[, set.name := sub("/.*", "", prob.dir)]
  dt[, problem := sub(".*/", "", prob.dir)]
  dt[folds.dt, on=list(set.name, problem)]
}
errors.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_errors.csv"))
possible.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_possible_errors.csv"))

errors.dt[, min.log.lambda := min.log.penalty]
errors.dt[, max.log.lambda := max.log.penalty]
test.fold.targets <- penaltyLearning::targetIntervals(
  errors.dt, "prob.dir")
errors.dt[, min.lambda := exp(min.log.lambda)]
errors.dt[, example := prob.dir]
finite.targets <- test.fold.targets[
  is.finite(min.log.lambda) | is.finite(max.log.lambda)
][order(prob.dir)]
finite.targets[, pred.in.interval := data.table::fcase(
  min.log.lambda == -Inf, max.log.lambda-1,
  max.log.lambda == Inf, min.log.lambda+1,
  -Inf < min.log.lambda & max.log.lambda < Inf, (min.log.lambda+max.log.lambda)/2)]
finite.targets[is.na(pred.in.interval)]
set.seed(1)
finite.targets[, pred.rnorm := rnorm(.N)]
pred.names <- finite.targets[, prob.dir]
diff.dt <- aum::aum_diffs_penalty(errors.dt[order(prob.dir)], pred.names)

## if(!file.exists("signal.list.annotation.sets.RData")){
##   download.file("https://rcdata.nau.edu/genomic-ml/cbio/neuroblastoma/signal.list.annotation.sets.RData", "signal.list.annotation.sets.RData")
## }
## (objs <- load("signal.list.annotation.sets.RData"))

squared.hinge.fun.list <- list(
  loss=function(x, e=1)ifelse(x<e,(x-e)^2,0),
  grad=function(x,e=1)ifelse(x<e,2*(x-e),0))
max.N <- length(pred.names)
N.vec <- as.integer(10^seq(1, log10(max.N), l=10))
ex.counts <- table(diff.dt[["example"]])
pred.type.vec <- grep("pred", names(finite.targets), value=TRUE)
timing.dt.list <- list()
for(N in N.vec){
  diff.N <- diff.dt[example < N]
  ex.N <- ex.counts[1:N]
  targets.N <- finite.targets[1:N]
  for(pred.type in pred.type.vec){
    pred.vec <- targets.N[[pred.type]]
    timing.df <- microbenchmark::microbenchmark(sort={
      sort(rep(pred.vec, ex.N)-diff.N[["pred"]])
    }, aum={
      aum::aum(diff.N, pred.vec)
    }, squared.hinge.each.example={
      arg.mat <- cbind(
        min=pred.vec-targets.N[["min.log.lambda"]],
        max=targets.N[["max.log.lambda"]]-pred.vec)
      result.mat.list <- list()
      for(fun.name in names(squared.hinge.fun.list)){
        fun <- squared.hinge.fun.list[[fun.name]]
        result.mat.list[[fun.name]] <- fun(arg.mat)
      }
      with(result.mat.list, list(
        loss=mean(loss[,"min"]+loss[,"max"]),
        grad=grad[,"min"]-grad[,"max"]))
    }, times=10)
    timing.dt.list[[paste(N, pred.type)]] <- with(timing.df, data.table(
      N, pred.type, seconds=time/1e9, algorithm=expr))
  }
}
(timing.dt <- do.call(rbind, timing.dt.list))
data.table::fwrite(timing.dt, "figure-aum-grad-speed-data.csv")

