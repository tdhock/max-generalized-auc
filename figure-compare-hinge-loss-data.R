library(data.table)
library(ggplot2)

pred.vec <- seq(-3, 3, by=0.25)
grid.dt <- data.table(expand.grid(
  pos=pred.vec,
  neg=pred.vec))
positive.part <- function(x)ifelse(0<x, x, 0)
sigmoid <- function(x)1/(1+exp(-x))
hinge <- function(pred, label)positive.part(1-pred*label)
hingeSet <- function(DT){
  DT[, hinge.loss := hinge(pred.real, label)]
  DT[, squared.hinge := hinge.loss^2 ]
}
err.dt <- rbind(
  data.table(
    min.log.lambda=c(-Inf, 0),
    max.log.lambda=c(0, Inf),
    fn=c(1, 0), possible.fn=1,
    fp=c(0, 0), possible.fp=0,
    label=1, obs="pos"),
  data.table(
    min.log.lambda=c(-Inf, 0),
    max.log.lambda=c(0, Inf),
    fn=c(0, 0), possible.fn=0,
    fp=c(0, 1), possible.fp=1,
    label=-1, obs="neg"))
err.dt[, labels := 1]
err.dt[, errors := fp+fn]
lab.dt <- unique(err.dt[, .(label, obs)])
loss.wide.list <- list()
for(pred.i in 1:nrow(grid.dt)){
  pred.wide <- grid.dt[pred.i]
  pred.wide[, cat(sprintf(
    "%4d / %4d %f %f\n",
    pred.i, nrow(grid.dt),
    pos, neg))]
  pred.tall <- melt(
    pred.wide,
    measure.vars=lab.dt$obs,
    value.name="pred.real",
    variable.name="obs"
  )[lab.dt, on="obs"]
  hingeSet(pred.tall)
  pred.tall[, logistic.loss := log(1+exp(-label*pred.real))]
  pred.tall[, pred.label := ifelse(0<pred.real, 1, -1)]
  pred.tall[, `01.loss` := ifelse(pred.label == label, 0, 1)]
  pred.tall[, pred.prob := sigmoid(pred.real)]
  pairwise <- dcast(
    melt(pred.tall, measure.vars=c("pred.prob","pred.real")),
    variable ~ obs)
  pairwise[, pred.real := pos-neg]
  pairwise[, label := 1]
  hingeSet(pairwise)
  pred.tall[, pred.log.lambda := pred.real]
  roc.list <- penaltyLearning::ROChange(
    err.dt, pred.tall,
    problem.vars="obs")
  err.marg <- data.table(err.dt)
  err.marg[, `:=`(
    min.log.lambda = c(-Inf, 1, -Inf, -1),
    max.log.lambda = c(1, Inf, -1, Inf)
  )]
  roc.marg <- penaltyLearning::ROChange(
    err.marg, pred.tall,
    problem.vars="obs")
  out.row <- data.table(
    pred.i, pred.wide,
    AUM=roc.list$aum,
    AUC=roc.list$auc,
    AUM.margin=roc.marg$aum,
    pairwise=dcast(
      pairwise, NULL ~ variable, value.var=c("hinge.loss", "squared.hinge")))
  for(out.col in c("01.loss","logistic.loss","hinge.loss","squared.hinge")){
    out.row[[out.col]] <- sum(pred.tall[[out.col]])
  }
  loss.wide.list[[pred.i]] <- out.row
}
loss.wide <- do.call(rbind, loss.wide.list)
data.table::fwrite(loss.wide, "figure-compare-hinge-loss-data.csv")
