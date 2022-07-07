library(data.table)
library(ggplot2)

pred.vec <- seq(-3, 3, by=0.25)
grid.dt <- data.table(expand.grid(
  pos=pred.vec,
  neg=pred.vec))
positive.part <- function(x)ifelse(0<x, x, 0)
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
  pred.tall <- melt(
    pred.wide,
    measure.vars=lab.dt$obs,
    value.name="pred.log.lambda",
    variable.name="obs"
  )[lab.dt, on="obs"]
  pred.tall[, hinge.loss := positive.part(1-label*pred.log.lambda)]
  pred.tall[, logistic.loss := log(1+exp(-label*pred.log.lambda))]
  pred.tall[, pred.label := ifelse(0<pred.log.lambda, 1, -1)]
  pred.tall[, `01.loss` := ifelse(pred.label == label, 0, 1)]
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
  loss.wide.list[[pred.i]] <- data.table(
    pred.i, pred.wide,
    AUM=roc.list$aum,
    AUC=roc.list$auc,
    AUM.margin=roc.marg$aum,
    `01.loss`=sum(pred.tall[["01.loss"]]),
    logistic.loss=sum(pred.tall$logistic.loss),
    hinge.loss=sum(pred.tall$hinge.loss))
}
loss.wide <- do.call(rbind, loss.wide.list)
data.table::fwrite(loss.wide, "figure-compare-hinge-loss-data.csv")
