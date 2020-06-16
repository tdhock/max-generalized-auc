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
  pred.tall[, pred.label := ifelse(0<pred.log.lambda, 1, -1)]
  pred.tall[, `01.loss` := ifelse(pred.label == label, 0, 1)]
  roc.list <- penaltyLearning::ROChange(
    err.dt, pred.tall,
    problem.vars="obs")
  loss.wide.list[[pred.i]] <- data.table(
    pred.i, pred.wide,
    aum=roc.list$aum,
    auc=roc.list$auc,
    `01.loss`=sum(pred.tall[["01.loss"]]),
    hinge.loss=sum(pred.tall$hinge.loss))
}
loss.wide <- do.call(rbind, loss.wide.list)

loss.wide[, neg.auc := -auc]
loss.tall <- melt(
  loss.wide,
  measure.vars = c("auc", "aum", "hinge.loss"),
  variable.name="loss.name",
  value.name="loss.value")
normalize <- function(x)(x-min(x))/(max(x)-min(x))
loss.tall[, loss.norm := normalize(loss.value), by=loss.name]
rect.dt <- data.table(
  xmin=0, xmax=Inf,
  ymin=-Inf, ymax=0)
gg <- ggplot()+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
    data=loss.tall)+
  ## geom_contour(aes(
  ##   pos, neg, z=loss.norm),
  ##   data=loss.tall)+
  geom_rect(aes(
    xmin=xmin, xmax=xmax,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="black",
    data=rect.dt)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  scale_fill_gradient(
    "Relative\nvalues",
    low="white",
    high="red")+
  coord_equal()+
  xlab("Real-valued prediction for positive label")+
  ylab("Real-valued prediction\nfor negative label")
png("figure-compare-hinge-loss.png", width=6, height=2, res=100, units="in")
print(gg)
dev.off()
