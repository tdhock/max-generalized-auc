source("packages.R")

before.dt <- data.table(
  tp=0,
  fp=0,
  possible.tp=1,
  possible.fp=1)
rep.dt <- data.table(
  tp=c(1, 1, 0, 0),
  fp=c(0, 1, 1, 0),
  possible.tp=1,
  possible.fp=1)
after.dt <- data.table(
  tp=c(1, 1),
  fp=c(0, 1),
  possible.tp=1,
  possible.fp=1)

rep.list <- replicate(1, rep.dt, simplify=FALSE)
several.dt <- do.call(rbind, rep.list)
segs.dt <- rbind(before.dt, several.dt, after.dt)[.N:1]
n.breaks <- nrow(segs.dt)-1L
break.vec <- 1:n.breaks
segs.dt[, min.log.lambda := c(-Inf, break.vec)]
segs.dt[, max.log.lambda := c(break.vec, Inf)]
print(segs.dt)
segs.dt[, problem := 1]
segs.dt[, fn := possible.tp-tp]
segs.dt[, possible.fn := possible.tp]
segs.dt[, errors := fp+fn]
segs.dt[, labels := 2]
pred.dt <- data.table(pred.log.lambda=1.5, problem=1)
(L <- penaltyLearning::ROChange(segs.dt, pred.dt, "problem"))

ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=L$roc)

segs.tall <- melt(
  segs.dt,
  measure.vars=c("fp", "tp"))

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ .)+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value),
    data=segs.tall)
