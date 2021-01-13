source("packages.R")

d <- function(min.log.lambda, fp, fn){
  data.table(min.log.lambda, fp, fn)
}
profile <- function(..., possible.fp, possible.fn, errors, labels){
  dt <- do.call(rbind, list(...))
  if(missing(possible.fp))possible.fp <- max(dt$fp)
  if(missing(possible.fn))possible.fn <- max(dt$fn)
  errors <- dt[, fp+fn]
  if(missing(labels))labels <- max(errors)
  dt[, data.table(
    min.log.lambda,
    max.log.lambda=c(min.log.lambda[-1], Inf),
    fp, fn, errors, possible.fp, possible.fn, labels)]
}
profile.list <- list(
  less=profile(
    d(-Inf, 0, 1),
    d(0, 1, 0)))
pred.dt <- data.table(problem=1, pred.log.lambda=0)
roc.dt.list <- list()
auc.dt.list <- list()
for(profile.i in seq_along(profile.list)){
  p <- data.table(profile.list[[profile.i]], problem=1)
  roc.list <- penaltyLearning::ROChange(p, pred.dt, problem.vars="problem")
  model <- names(profile.list)[[profile.i]]
  roc.dt.list[[profile.i]] <- data.table(model, roc.list$roc)
  auc.dt.list[[profile.i]] <- with(roc.list, data.table(model, auc, aum))
}
roc.dt <- do.call(rbind, roc.dt.list)
auc.dt <- do.call(rbind, auc.dt.list)
fp.fn.dt <- data.table::melt(roc.dt, measure.vars=c("fp", "fn", "min.fp.fn"))
err.sizes <- c(
  fp=3,
  fn=2,
  min.fp.fn=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  min.fp.fn="black")
ggplot()+
  facet_grid(model ~ ., labeller=label_both)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_rect(aes(
    xmin=min.thresh, xmax=max.thresh,
    ymin=0, ymax=value),
    color="grey",
    fill="grey",
    data=fp.fn.dt[variable=="min.fp.fn"])+
  geom_segment(aes(
    min.thresh, value,
    color=variable, size=variable,
    xend=max.thresh, yend=value),
    data=fp.fn.dt)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)+
  scale_x_continuous("prediction threshold")

ggplot()+
  geom_path(aes(
    FPR, TPR, color=model, size=model, group=model),
    data=roc.dt)

