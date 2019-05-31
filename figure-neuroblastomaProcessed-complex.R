source("packages.R")

data(neuroblastomaProcessed, package="penaltyLearning")

counts <- neuroblastomaProcessed$errors[, {
  diff.tab <- table(factor(diff(errors), c("-1", "0", "1")))
  L <- as.list(diff.tab)
  size <- max.log.lambda-min.log.lambda
  for(fun.name in c("min", "max")){
    fun <- get(fun.name)
    L[[paste0(fun.name, ".size")]] <- min(size[errors==fun(errors)])
  }
  L$mean.size <- with(L, (min.size+max.size)/2)
  L
}, by=list(profile.id, chromosome)]
two.changes <- counts[1 < `-1` | 1 < `1`]
two.changes <- counts[order(-`-1`, -`1`, -mean.size)][profile.id != 481][1:8]
two.changes[, panel := paste0(
  ifelse(`-1`==2, "p", "n"), #positive or negative label
  profile.id, ".", chromosome)]
some.err <- neuroblastomaProcessed$errors[two.changes, on=list(
  profile.id, chromosome)]
err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
some.err.tall <- melt(
  some.err,
  measure.vars=names(err.colors))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(profile.id + chromosome ~ .)+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value,
    color=variable, size=variable),
    data=some.err.tall)+
  scale_y_continuous(
    "errors",
    breaks=c(0,1),
    limits=c(-0.2, 1.2))+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)

some.err.tall[, value.i := cumsum(
  c(FALSE, diff(value) != 0)
), by=list(panel, profile.id, chromosome, variable)]
segs.err.tall <- some.err.tall[, list(
  min.log.lambda=min(min.log.lambda),
  max.log.lambda=max(max.log.lambda),
  value=value[1]
), by=list(panel, profile.id, chromosome, variable, value.i)]
segs.min.tall <- segs.err.tall[, {
  .SD[value==min(value)]
}, by=list(panel, profile.id, chromosome, variable)]
segs.min.err <- segs.min.tall[variable=="errors"]
segs.min.err[, mid.log.lambda := (min.log.lambda+max.log.lambda)/2]
set.seed(1)
size <- segs.min.err[is.finite(mid.log.lambda), log(1+rexp(.N, mean(
  max.log.lambda-min.log.lambda)))]
size <- 0.1
segs.min.err[, pred.log.lambda := ifelse(
  min.log.lambda == -Inf, max.log.lambda-size, ifelse(
    max.log.lambda == Inf, min.log.lambda+size, mid.log.lambda))]
segs.min.err[, interval := ifelse(
  is.finite(mid.log.lambda), "finite", "infinite")]
auc.dt <- segs.min.err[, {
  L <- penaltyLearning::ROChange(
    some.err, .SD, c("panel"))
  with(L, data.table(auc, thresholds[threshold=="predicted"]))
}, by=list(interval)]
roc.dt <- segs.min.err[, {
  L <- penaltyLearning::ROChange(
    some.err, .SD, c("panel"))
  L$roc
}, by=list(interval)]
roc.dt[, row.i := 1:.N, by=interval]
ggplot()+
  theme_bw()+
  geom_path(aes(
    FPR, TPR, color=interval, size=interval),
    data=roc.dt)+
  geom_point(aes(
    FPR, TPR, color=interval, size=interval),
    fill=NA,
    shape=21,
    data=auc.dt)+
  scale_size_manual(values=c(
    finite=2,
    infinite=1))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, 'lines'))+
  facet_grid(. ~ interval)+
  geom_path(aes(
    FPR, TPR),
    data=roc.dt)+
  geom_point(aes(
    FPR, TPR),
    fill=NA,
    shape=21,
    data=auc.dt)+
  geom_text(aes(
    FPR, TPR, label=row.i),
    data=roc.dt)
roc.dt[, .(interval, row.i, min.thresh, max.thresh, fp, fn)]
roc.dt[, min.fp.fn := ifelse(fp<fn, fp, fn)]
roc.dt[, width.thresh := max.thresh-min.thresh]
roc.dt[!(width.thresh==Inf & min.fp.fn==0), list(
  aub=sum(min.fp.fn*width.thresh)
), by=list(interval)]
## AUB=0 for good curve / predictions in infinite interval, AUB=2.58
## for bad cruve / predictions in finite interval.

roc.tall <- melt(
  roc.dt,
  measure.vars=names(err.colors))
vert.dt <- roc.tall[, {
  data.table(
    thresh=min.thresh[-.N],
    value=value[-.N],
    next.value=value[-1])
}, by=list(interval, variable)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, 'lines'))+
  facet_grid(. ~ interval)+
  geom_segment(aes(
    min.thresh, value,
    xend=max.thresh, yend=value,
    size=variable, color=variable),
    data=roc.tall)+
  geom_segment(aes(
    thresh, value,
    xend=thresh, yend=next.value,
    color=variable),
    data=vert.dt)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)

data(neuroblastoma, package="neuroblastoma")
some.profiles <- data.table(neuroblastoma$profiles)[two.changes, on=list(
  profile.id, chromosome)]
some.labels <- data.table(neuroblastoma$annotations)[two.changes, on=list(
  profile.id, chromosome)]
some.changes.list <- list()
some.models.list <- list()
for(profile.i in 1:nrow(two.changes)){
  profile.info <- two.changes[profile.i]
  one.profile <- some.profiles[profile.info, on=list(
    profile.id, chromosome, panel)]
  max.segments <- 20L
  fit <- jointseg::Fpsn(one.profile$logratio, max.segments)
  some.models.list[[profile.i]] <- data.table(
    profile.info, n.segments=1:max.segments, loss=fit$J.est)
  for(n.segments in 2:max.segments){
    end.i <- fit$t.est[n.segments, 1:n.segments]
    before <- end.i[-length(end.i)]
    after <- before+1
    change.pos <- one.profile[, (position[before]+position[after])/2]
    some.changes.list[[paste(profile.i, n.segments)]] <- data.table(
      profile.info, n.segments, change.pos)
  }
}
some.models <- do.call(rbind, some.models.list)
some.changes <- do.call(rbind, some.changes.list)
err.list <- penaltyLearning::labelError(
  some.models, some.labels, some.changes,
  change.var="change.pos",
  label.vars=c("min", "max"),
  model.vars="n.segments",
  problem.vars="panel")

thresh.dt <- roc.dt[, {
  data.table(thresh=seq(
    floor(min(max.thresh)), ceiling(max(min.thresh)), by=0.04))
}]
thresh.pred <- data.table(id=1, thresh.dt)[data.table(
  id=1, segs.min.err), on=list(id), allow.cartesian=TRUE]
thresh.pred[, pred.plus.thresh := pred.log.lambda+thresh]
thresh.pred[, pred0 := pred.plus.thresh]
setkey(thresh.pred, panel, pred0, pred.plus.thresh)
setkey(some.err, panel, min.log.lambda, max.log.lambda)
show.pred <- foverlaps(
  some.err,
  thresh.pred[, .(thresh, interval, panel, pred0, pred.plus.thresh)],
  nomatch=0L)
show.pred[, value.fac := ifelse(errors==0, "correct", "error")]
show.label.errors <- show.pred[err.list$label.errors, nomatch=0L, on=list(
  panel, n.segments)]
show.changes <- show.pred[some.changes, nomatch=0L, on=list(
  panel, n.segments), allow.cartesian=TRUE]
thresh.dt[, thresh0 := thresh]
setkey(thresh.dt, thresh, thresh0)
setkey(roc.dt, min.thresh, max.thresh)
roc.points <- foverlaps(roc.dt, thresh.dt, nomatch=0L)
tail.size <- 0.5
roc.points[, tail.thresh := thresh-tail.size]
setkey(roc.points, interval, tail.thresh, thresh)
setkey(roc.dt, interval, min.thresh, max.thresh)
roc.tails <- foverlaps(
  roc.points[, .(thresh, interval, tail.thresh)],
  roc.dt[, .(interval, min.thresh, max.thresh, FPR, TPR)],
  nomatch=0L)
err.sizes.animint <- c(
  fp=5,
  fn=2.5,
  errors=1.5)
roc.segs <- roc.dt[order(interval, min.thresh), {
  dt <- data.table(
    from.FPR=FPR[-.N],
    to.FPR=FPR[-1],
    from.TPR=TPR[-.N],
    to.TPR=TPR[-1])
  print(dt)
  ##browser()
  dt
}, by=list(interval)]
yexp <- 1
roc.u <- roc.points[, .SD[1], by=list(FPR, TPR, interval)]
roc.color <- "violet"
viz <- animint(
  title="ROC curves for neuroblastoma data with several minima",
  duration=list(thresh=250),
  time=list(variable="thresh", ms=300),
  samples=ggplot()+
    ggtitle("Sample label error curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=400)+
    facet_grid(panel ~ .)+
    scale_fill_manual(
      "Prediction",
      values=c(correct="white", error="black"))+
    geom_tallrect(aes(
      xmin=pred.plus.thresh-tail.size, xmax=pred.plus.thresh,
      key=interval),
      fill=roc.color,
      color=NA,
      alpha=0.2,
      showSelected="thresh",
      data=show.pred)+
    geom_vline(aes(
      xintercept=pred.plus.thresh,
      key=interval,
      linetype=interval),
      color=roc.color,
      showSelected="thresh",
      data=show.pred)+
    geom_segment(aes(
      min.log.lambda, value,
      xend=max.log.lambda, yend=value,
      color=variable, size=variable),
      data=segs.err.tall)+
    geom_point(aes(
      pred.plus.thresh, errors,
      key=interval,
      fill=value.fac),
      color=roc.color,
      size=4,
      showSelected=c("interval", "thresh"),
      data=show.pred)+
    scale_y_continuous(
      "errors",
      breaks=c(0,1),
      limits=c(0-yexp, 1+yexp))+
    scale_x_continuous("log(penalty)")+
    scale_color_manual(values=err.colors)+
    scale_size_manual(values=err.sizes.animint),
  thresholds=ggplot()+
    ggtitle("Total label error curves, select threshold")+
    theme_bw()+
    theme_animint(width=600)+
    facet_grid(. ~ interval)+
    geom_segment(aes(
      min.thresh, value,
      xend=max.thresh, yend=value,
      size=variable, color=variable),
      data=roc.tall)+
    geom_segment(aes(
      thresh, value,
      xend=thresh, yend=next.value,
      color=variable),
      data=vert.dt)+
    geom_tallrect(aes(
      xmin=thresh-tail.size, xmax=thresh, key=1),
      data=thresh.dt,
      fill=roc.color,
      color=NA,
      showSelected="thresh",
      alpha=0.2)+
    make_tallrect(thresh.dt, "thresh")+
    scale_y_continuous(
      "errors",
      breaks=seq(0, 20, by=2))+
    scale_x_continuous("Threshold = constant added to predicted values")+
    scale_color_manual(values=err.colors)+
    scale_size_manual(values=err.sizes.animint),
  roc=ggplot()+
    ggtitle("ROC curves, select threshold")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, 'lines'))+
    facet_grid(. ~ interval)+
    geom_text(aes(
      0.2, 0.7, label=sprintf("AUC=%.2f", auc)),
      data=auc.dt)+
    ## geom_path(aes(
    ##   FPR, TPR),
    ##   alpha=0.2,
    ##   data=roc.dt)+
    geom_segment(aes(
      from.FPR, from.TPR,
      xend=to.FPR, yend=to.TPR),
      data=roc.segs,
      alpha=0.2)+
    geom_path(aes(
      FPR, TPR, key=1),
      color=roc.color,
      showSelected="thresh",
      size=3,
      data=roc.tails)+
    coord_equal()+
    scale_x_continuous(
      "False Positive Rate",
      breaks=seq(0, 1, by=0.2))+
    scale_y_continuous(
      "True Positive Rate",
      breaks=seq(0, 1, by=0.2))+
    geom_point(aes(
      FPR, TPR, label=row.i),
      clickSelects="thresh",
      size=4,
      alpha=0.7,
      data=roc.u)+
    geom_point(aes(
      FPR, TPR, key=1, label=row.i),
      showSelected="thresh",
      color=roc.color,
      size=4,
      data=roc.points),
  profiles=ggplot()+
    ggtitle("Noisy data with predicted changes and label errors")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=500, width=1400)+
    facet_grid(panel ~ interval, scales="free")+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6, fill=annotation),
      data=some.labels,
      alpha=0.5,
      color="grey")+
    scale_linetype_manual(
      "Error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6,
      key=paste(min, max),
      linetype=status),
      data=show.label.errors,
      showSelected=c("thresh"),
      fill=NA,
      size=2,
      color="black")+
    scale_fill_manual(values=penaltyLearning::change.colors)+
    geom_point(aes(
      position/1e6, logratio),
      color="grey50",
      data=some.profiles)+
    xlab("Position on chromosome (Mb = mega bases)")+
    geom_vline(aes(
      xintercept=change.pos/1e6, key=change.pos),
      data=show.changes,
      showSelected="thresh",
      color="green"))
viz$roc
unlink("figure-neuroblastomaProcessed-complex", recursive=TRUE)
animint2dir(viz, "figure-neuroblastomaProcessed-complex")
