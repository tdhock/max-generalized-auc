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
(expand.args <- sapply(two.changes$panel, function(L){
  c("finite", "infinite")
}, USE.NAMES=TRUE, simplify=FALSE))
combos.wide <- data.table(do.call(expand.grid, expand.args))
combos.wide[, combo.i := 1:.N]
combos.tall <- melt(
  combos.wide,
  id.vars="combo.i",
  variable.name="panel",
  value.name="interval")

auc.dt.list <- list()
size.vec <- 10^seq(-4, 2)
##size.vec <- 1e-4
pred.dt.list <- list()
for(size in size.vec){
  print(size)
  segs.min.err[, pred.log.lambda := ifelse(
    min.log.lambda == -Inf, max.log.lambda-size, ifelse(
      max.log.lambda == Inf, min.log.lambda+size, mid.log.lambda))]
  segs.min.err[, interval := ifelse(
    is.finite(mid.log.lambda), "finite", "infinite")]
  pred.dt <- segs.min.err[combos.tall, on=list(panel, interval)]
  pred.dt.list[[paste(size)]] <- pred.dt
  auc.dt.list[[paste(size)]] <- pred.dt[, {
    L <- penaltyLearning::ROChange(
      some.err, .SD, c("panel"))
    L$roc[, min.fp.fn := ifelse(fp<fn, fp, fn)]
    L$roc[, width.thresh := max.thresh-min.thresh]
    aub <- L$roc[!(width.thresh==Inf & min.fp.fn==0), {
      sum(min.fp.fn*width.thresh)
    }]
    with(L, data.table(
      auc, aub, size, roc=list(list(roc)),
      n.finite=sum(interval=="finite"),
      thresholds[threshold=="predicted"]))
  }, by=list(combo.i)]
}

neuroblastomaProcessed.combinations <- list(
  problems=two.changes,
  segs.min.err=segs.min.err,
  some.err=some.err,
  combos=combos.tall,
  pred=do.call(rbind, pred.dt.list),
  auc=do.call(rbind, auc.dt.list))
saveRDS(
  neuroblastomaProcessed.combinations,
  "neuroblastomaProcessed.combinations.rds")

