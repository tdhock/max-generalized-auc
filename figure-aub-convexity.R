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
two.changes <- counts[order(-`-1`, -`1`, -mean.size)][profile.id != 481][2:3]
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

grid.by <- 0.05
pred.wide <- data.table(expand.grid(
  p4.2=seq(0, 3, by=grid.by),
  n513.3=seq(-5, -2, by=grid.by)))
mvars <- paste(names(pred.wide))
pred.wide[, combo.i := 1:.N]
pred.tall <- melt(
  pred.wide,
  id.vars="combo.i",
  measure.vars=mvars,
  variable.name="panel",
  value.name="pred.log.lambda")[two.changes, nomatch=0L, on=list(panel)]
aub.dt <- pred.tall[order(combo.i)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "panel")
  roc.dt <- data.table(L$roc)
  roc.dt[, min.fp.fn := ifelse(fp<fn, fp, fn)]
  roc.dt[, width.thresh := max.thresh-min.thresh]
  aub <- roc.dt[!(width.thresh==Inf & min.fp.fn==0), sum(min.fp.fn*width.thresh)]
  pred.errors <- L$thresholds[threshold=="predicted", errors]
  min.errors <- L$thresholds[threshold=="min.error", errors]
  data.table(aub, auc=L$auc, pred.errors, min.errors)
}, by=list(combo.i)]

metrics.tall <- melt(
  aub.dt,
  measure.vars=c("aub", "auc", "pred.errors", "min.errors")
)[pred.wide, on=list(combo.i)]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ variable)+
  geom_tile(aes(
    p4.2, n513.3, fill=value),
    data=metrics.tall)+
  scale_fill_gradient2(low="blue", high="red")+
  coord_equal()
png("figure-aub-convexity-heatmap.png", 6, 6, units="in", res=100)
print(gg)
dev.off()

p0 <- metrics.tall[p4.2==2]
vline.dt <- p0[variable=="aub" & value==0, list(n513.3)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ variable)+
  geom_vline(aes(
    xintercept=n513.3),
    color="grey",
    size=2,
    data=vline.dt)+
  geom_point(aes(
    n513.3, value),
    shape=1,
    data=p0)+
  xlab("predicted log(penalty)")
png("figure-aub-convexity.png", 10, 2, units="in", res=100)
print(gg)
dev.off()

