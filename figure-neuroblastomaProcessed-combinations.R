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
for(size in size.vec){
  print(size)
  segs.min.err[, pred.log.lambda := ifelse(
    min.log.lambda == -Inf, max.log.lambda-size, ifelse(
      max.log.lambda == Inf, min.log.lambda+size, mid.log.lambda))]
  segs.min.err[, interval := ifelse(
    is.finite(mid.log.lambda), "finite", "infinite")]
  pred.dt <- segs.min.err[combos.tall, on=list(panel, interval)]
  auc.dt.list[[paste(size)]] <- pred.dt[, {
    L <- penaltyLearning::ROChange(
      some.err, .SD, c("panel"))
    L$roc[, min.fp.fn := ifelse(fp<fn, fp, fn)]
    L$roc[, width.thresh := max.thresh-min.thresh]
    aub <- L$roc[!(width.thresh==Inf & min.fp.fn==0), {
      sum(min.fp.fn*width.thresh)
    }]
    with(L, data.table(
      auc, aub, size,
      n.finite=sum(interval=="finite"),
      thresholds[threshold=="predicted"]))
  }, by=list(combo.i)]
}
auc.dt <- do.call(rbind, auc.dt.list)

worst <- auc.dt[which.max(auc)]
worst.combo <- combos.tall[worst, .(panel, interval), on=list(combo.i)]
segs.min.err[, pred.log.lambda := ifelse(
  min.log.lambda == -Inf, max.log.lambda-worst$size, ifelse(
    max.log.lambda == Inf, min.log.lambda+worst$size, mid.log.lambda))]
segs.min.err[, interval := ifelse(
  is.finite(mid.log.lambda), "finite", "infinite")]
pred.dt <- segs.min.err[worst.combo, on=list(panel, interval)]
L <- penaltyLearning::ROChange(
  some.err, pred.dt, c("panel"))
L$auc

L$auc.polygon[, row := 1:.N]
ggplot()+
  geom_polygon(aes(
    FPR, TPR),
    fill="red",
    color="black",
    alpha=0.5,
    data=L$auc.polygon)+
  geom_text(aes(
    FPR, TPR, label=row),
    data=L$auc.polygon)

sel.dt <- L$auc.polygon[row>1, .(row, first=1)]
setkey(sel.dt, first, row)
L$auc.polygon[, row0 := row]
setkey(L$auc.polygon, row, row0)
cum.poly <- foverlaps(sel.dt, L$auc.polygon, nomatch=0L)
cum.poly[, added := ifelse(i.row==row, "new", "old")]
lim <- c(-0.2, 1.2)
gg <- ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=cum.poly)+
  geom_text(aes(
    FPR, TPR, label=row, color=added),
    data=cum.poly)+
  scale_color_manual(values=c("new"="red", old="black"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("i.row", nrow=2)+
  ##facet_grid(. ~ i.row)+
  coord_equal(xlim=lim, ylim=lim)+
  scale_x_continuous(breaks=seq(0, 1, by=0.5), labels=c("0", "0.5", "1"))+
  scale_y_continuous(breaks=seq(0, 1, by=0.5))
png("figure-neuroblastomaProcessed-combinations-worst.png", 12, 3, units="in", res=100)
print(gg)
dev.off()

gg <- ggplot()+
  geom_point(aes(
    aub, auc),
    color="black",
    shape=21,
    size=5,
    fill=NA,
    data=auc.dt)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ size)
print(gg)
auc.dt[order(aub), .(auc, aub, size, combo.i)]

aub.count <- auc.dt[, list(
  combos=.N
), by=list(aub=round(aub, 4), size, auc=round(auc, 4))]
gg <- ggplot()+
  geom_hline(aes(
    yintercept=yint),
    data=data.table(yint=1),
    color="grey50")+
  geom_point(aes(
    aub, auc, fill=combos),
    shape=21,
    size=5,
    data=aub.count)+
  scale_fill_gradient(low="white", high="red")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(size ~ .)+
  geom_text(aes(
    aub, auc, label=combos),
    size=3,
    data=aub.count)+
  scale_y_continuous(
    "Area under ROC curve",
    breaks=seq(0, 1.2, by=0.2))+
  scale_x_continuous(
    "Area under both TP and FP curves")
print(gg)
png("figure-neuroblastomaProcessed-combinations-scatter.png", 12, 9, units="in", res=100)
print(gg)
dev.off()

auc.count <- auc.dt[, list(
  combos=.N
), by=list(n.finite, size, auc=round(auc, 4))]
gg <- ggplot()+
  geom_tile(aes(
    n.finite, auc, fill=combos),
    data=auc.count)+
  geom_point(aes(
    n.finite, auc),
    color="black",
    shape=21,
    size=5,
    fill=NA,
    data=worst)+
  scale_fill_gradient(low="white", high="red")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ size)+
  geom_text(aes(
    n.finite, auc, label=combos),
    size=3,
    data=auc.count)+
  scale_x_continuous(
    "Number of predictions in finite min error interval (other predictions in the infinite min error interval)",
    breaks=unique(auc.count$n.finite))
png("figure-neuroblastomaProcessed-combinations.png", 12, 3, units="in", res=100)
print(gg)
dev.off()

