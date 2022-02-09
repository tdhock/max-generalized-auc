library(ggplot2)
library(data.table)
out.loss.list <- list()
for(seed.csv in Sys.glob("figure-sonar-comparisons-data-seed*.csv")){
  out.loss.list[[seed.csv]] <- data.table::fread(seed.csv)
}
out.loss <- do.call(rbind, out.loss.list)

(max.valid.auc <- out.loss[
  set.name=="validation",
  .SD[which.max(auc), .(iteration, auc, set.name, step.size)],
  by=.(seed, loss.name)])
extremes.selected <- sum(max.valid.auc$step.size %in% range(out.loss$step.size))
if(0 < extremes.selected){
  stop("some extreme step sizes selected, should increase grid")
}
show.loss <- out.loss[set.name != "test"][
  max.valid.auc[,.(seed,loss.name,step.size)],
  on=.(seed, loss.name, step.size)
]

## lots of data.
gg <- ggplot(,aes(
  iteration, auc, color=set.name))+
  geom_line(aes(
    group=paste(step.size, set.name)),
    data=out.loss)+
  facet_grid(
    seed ~ loss.name,
    labeller="label_both",
    scales='free',
    space='fixed')

ggplot(,aes(
  iteration, auc, color=set.name))+
  geom_line(
    data=show.loss)+
  geom_point(
    shape=1,
    data=max.valid.auc)+
  facet_grid(
    seed ~ loss.name,
    labeller="label_both",
    scales='free',
    space='fixed')

one.seed <- 6
point.dt <- max.valid.auc[seed==one.seed]
line.dt <- show.loss[seed==one.seed & iteration<200]
gg <- ggplot(,aes(
  iteration, auc, color=loss.name))+
  facet_grid(
    set.name ~ .,
    labeller="label_both",
    scales='free',
    space='fixed')+
  geom_point(
    shape=1,
    data=point.dt)+
  geom_line(
    size=1,
    data=line.dt)
directlabels::direct.label(gg, "top.polygons")

dl.dt <- rbind(
  point.dt,
  line.dt[, .SD[
    which.max(iteration)
  ], by=loss.name][, names(point.dt),with=FALSE])
ggplot(,aes(
  iteration, auc, color=loss.name))+
  ## directlabels::geom_dl(aes(
  ##   label=loss.name),
  ##   method=list("top.polygons",directlabels::dl.trans(y=y+0.1)),
  ##   data=dl.dt)+
  facet_grid(
    . ~ set.name,
    labeller="label_both",
    scales='free',
    space='fixed')+
  geom_point(
    shape=1,
    data=point.dt)+
  geom_line(
    size=1,
    data=line.dt)

gg <- ggplot(,aes(
  iteration, auc, color=loss.name))+
  geom_point(
    shape=21,
    size=2,
    fill="black",
    data=point.dt)+
  geom_line(aes(
    linetype=set.name),
    size=0.5,
    data=line.dt)+
  ## directlabels::geom_dl(aes(
  ##   y=auc+0.01,
  ##   label=loss.name),
  ##   method=list(cex=0.5,"top.polygons"),
  ##   data=point.dt)+
  scale_linetype_manual(values=c(
    subtrain="dashed",
    validation="solid"))
png(
  "figure-sonar-comparisons-iterations.png",
  width=6, height=4, res=200, units="in")
print(gg)
dev.off()

test.loss <- out.loss[set.name=="test"]
test.selected <- test.loss[
  max.valid.auc[, .(seed, loss.name, step.size, iteration)],
  on=.(seed, loss.name, step.size, iteration)]
ggplot()+
  geom_point(aes(
    auc, loss.name),
    data=test.selected)

test.selected.stats <- test.selected[, .(
  median=median(auc),
  mean=mean(auc),
  sd=sd(auc),
  q25=quantile(auc, 0.25),
  q75=quantile(auc, 0.75)
), by=.(loss.name)][order(mean)]
levs <- test.selected.stats$loss.name
test.selected.stats[, Loss := factor(loss.name, levs)]
test.selected[, Loss := factor(loss.name, levs)]
ggplot()+
  geom_segment(aes(
    q25, Loss,
    xend=q75, yend=loss.name),
    data=test.selected.stats)+
  geom_label(aes(
    median, Loss, label="median"),
    data=test.selected.stats)+
  geom_point(aes(
    auc, Loss),
    data=test.selected)

test.wide <- dcast(
  test.selected,
  seed ~ loss.name,
  value.var="auc")
test.compare <- melt(
  test.wide,
  id.vars=c("seed", "aum.count"),
  variable.name="other.loss.name",
  value.name="other.loss.auc")
test.dt <- test.compare[, {
  L <- t.test(
    aum.count, other.loss.auc,
    alternative="greater",
    paired=TRUE)
  L[c("statistic","p.value")]
}, by=other.loss.name][order(p.value)]
ggplot()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  geom_point(aes(
    other.loss.auc, aum.count),
    data=test.compare)+
  coord_equal()+
  facet_grid(. ~ other.loss.name)

test.dt[, Loss := factor(other.loss.name, levs)]
text.x <- Inf
text.size <- 3
text.hjust <- 1
gg <- ggplot()+
  geom_segment(aes(
    mean-sd, Loss,
    xend=mean+sd, yend=Loss),
    data=test.selected.stats)+
  geom_point(aes(
    mean, Loss),
    data=test.selected.stats)+
  geom_text(aes(
    text.x, Loss, label=sprintf("Diff=%.1f p=%.3f", statistic, p.value)),
    vjust=-0.2,
    size=text.size,
    hjust=text.hjust,
    data=test.dt)+
  scale_x_continuous(
    "Test AUC")
png("figure-sonar-comparisons.png", width=4, height=1.5, res=200, units="in")
print(gg)
dev.off()
