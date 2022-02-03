library(ggplot2)
library(data.table)
out.loss <- data.table::fread("figure-sonar-comparisons-data.csv")

(max.valid.auc <- out.loss[
  set.name=="validation",
  .SD[which.max(auc), .(iteration, auc, set.name, step.size)],
  by=.(seed, loss.name)])
show.loss <- out.loss[set.name != "test"][
  max.valid.auc[,.(seed,loss.name,step.size)],
  on=.(seed, loss.name, step.size)
]

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
  q25=quantile(auc, 0.25),
  q75=quantile(auc, 0.75)
), by=.(loss.name)][order(median)]
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

gg <- ggplot()+
  geom_segment(aes(
    q25, Loss,
    xend=q75, yend=loss.name),
    data=test.selected.stats)+
  geom_point(aes(
    median, Loss),
    data=test.selected.stats)+
  scale_x_continuous(
    "Test AUC")
png("figure-sonar-comparisons.png", width=4, height=1.3, res=200, units="in")
print(gg)
dev.off()
