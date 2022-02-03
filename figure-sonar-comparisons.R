library(ggplot2)
library(data.table)
out.loss <- data.table::fread("figure-sonar-comparisons.csv")

(max.valid.auc <- out.loss[
  set.name=="validation",
  .SD[which.max(auc), .(iteration, auc, set.name)],
  by=.(seed, loss.name)])
ggplot(,aes(
  iteration, auc, color=set.name))+
  geom_line(
    data=out.loss)+
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
  max.valid.auc[, .(seed, loss.name, iteration)],
  on=.(seed, loss.name, iteration)]
ggplot()+
  geom_point(aes(
    auc, loss.name),
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
ggplot()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  geom_point(aes(
    other.loss.auc, aum.count),
    data=test.compare)+
  coord_equal()+
  facet_grid(. ~ other.loss.name)

