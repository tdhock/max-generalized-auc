source("packages.R")

data.list <- readRDS("figure-unbalanced-grad-desc-data.rds")
x.lab <- "Test AUC, median and quartiles over 10 random train sets"
seeds.wide <- dcast(
  data.list[["result"]],
  prop.pos + seed ~ model,
  value.var="auc")
seeds.tall <- melt(
  seeds.wide,
  id.vars=c("prop.pos","seed","aum.count"),
  variable.name="baseline.name",
  value.name="baseline.auc")
p.tall <- seeds.tall[, {
  result <- t.test(aum.count, baseline.auc, alternative="greater", paired=TRUE)
  with(result, data.table(p.value))
}, by=.(prop.pos, baseline.name)]
dcast(p.tall, baseline.name ~ prop.pos, value.var="p.value")

result.tall <- melt(data.list[["result"]], measure.vars=c("accuracy", "auc"))
result.tall[, percent.positive.labels := factor(prop.pos*100)]
ggplot()+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_point(aes(
    percent.positive.labels, value, color=model),
    data=result.tall)

result.tall[variable=="auc"]

result.stats <- result.tall[, .(
  max=max(value),
  q75=quantile(value, 0.75),
  median=median(value),
  q25=quantile(value, 0.25),
  min=min(value),
  seeds=.N
), by=.(variable, prop.pos, `percent\npositive\nlabels`=percent.positive.labels, model=sub("aum", "AUM", model))]
dcast(result.stats[variable=="auc"], model ~ prop.pos, value.var="median")

glmnet.stats <- result.stats[grepl("glmnet", model)]
gg <- ggplot()+
  ggtitle(paste0(
    "cv.glmnet run on data sets with same number of observations, N=",
    data.list[["N.obs"]],
    "\nand with different proportions of positive labels"))+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_ribbon(aes(
    prop.pos, ymin=min, ymax=max, fill=model),
    alpha=0.5,
    data=glmnet.stats)+
  geom_line(aes(
    prop.pos, median, color=model),
    data=glmnet.stats)+
  scale_x_continuous(
    "Proportion positive labels in train set",
    breaks=unique(result.stats[["prop.pos"]]))+
  ylab("Accuracy or AUC of predictions
on a test set of 50% positive
and 50% negative labels")
gg

logistic.stats <- result.stats[grepl("logistic", model)]
gg <- ggplot()+
  ggtitle(paste0(
    "Logistic regression grad descent run on data sets with same number of observations, N=",
    data.list[["N.obs"]],
    "\nand with different proportions of positive labels"))+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_ribbon(aes(
    prop.pos, ymin=min, ymax=max, fill=model),
    alpha=0.5,
    data=logistic.stats)+
  geom_line(aes(
    prop.pos, median, color=model),
    data=logistic.stats)+
  scale_x_continuous(
    "Proportion positive labels in train set",
    breaks=unique(result.stats[["prop.pos"]]))+
  ylab("Accuracy or AUC of predictions
on a test set of 50% positive
and 50% negative labels")
gg

log.glm.stats <- result.stats[variable=="auc" & grepl("logistic|glmnet", model)]
log.glm.stats[, regularization := ifelse(
  grepl("logistic", model), "early stopping", "L2 norm")]
log.glm.stats[, weights := ifelse(
  grepl("weighted|balanced", model), "balanced", "identity")]
x.lo <- 0.984
gg <- ggplot()+
  ggtitle("Comparing logistic regression models (control experiment)")+
  scale_x_continuous(
    x.lab)+
  coord_cartesian(xlim=c(x.lo, 1))+
  geom_point(aes(
    ifelse(median<x.lo, -Inf, median), weights),
    shape=1,
    data=log.glm.stats)+
  geom_segment(aes(
    q25, weights,
    xend=q75, yend=weights),
    data=log.glm.stats)+
  facet_grid(`percent\npositive\nlabels` ~ regularization, labeller=label_both, scales="free")
png("figure-unbalanced-grad-desc-logistic.png", width=6, height=3, res=200, units="in")
print(gg)
dev.off()

aum.stats <- result.stats[grepl("AUM", model)]
gg <- ggplot()+
  ggtitle(paste0(
    "AUM gradient descent with early stopping run on data sets
with same number of observations, N=",
    data.list[["N.obs"]],
    "\nand with different proportions of positive labels"))+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_ribbon(aes(
    prop.pos, ymin=min, ymax=max, fill=model),
    alpha=0.5,
    data=aum.stats)+
  geom_line(aes(
    prop.pos, median, color=model),
    data=aum.stats)+
  scale_x_continuous(
    "Proportion positive labels in train set",
    breaks=unique(result.stats[["prop.pos"]]))+
  ylab("Accuracy or AUC of predictions
on a test set of 50% positive
and 50% negative labels")
gg

aum.stats.auc <- aum.stats[variable=="auc"]
gg <- ggplot()+
  ggtitle("(a) Comparing AUM variants")+
  scale_x_continuous(
    x.lab,
    limits=c(0.985, 1))+
  ylab("Loss function")+
  geom_point(aes(
    median, model),
    shape=1,
    data=aum.stats.auc)+
  geom_segment(aes(
    q25, model,
    xend=q75, yend=model),
    data=aum.stats.auc)+
  facet_grid(`percent\npositive\nlabels` ~ ., labeller=label_both, scales="free")
png("figure-unbalanced-grad-desc-aum.png", width=5, height=2.5, units="in", res=200)
print(gg)
dev.off()

levs <- c("AUM.count", "squared.hinge.all.pairs", "logistic.weighted")
compare.stats <- result.stats[model %in% levs]
compare.stats[, model.fac := factor(model, levs)]
gg <- ggplot()+
  ggtitle(paste0(
    "AUM gradient descent with early stopping run on data sets with same number of observations, N=",
    data.list[["N.obs"]],
    "\nand with different proportions of positive labels"))+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_ribbon(aes(
    prop.pos, ymin=min, ymax=max, fill=model),
    alpha=0.5,
    data=compare.stats)+
  geom_line(aes(
    prop.pos, median, color=model),
    data=compare.stats)+
  scale_x_continuous(
    "Proportion positive labels in train set",
    breaks=unique(result.stats[["prop.pos"]]))+
  ylab("Accuracy or AUC of predictions
on a test set of 50% positive
and 50% negative labels")
gg

compare.stats.auc <- compare.stats[variable=="auc"]
gg <- ggplot()+
  ggtitle("(b) AUM compared to baselines")+
  scale_x_continuous(
    x.lab,
    limits=c(0.990, 1))+
  ylab("Loss function")+
  geom_point(aes(
    median, model.fac),
    shape=1,
    data=compare.stats.auc)+
  geom_segment(aes(
    q25, model.fac,
    xend=q75, yend=model.fac),
    data=compare.stats.auc)+
  facet_grid(`percent\npositive\nlabels` ~ ., labeller=label_both, scales="free")
png("figure-unbalanced-grad-desc.png", width=5, height=2.5, units="in", res=200)
print(gg)
dev.off()
