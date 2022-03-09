source("packages.R")
result.list <- readRDS("figure-aum-optimized-data.rds")
iterations.tall <- melt(result.list$iterations, id="step.number")
iterations.tall[, Variable := ifelse(
  variable=="min.errors", "label errors", toupper(variable))]
gg <- ggplot()+
  geom_line(aes(
    step.number, value),
    data=iterations.tall)+
  facet_grid(Variable ~ ., scales="free")+
  ylab("")+
  xlab("Iteration of gradient descent algorithm")
png("figure-aum-optimized-iterations.png", width=3, height=3, units="in", res=200)
print(gg)
dev.off()

two.it <- iterations.tall[step.number %in% range(step.number)]
two.it[, emph := ifelse(step.number==1, "initial", "optimized")]
emph <- gg+
  geom_point(aes(
    step.number, value, color=emph),
    data=two.it)+
  theme(legend.position = "none")
png("figure-aum-optimized-iterations-emph.png", width=3, height=3, units="in", res=200)
print(emph)
dev.off()

result.list$auc[, `:=`(x=c(0.25), y=c(0.75, 0.5))]
result.list$roc[, is.monotonic := c(NA, (diff(fp)<=0) & (diff(fn)>=0)), by=pred.name]
result.list$roc[, line.i := 1:.N, by=pred.name]
result.list$roc[, range(line.i), by=pred.name]
## same number of monotonic moves, makes sense because all the same
## error diffs but happening at different thresholds.
result.list$roc[, .(n.monotonic=sum(is.monotonic[-1])), by=pred.name]
ggplot()+
  geom_path(aes(
    FPR, TPR, color=pred.name),
    data=result.list$roc)+
  geom_point(aes(
    FPR, TPR, color=pred.name),
    data=result.list$roc[is.monotonic==FALSE])+
  geom_text(aes(
    FPR, TPR, color=pred.name, label=line.i),
    hjust=1,vjust=1,
    data=result.list$roc[is.monotonic==FALSE])+
  geom_point(aes(
    FPR, TPR, color=pred.name),
    fill="white",
    shape=21,
    data=result.list$auc)+
  geom_segment(aes(
    x, y,
    xend=FPR, yend=TPR,
    color=pred.name),
    data=result.list$auc)+
  geom_label(aes(
    x, y, color=pred.name,
    label=sprintf(
      "%s errors=%d auc=%.2f",
      pred.name, errors, auc)),
    size=3,
    hjust=0,
    data=result.list$auc)+
  coord_equal()+
  guides(color="none")

gg <- ggplot()+
  geom_path(aes(
    FPR, TPR, color=pred.name),
    data=result.list$roc)+
  geom_point(aes(
    FPR, TPR, color=pred.name),
    fill="white",
    shape=21,
    data=result.list$auc)+
  geom_segment(aes(
    x, y,
    xend=FPR, yend=TPR,
    color=pred.name),
    data=result.list$auc)+
  geom_label(aes(
    x, y, color=pred.name,
    label=sprintf(
      "%s errors=%d auc=%.2f",
      pred.name, errors, auc)),
    size=3,
    hjust=0,
    data=result.list$auc)+
  coord_equal()+
  guides(color="none")
png("figure-aum-optimized.png", width=3, height=3, units="in", res=200)
print(gg)
dev.off()
