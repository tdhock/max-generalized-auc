source("packages.R")

data(neuroblastomaProcessed, package="penaltyLearning")

e <- function(label, profile.id, chromosome){
  data.table(label, profile.id=factor(profile.id), chromosome=factor(chromosome))
}
select.dt <- rbind(
  e("positive", 4, 2),
  e("negative", 513, 3))
some.err <- neuroblastomaProcessed$errors[select.dt, .(
  fp, fn, possible.fp, possible.fn,
  min.log.lambda=-max.log.lambda,
  max.log.lambda=-min.log.lambda,
  errors, labels,
  label
), on=list(profile.id, chromosome)]
err.sizes <- c(
  fp=3,
  fn=2)
err.colors <- c(
  fp="red",
  fn="deepskyblue")
some.err.tall <- melt(
  some.err,
  measure.vars=names(err.colors))
leg <- "Error type"
some.err.tall[, Label := paste0("\n", label)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(Label ~ ., labeller=label_both)+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value,
    color=variable, size=variable),
    data=some.err.tall)+
  scale_y_continuous(
    "Label errors",
    breaks=c(0,1),
    limits=c(-0.2, 1.2))+
  scale_color_manual(leg,values=err.colors)+
  scale_size_manual(leg,values=err.sizes)+
  scale_x_continuous(
    "Predicted value f(x)",
    breaks=seq(-2, 6, by=2))+
  coord_cartesian(xlim=c(-3, 5))
png("figure-aum-convexity-profiles.png", 3.5, 2, units="in", res=200)
print(gg)
dev.off()

dmin <- 4.5
dmax <- 6.5
some.err[, fp.diff := c(NA, diff(fp))]
some.err[, fn.diff := c(NA, diff(fn))]
pred.dt <- some.err[fp.diff != 0 | fn.diff != 0, list(
  plist=list(min.log.lambda)
), by=label]
plist <- with(pred.dt, structure(plist, names=label))
grid.dt <- data.table(do.call(expand.grid, plist))
grid.dt[, pred.diff := negative-positive]
border.pred <- grid.dt[
  dmin < pred.diff & pred.diff < dmax]
grid.pred <- data.table(
  pred.diff=seq(dmin, dmax, by=0.02))
grid.pred[, positive := 0]
grid.pred[, negative := pred.diff]
both.pred <- rbind(
  border.pred[, .(positive, negative, pred.diff, differentiable=FALSE)],
  grid.pred[, .(positive, negative, pred.diff, differentiable=TRUE)])
pred.tall <- melt(
  both.pred,
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda")[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[order(pred.diff)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "label")
  with(L, data.table(aum, auc, roc=list(roc)))
}, by=list(pred.diff, differentiable)]
metrics.tall <- melt(
  metrics.wide,
  measure.vars=c("aum", "auc")
)

gg <- ggplot()+
  facet_grid(variable ~ ., scales="free", space="free")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  geom_point(aes(
    pred.diff, value, fill=differentiable),
    size=1,
    shape=21,
    data=metrics.tall[order(-differentiable)])+
  xlab("Difference in predicted values, f(negative) - f(positive)")+
  scale_y_continuous("", breaks=seq(0, 2, by=0.5))
png("figure-aum-convexity.png", 5, 2, units="in", res=200)
print(gg)
dev.off()

