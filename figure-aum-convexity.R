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
  "min(FP,FN)"=1,
  FP=3,
  FN=2)
err.colors <- c(
  "min(FP,FN)"="black",
  FP="red",
  FN="deepskyblue")
some.err.tall <- melt(
  some.err,
  measure.vars=c("fp","fn"),
  variable.name="var.lower")
some.err.tall[, variable := toupper(var.lower)]
leg <- "Error type"
some.err.tall[, Label := paste0("\n", label)]
gg.err <- ggplot()+
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
print(gg.err)
dev.off()

dmin <- 3.5
dmax <- 7.5
some.err[, fp.diff := c(NA, diff(fp)), by=label]
some.err[, fn.diff := c(NA, diff(fn)), by=label]
some.diff <- some.err[fp.diff != 0 | fn.diff != 0, .(
  id=1, label, fp.diff, fn.diff, pred.log.lambda=min.log.lambda)]
some.diff[, fp.cum := cumsum(fp.diff), by=label]
some.diff[, fn.cum := rev(cumsum(rev(-fn.diff))), by=label]
dlist <- split(some.diff, some.diff[["label"]])
grid.dt <- with(dlist, positive[negative, on="id", allow.cartesian=TRUE])
grid.dt[, negative := i.pred.log.lambda]
grid.dt[, positive := pred.log.lambda]
grid.dt[, pred.diff := negative - positive]
grid.sorted <- grid.dt[order(pred.diff), .(
  pred.diff, fn=fn.cum, fp=i.fp.cum)]
grid.sorted[, min.fp.fn := pmin(fp,fn)]
border.pred <- grid.dt[
  dmin < pred.diff & pred.diff < dmax]
grid.pred <- data.table(
  pred.diff=seq(dmin, dmax, by=0.02))
grid.pred[, positive := 0]
grid.pred[, negative := pred.diff]
both.pred <- rbind(
  border.pred[, .(positive, negative, pred.diff, differentiable=FALSE)],
  grid.pred[, .(positive, negative, pred.diff, differentiable=TRUE)])
##positive=0, negative=pred.diff.
pred.tall <- melt(
  both.pred,
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda")[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[order(pred.diff)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "label")
  with(L, data.table(
    aum, auc,
    SM=L$roc[min.thresh < max.thresh, sum(min.fp.fn)],
    roc=list(roc)))
}, by=list(pred.diff, differentiable)]
metrics.wide[auc==max(auc)] #max auc => aum>0.

pred.diff.vec <- c(4.5, 5, 5.14)
vline.dt <- rbind(
  data.table(pred=pred.diff.vec, Label="\nnegative"),
  data.table(pred=0, Label="\npositive"))
vline.dt[, pred.value := pred-1]
gg.err+
  geom_vline(aes(
    xintercept=pred.value),
    data=vline.dt)
## TODO three or more slides showing alignment.

show.roc.dt.list <- list()
for(pdiff in pred.diff.vec){
  select.dt <- data.table(pred.diff=pdiff)
  pdiff.metrics <- metrics.wide[select.dt, on="pred.diff", roll="nearest"]
  pdiff.roc <- pdiff.metrics[["roc"]][[1]]
  show.roc.dt.list[[paste(pdiff)]] <- data.table(
    pred.diff=pdiff,
    pdiff.metrics[, .(AUC=auc, AUM=round(aum,3), SM)],
    pdiff.roc)
}
(show.roc.dt <- do.call(rbind, show.roc.dt.list))
show.roc.dt[, min.fp.fn := pmin(fp, fn)]
show.roc.tall <- melt(
  show.roc.dt,
  measure=c("fp","fn","min.fp.fn"),
  variable.name="lower.var")
show.roc.tall[, variable := ifelse(
  lower.var=="min.fp.fn", "min(FP,FN)", toupper(lower.var))]
gg <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  facet_grid(pred.diff + AUC + AUM ~ ., labeller=label_both)+
  geom_rect(aes(
    xmin=min.thresh, xmax=max.thresh,
    ymin=0, ymax=min.fp.fn),
    fill="grey",
    color=NA,
    alpha=0.5,
    show.roc.dt)+
  geom_segment(aes(
    min.thresh, value,
    xend=max.thresh, yend=value,
    color=variable, size=variable),
    data=show.roc.tall)+
  scale_y_continuous(
    "Label errors",
    breaks=c(0,1))+
  scale_color_manual(leg,values=err.colors)+
  scale_size_manual(leg,values=err.sizes)+
  geom_blank(aes(
    x, y),
    data=data.table(x=0, y=c(-0.4,1.4)))+
  scale_x_continuous(
    "Constant added to predicted values")
png("figure-aum-convexity-thresholds.png", 5, 3.5, units="in", res=200)
print(gg)
dev.off()

metrics.tall <- melt(
  metrics.wide,
  measure.vars=c("aum", "auc", "SM"),
  variable.name="var.lower"
)
metrics.tall[, variable := toupper(var.lower)]
gg <- ggplot()+
  theme(panel.spacing=grid::unit(1, "lines"))+
  theme(text=element_text(size = 15))+
  ##theme(legend.position=c(0.8, 0.15))+
  theme(legend.position="bottom")+
  facet_grid(variable ~ ., scales="free")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  geom_point(aes(
    pred.diff, value, fill=differentiable),
    size=1,
    shape=21,
    data=metrics.tall[order(-differentiable)])+
  xlab("Prediction difference, f(negative) - f(positive)")+
  coord_cartesian(xlim=c(4,7))+
  scale_y_continuous("", breaks=seq(0, 3, by=1))
gg

gg.emph <- gg+
  theme_bw()+
  geom_vline(aes(
    xintercept=pred.diff),
    color="grey50",
    data=data.table(pred.diff=pred.diff.vec))
png("figure-aum-convexity-emph.png", 5, 3, units="in", res=200)
print(gg.emph)
dev.off()

metrics.no.SM <- metrics.tall[variable != "SM"]
gg.no.SM <- ggplot()+
  theme(panel.spacing=grid::unit(1, "lines"))+
  theme(text=element_text(size = 15))+
  theme(legend.position="bottom")+
  facet_grid(variable ~ ., scales="free")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  geom_point(aes(
    pred.diff, value, fill=differentiable),
    size=1,
    shape=21,
    data=metrics.no.SM[order(-differentiable)])+
  xlab("Prediction difference, f(negative) - f(positive)")+
  coord_cartesian(xlim=c(4,7))+
  scale_y_continuous("", breaks=seq(0, 3, by=1))
png("figure-aum-convexity-no-SM.png", 4.2, 3, units="in", res=200)
print(gg.no.SM)
dev.off()

png("figure-aum-convexity.png", 4.2, 3, units="in", res=200)
print(gg)
dev.off()

