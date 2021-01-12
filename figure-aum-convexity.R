source("packages.R")

data(neuroblastomaProcessed, package="penaltyLearning")

e <- function(label, profile.id, chromosome){
  data.table(label, profile.id=factor(profile.id), chromosome=factor(chromosome))
}
select.dt <- rbind(
  e("positive", 4, 2),
  e("negative", 513, 3))
some.err <- neuroblastomaProcessed$errors[select.dt, on=list(
  profile.id, chromosome)]
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
  scale_x_continuous("Predicted value f(x)")
png("figure-aum-convexity-profiles.png", 3.5, 2, units="in", res=200)
print(gg)
dev.off()

pred.wide <- data.table(
  pred.diff=seq(4.5, 6.5, by=0.01))
pred.wide[, positive := pred.diff]
pred.wide[, negative := 0]
pred.tall <- melt(
  pred.wide,
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda")[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[order(pred.diff)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "label")
  with(L, data.table(aum, auc, roc=list(roc)))
}, by=list(pred.diff)]
metrics.tall <- melt(
  metrics.wide,
  measure.vars=c("aum", "auc")
)

gg <- ggplot()+
  facet_grid(variable ~ ., scales="free", space="free")+
  geom_point(aes(
    pred.diff, value),
    data=metrics.tall)+
  xlab("Difference in predicted values, f(positive) - f(negative)")+
  scale_y_continuous("", breaks=seq(0, 2, by=0.5))
png("figure-aum-convexity.png", 4, 2, units="in", res=200)
print(gg)
dev.off()

