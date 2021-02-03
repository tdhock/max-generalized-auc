source("packages.R")

d <- function(min.log.lambda, fp, fn){
  data.table(min.log.lambda, fp, fn)
}
profile <- function(..., possible.fp, possible.fn, errors, labels){
  dt <- do.call(rbind, list(...))
  if(missing(possible.fp))possible.fp <- max(dt$fp)
  if(missing(possible.fn))possible.fn <- max(dt$fn)
  errors <- dt[, fp+fn]
  if(missing(labels))labels <- max(errors)
  dt[, data.table(
    min.log.lambda,
    max.log.lambda=c(min.log.lambda[-1], Inf),
    fp, fn, errors, possible.fp, possible.fn, labels)]
}
profile.list <- list(
  negative=profile(
    d(-Inf, 0, 0),
    d(0, 1, 0)),
  positive=profile(
    d(-Inf, 0, 1),
    d(0, 0, 0)))
profile.wide <- data.table(
  label=names(profile.list)
)[, profile.list[[label]], by=label]
err.sizes <- c(
  fp=3,
  fn=2)
err.colors <- c(
  fp="red",
  fn="deepskyblue")
profile.tall <- data.table::melt(
  profile.wide,
  measure=c("fp", "fn"))
leg <- "Error type"
gg <- ggplot()+
  facet_grid(. ~ label, labeller=label_both)+
  ## theme_bw()+
  ## theme(panel.spacing=grid::unit(0, "lines"))+
  geom_segment(aes(
    min.log.lambda, value,
    color=variable, size=variable,
    xend=max.log.lambda, yend=value),
    data=profile.tall)+
  scale_y_continuous(
    "Label errors",
    breaks=c(0,1),
    limits=c(-0.2, 1.2))+
  scale_color_manual(leg, values=err.colors)+
  scale_size_manual(leg, values=err.sizes)+
  scale_x_continuous(
    "Predicted value f(x)",
    limits=c(-1.8, 1.8))
png("figure-binary-class.png", width=4, height=2, res=200, units="in")
print(gg)
dev.off()
