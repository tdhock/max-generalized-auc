library(data.table)
library(ggplot2)

loss.wide <- data.table::fread("figure-compare-hinge-loss-data.csv")

loss.wide[, neg.auc := -auc]
loss.tall <- melt(
  loss.wide,
  measure.vars = c("auc", "aum", "hinge.loss", "aum.marg"),
  variable.name="loss.name",
  value.name="loss.value")
normalize <- function(x)(x-min(x))/(max(x)-min(x))
loss.tall[, loss.norm := normalize(loss.value), by=loss.name]
rect.dt <- data.table(
  xmin=0, xmax=Inf,
  ymin=-Inf, ymax=0)
loss.contour <- loss.tall[loss.name != "auc"]
show.breaks <- seq(1, 7, by=1)
gg <- ggplot()+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_tile(aes(
    pos, neg, fill=loss.value),
    data=loss.contour)+
  geom_rect(aes(
    xmin=xmin, xmax=xmax,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="black",
    data=rect.dt)+
  geom_contour(aes(
    pos, neg, z=loss.value),
    breaks=show.breaks,
    size=1,
    data=loss.contour)+
  metR::geom_label_contour(aes(
    pos, neg, z=loss.value),
    breaks=show.breaks,
    size=2,
    data=loss.contour)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  scale_fill_gradient(
    "Loss\nvalues",
    low="white",
    high="red")+
  coord_equal()+
  xlab("Real-valued prediction for positive label")+
  ylab("Real-valued prediction\nfor negative label")
png(
  "figure-compare-hinge-loss-contours.png", 
  width=5, height=2, res=200, units="in")
print(gg)
dev.off()

gg <- ggplot()+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
    data=loss.tall)+
  ## geom_contour(aes(
  ##   pos, neg, z=loss.norm),
  ##   data=loss.tall)+
  geom_rect(aes(
    xmin=xmin, xmax=xmax,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="black",
    data=rect.dt)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  scale_fill_gradient(
    "Relative\nvalues",
    low="white",
    high="blue")+
  coord_equal()+
  xlab("Real-valued prediction for positive label")+
  ylab("Real-valued prediction\nfor negative label")
png("figure-compare-hinge-loss.png", width=6, height=2, res=200, units="in")
print(gg)
dev.off()
