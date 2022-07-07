library(data.table)
library(ggplot2)

loss.wide <- data.table::fread("figure-compare-hinge-loss-data.csv")
loss.tall <- melt(
  loss.wide,
  measure.vars = c("AUC", "AUM", "hinge.loss", "AUM.margin", "logistic.loss"),
  variable.name="loss.name",
  value.name="loss.value")
normalize <- function(x)(x-min(x))/(max(x)-min(x))
loss.tall[, loss.norm := normalize(loss.value), by=loss.name]
rect.dt <- data.table(
  xmin=0, xmax=Inf,
  ymin=-Inf, ymax=0)
loss.contour <- loss.tall[loss.name %in% c("AUM", "AUM.margin", "hinge.loss")]
show.breaks <- seq(1, 7, by=1)
gg <- ggplot()+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_tile(aes(
    pos, neg, fill=loss.value),
    data=loss.contour)+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=-Inf,y=-Inf,label="   correct rank"),
    color="grey50",
    vjust=1,
    angle=45,
    hjust=0)+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=Inf,y=-Inf,label="correct \nlabel "),
    hjust=1,
    vjust=-0.2)+
  geom_rect(aes(
    xmin=xmin, xmax=xmax,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="black",
    data=rect.dt)+
  metR::geom_contour2(aes(
    pos, neg, z=loss.value, label=stat(level)),
    breaks=show.breaks,
    color="blue",
    size=1,
    data=loss.contour)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  scale_fill_gradient(
    "Loss\nvalues",
    low="white",
    high="red")+
  coord_equal()+
  geom_abline(aes(
    intercept=intercept, slope=slope),
    color="grey",
    data=data.table(intercept=0, slope=1))+
  xlab("Real-valued prediction for positive label")+
  ylab("Real-valued prediction for negative label")
png(
  "figure-compare-hinge-loss-contours.png", 
  width=8, height=3, res=200, units="in")
print(gg)
dev.off()

log.aum <- loss.tall[loss.name %in% c("AUM","logistic.loss")]
gg <- ggplot()+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_tile(aes(
    pos, neg, fill=loss.value),
    data=log.aum)+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=-Inf,y=-Inf,label="   correct rank"),
    color="grey50",
    vjust=1,
    angle=45,
    hjust=0)+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=Inf,y=-Inf,label="correct \nlabel "),
    hjust=1,
    vjust=-0.2)+
  geom_rect(aes(
    xmin=xmin, xmax=xmax,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="black",
    data=rect.dt)+
  metR::geom_contour2(aes(
    pos, neg, z=loss.value, label=stat(level)),
    breaks=show.breaks,
    color="blue",
    size=1,
    data=log.aum)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  scale_fill_gradient(
    "Loss\nvalues",
    low="white",
    high="red")+
  coord_equal()+
  geom_abline(aes(
    intercept=intercept, slope=slope),
    color="grey",
    data=data.table(intercept=0, slope=1))+
  xlab("Predicted score for positive label")+
  ylab("Predicted score for negative label")
png(
  "figure-compare-hinge-loss-contours-logistic.png", 
  width=5.5, height=3, res=200, units="in")
print(gg)
dev.off()

gg <- ggplot()+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
    data=loss.tall)+
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
