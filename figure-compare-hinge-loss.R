library(data.table)
library(ggplot2)

loss.wide <- data.table::fread("figure-compare-hinge-loss-data.csv")
loss.tall <- melt(
  loss.wide,
  measure.vars = c(
    "AUC", "AUM", "AUM.margin", "logistic.loss",
    grep("pairwise.*_", names(loss.wide), value=TRUE),
    "squared.hinge", "hinge.loss"),
  variable.name="loss.name",
  value.name="loss.value")
normalize <- function(x)(x-min(x))/(max(x)-min(x))
loss.tall[, loss.norm := normalize(loss.value), by=loss.name]
rect.dt <- data.table(
  xmin=0, xmax=Inf,
  ymin=-Inf, ymax=0)
levs <- c("AUM", "AUM.margin", "hinge.loss")
loss.contour <- loss.tall[loss.name %in% levs]
loss.contour[, loss.fac := factor(loss.name, levs)]
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
  facet_grid(. ~ loss.fac)+
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

loss.line <- loss.tall[pos+neg==0 & grepl("pairwise",loss.name)]
loss.line[, loss := sub("pred.prob","last.act=sigmoid",sub("pred.real","last.act=none",sub("pairwise[.]","",sub("_","\n",loss.name))))]
gg <- ggplot()+
  ggtitle("Comparing pairwise loss functions")+
  geom_line(aes(
    pos, loss.value, color=loss),
    size=1,
    data=loss.line)+
  scale_x_continuous(
    "Predicted score yhat for positive example",
    breaks=seq(-3, 3, by=1))+
  coord_cartesian(xlim=c(-5, 3),ylim=c(0,9))+
  scale_y_continuous(
    "Pairwise loss after last activation (a), Loss[a(yhat)-a(-yhat)]",
    limits=c(0,10),
    breaks=seq(0,8,by=2))
dl <- directlabels::direct.label(gg, "left.polygons")
png(
  "figure-compare-hinge-loss-pairwise-line.png", 
  width=5, height=5, res=200, units="in")
print(dl)
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
    color="violet",
    size=1,
    data=log.aum)+
  geom_label(aes(
    x,y,label=label),
    data=data.table(x=-Inf,y=-Inf,label="   correct rank"),
    color="grey50",
    alpha=0.5,
    vjust=1,
    angle=45,
    hjust=0)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ loss.name)+
  scale_fill_gradient(
    "Loss\nvalues",
    low="white",
    high=rgb(0.3,0.5,0.3))+
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

log.aum <- loss.tall[loss.name %in% c(
  "pairwise.squared.hinge_pred.real",
  "pairwise.squared.hinge_pred.prob")]
log.aum[, last.layer.output := sub(".*[.]", "", loss.name)]
show.breaks <- seq(0, 3.5, by=0.5)
gg <- ggplot()+
  ggtitle("Pairwise squared hinge loss functions")+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
  facet_grid(. ~ last.layer.output, labeller=label_both)+
  scale_fill_gradient(
    "Relative\nloss\nvalues",
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
  "figure-compare-hinge-loss-squared-pairwise-relative.png", 
  width=5.5, height=3, res=200, units="in")
print(gg)
dev.off()

show.breaks <- seq(0, 3.5, by=0.5)
gg <- ggplot()+
  ggtitle("Pairwise squared hinge loss functions")+
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
  facet_grid(. ~ last.layer.output, labeller=label_both)+
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
  "figure-compare-hinge-loss-squared-pairwise.png", 
  width=5.5, height=3, res=200, units="in")
print(gg)
dev.off()

log.aum <- loss.tall[loss.name %in% c(
  "pairwise.hinge.loss_pred.real",
  "pairwise.hinge.loss_pred.prob")]
log.aum[, last.layer.output := sub(".*[.]", "", loss.name)]
show.breaks <- seq(0, 1.75, by=0.25)
gg <- ggplot()+
  ggtitle("Pairwise linear hinge loss functions")+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
  facet_grid(. ~ last.layer.output, labeller=label_both)+
  scale_fill_gradient(
    "Relative\nloss\nvalues",
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
  "figure-compare-hinge-loss-hinge-pairwise-relative.png", 
  width=5.5, height=3, res=200, units="in")
print(gg)
dev.off()

show.breaks <- seq(0, 1.75, by=0.25)
gg <- ggplot()+
  ggtitle("Pairwise linear hinge loss functions")+
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
    data=data.table(x=Inf,y=-Inf,label="correct \nlabels "),
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
  facet_grid(. ~ last.layer.output, labeller=label_both)+
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
  "figure-compare-hinge-loss-hinge-pairwise.png", 
  width=5.5, height=3, res=200, units="in")
print(gg)
dev.off()

gg <- ggplot()+
  geom_tile(aes(
    pos, neg, fill=loss.norm),
    data=loss.tall[loss.name %in% c("AUC", "AUM", "logistic.loss")])+
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
    high="red")+
  coord_equal()+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=-3,y=-3,label="   correct rank"),
    color="grey50",
    vjust=1.3,
    angle=45,
    hjust=0)+
  geom_text(aes(
    x,y,label=label),
    data=data.table(x=0,y=0,label="correct \nlabels "),
    hjust=-0.1,
    vjust=1.1)+
  xlab("Real-valued prediction for positive label")+
  ylab("Real-valued prediction\nfor negative label")
png("figure-compare-hinge-loss.png", width=7, height=2.6, res=200, units="in")
print(gg)
dev.off()
