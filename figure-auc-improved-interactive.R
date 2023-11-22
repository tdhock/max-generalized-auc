library(data.table)
library(animint2)
auc.optimized <- readRDS("auc.improved.rds")[
  initialization=="min.error"
][
, predictions := sub("improved", "optimized", pred.name)
]
auc.optimized[, set.fold := paste0(set.name, "/", fold)]
roc.dt.list <- list()
for(test.fold.i in 1:nrow(auc.optimized)){
  one.fold <- auc.optimized[test.fold.i]
  roc.dt.list[[test.fold.i]] <- one.fold[, data.table(
    set.fold, predictions, roc[[1]])]
}
(roc.dt <- do.call(rbind, roc.dt.list))
roc.dt[, fn0 := fn-min(fn), by=.(set.fold)]
roc.dt[, min.fp.fn := ifelse(fp<fn0, fp, fn0)]
roc.dt[, min.diff := c(0,diff(min.fp.fn)), by=.(set.fold,predictions)]
roc.dt[, cum.abs.diff := cumsum(abs(min.diff)), by=.(set.fold,predictions)]
text.info <- roc.dt[
  cum.abs.diff==0
][
, data.table(text.thresh=max.thresh)[which.max(max.thresh)]
, by=.(set.fold,predictions)
]
roc.dt[, width := max.thresh-min.thresh]
roc.dt[, area := ifelse(min.fp.fn==0, 0, min.fp.fn*width)]
(aum.dt <- roc.dt[, .(
  aum=sum(area)
), by=.(set.fold, predictions)
][
  text.info, on=.(set.fold,predictions)
][order(aum)])
aum.dt[, log.aum := log10(aum+1)]
aum.wide <- dcast(aum.dt, set.fold ~ predictions, value.var="log.aum")
aum.wide[, status := ifelse(
  initial==optimized, "same", ifelse(
    initial>optimized, "better", "worse"))]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    auc, set.name, color=predictions),
    data=auc.optimized)

auc.optimized[, accuracy.percent := 100-error.percent]
auc.tall <- melt(auc.optimized, measure.vars=c("accuracy.percent", "auc"))
auc.stats <- auc.tall[, .(
  mean=mean(value),
  sd=sd(value)
), by=.(set.name, variable, predictions)]
auc.only <- auc.stats[variable=="auc" & predictions=="optimized"][order(mean)]
set.levs <- auc.only$set.name
auc.stats[, set.fac := factor(set.name, set.levs)]
auc.tall[, set.fac := factor(set.name, set.levs)]
auc.only.wide <- dcast(
  auc.stats[variable=="auc"], set.name ~ predictions, value.var="mean")
auc.only.wide[, diff := optimized - initial]
auc.only.wide[order(diff)]

roc.wide <- dcast(
  auc.optimized,
  set.name + fold + set.fold ~ predictions,
  value.var=c("auc", "accuracy.percent"))
roc.wide[, auc_status := ifelse(
  auc_initial==auc_optimized, "same", ifelse(
    auc_initial<auc_optimized, "better", "worse"))]
roc.wide[auc_initial>auc_optimized]
roc.wide[, accuracy.percent_status := ifelse(
  accuracy.percent_initial==accuracy.percent_optimized, "same", ifelse(
    accuracy.percent_initial<accuracy.percent_optimized, "better", "worse"))]

err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
roc.dt[, seg.i := 1:.N, by=.(set.fold, predictions)]
roc.dt[, mid.thresh := (min.thresh+max.thresh)/2]
roc.tall <- melt(
  roc.dt,
  measure.vars=names(err.sizes))
status.colors <- c(
  same="black",
  better="blue",
  worse="red")
if(FALSE){
  tallrect.dt <- data.table(
    mid.thresh=seq(-10, 5, by=0.2))
  roc.dots <- roc.dt[tallrect.dt, .(
    set.fold, predictions, mid.thresh=i.mid.thresh, FPR, TPR
  ), on=.(
    min.thresh<mid.thresh, max.thresh>mid.thresh)]
}
cum.dt <- melt(
  roc.dt,
  measure.vars=c("FPR","TPR"),
  id.vars=c("mid.thresh","set.fold","predictions")
)[
, diff := c(0, diff(value)), by=.(set.fold,variable,predictions)
][
, .(L2=sqrt(sum(diff^2))), keyby=.(set.fold,predictions,mid.thresh)
][
, cum.L2 := cumsum(L2), by=.(set.fold,predictions)
][]
n.grid <- 101
grid.dt <- cum.dt[, .(
  cum.grid=seq(0, max(cum.L2), l=n.grid),
  percent.curve=seq(0, 100, l=n.grid)
), by=.(set.fold,predictions)]
grid.thresh <- cum.dt[
  grid.dt,
  data.table(
    set.fold, predictions, mid.thresh, percent.curve
  )
, on=.(set.fold, predictions, cum.L2=cum.grid), roll="nearest"]
grid.roc <- roc.dt[
  grid.thresh,
  .(set.fold, predictions, percent.curve, 
    min.thresh, mid.thresh, max.thresh, FPR, TPR),
  on=.(set.fold, predictions, mid.thresh)
]
grid.thresh[, .(n.thresh=.N), by=.(set.fold, predictions)]
ggplot()+
  coord_equal()+
  geom_point(aes(
    FPR, TPR, color=predictions),
    data=grid.roc)+
  facet_wrap("set.fold")
ggplot()+
  geom_rect(aes(
    ymin=ymin,
    ymax=ymin+0.5,
    fill=predictions,
    xmin=min.thresh,
    xmax=max.thresh),
    data=grid.roc[, ymin := ifelse(predictions=="initial", 0, 0.5)])+
  facet_wrap("set.fold")
ggplot()+
  geom_line(aes(
    mid.thresh, percent.curve, color=predictions),
    data=grid.roc)+
  facet_wrap("set.fold")
viz <- animint(
  title="Initial/optimized AUM/AUC for change-point problems",
  out.dir="figure-auc-improved-interactive",
  ## ggplot()+
  ##   ggtitle("Data sets ordered by mean optimized AUC")+
  ##   guides(size="none", color="none")+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "lines"))+
  ##   theme_animint(width=800)+
  ##   facet_grid(. ~ variable, scales="free")+
  ##   geom_segment(aes(
  ##     mean+sd, set.fac,
  ##     color=predictions, size=predictions,
  ##     xend=mean-sd, yend=set.fac),
  ##     data=auc.stats)+
  ##   scale_size_manual(values=c(optimized=2, initial=3))+
  ##   geom_point(aes(
  ##     mean, set.fac,
  ##     color=predictions),
  ##     shape=21,
  ##     size=3,
  ##     fill="white",
  ##     data=auc.stats)+
  ##   geom_point(aes(
  ##     value, set.fac,
  ##     color=predictions),
  ##     clickSelects="set.fold",
  ##     alpha=0.6,
  ##     size=4,
  ##     data=auc.tall)+
  ##   xlab("")+
  ##   ylab("Data set"),
  labelErr=ggplot()+
    ggtitle("FP/FN curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(update_axes="y")+
    facet_grid(predictions ~ .)+
    xlab("Prediction threshold")+
    ylab("Incorrectly predicted labels")+
    geom_vline(aes(
      xintercept=(min.thresh+max.thresh)/2,
      key=1),
      data=auc.optimized,
      color="grey",
      showSelected="set.fold")+
    geom_line(aes(
      mid.thresh, value,
      key=variable, group=variable,
      color=variable, size=variable),
      data=roc.tall,
      showSelected="set.fold")+
    geom_text(aes(
      text.thresh-1, text.err,
      key=1,
      label=sprintf("AUM=%.1f", aum)),
      hjust=1,
      showSelected="set.fold",
      data=aum.dt[, text.err := 0],
      color="grey50")+
    geom_polygon(aes(
      mid.thresh, min.fp.fn, key=1),
      fill="grey",
      data=roc.dt,
      size=0,
      showSelected="set.fold")+
    ## make_tallrect(tallrect.dt, "mid.thresh")+
    geom_tallrect(aes(
      xmin=min.thresh,
      xmax=max.thresh,
      key=percent.curve),
      data=grid.roc,
      showSelected="set.fold",
      clickSelects="percent.curve",
      alpha=0.5)+
    scale_size_manual(values=err.sizes)+
    scale_color_manual(values=err.colors),
  selector.types=list(
    variable="multiple"),
  roc=ggplot()+
    ggtitle("ROC curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_path(aes(
      FPR, TPR,
      key=predictions,
      color=predictions, group=predictions),
      data=roc.dt,
      showSelected="set.fold")+
    geom_point(aes(
      FPR, TPR, color=predictions, key=predictions),
      fill="white",
      data=auc.optimized,
      showSelected="set.fold")+
    geom_point(aes(#for smooth transitions.
      FPR, TPR, color=predictions, key=predictions),
      data=grid.roc,
      showSelected=c("set.fold", "percent.curve"),
      size=4,
      alpha=0.5)+
    geom_point(aes(
      FPR, TPR, color=predictions, key=paste(predictions, percent.curve)),
      data=grid.roc,
      showSelected="set.fold",
      clickSelects="percent.curve",
      size=4,
      alpha=0.5),
    ## geom_point(aes(#for smooth transitions.
    ##   FPR, TPR, color=predictions, key=predictions),
    ##   data=roc.dots,
    ##   showSelected=c("set.fold", "mid.thresh"),
    ##   size=4,
    ##   alpha=0.5)+
    ## geom_point(aes(
    ##   FPR, TPR, color=predictions, key=paste(predictions, mid.thresh)),
    ##   data=roc.dots,
    ##   showSelected="set.fold",
    ##   clickSelects="mid.thresh",
    ##   size=4,
    ##   alpha=0.5),
  labelAcc=ggplot()+
    ggtitle("Percent correctly predicted labels")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    scale_color_manual(values=status.colors)+
    guides(color="none")+
    coord_equal()+
    geom_point(aes(
      accuracy.percent_initial, accuracy.percent_optimized,
      key=set.fold),
      clickSelects="set.fold",
      alpha_off=1,
      fill_off="transparent",
      fill="black",
      color="black",
      size=4,
      data=roc.wide),
  logAUM=ggplot()+
    ggtitle("Log[Area under Min(FP,FN) + 1]")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    guides(color="none")+
    scale_color_manual("Status",values=status.colors)+
    geom_point(aes(
      initial, optimized,
      key=set.fold),
      clickSelects="set.fold",
      size=4,
      alpha_off=1,
      fill_off="transparent",
      fill="black",
      color="black",
      data=aum.wide)+
    coord_equal(),
  auc=ggplot()+
    ggtitle("Area under the ROC curve")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    geom_point(aes(
      auc_initial, auc_optimized,
      key=set.fold),
      clickSelects="set.fold",
      size=4,
      alpha_off=1,
      fill_off="transparent",
      fill="black",
      color="black",
      data=roc.wide)+
    coord_equal(),
  duration=list(
    set.fold=500,
    percent.curve=500),
  time=list(
    ms=500,
    variable="percent.curve"),
  first=list(
    set.fold="H3K4me3_XJ_immune/4"
  ),
  source="https://github.com/tdhock/max-generalized-auc/blob/master/figure-auc-improved-interactive.R"
)
animint2dir(viz, viz$out.dir, open.browser=FALSE)
if(FALSE){
  animint2pages(viz, "2023-11-21-auc-improved")
}
