library(data.table)
library(animint2)
auc.improved <- readRDS("auc.improved.rds")

auc.improved[, set.fold := paste0(set.name, "/", fold)]
roc.dt.list <- list()
for(test.fold.i in 1:nrow(auc.improved)){
  one.fold <- auc.improved[test.fold.i]
  roc.dt.list[[test.fold.i]] <- one.fold[, data.table(
    set.fold, pred.name, roc[[1]])]
}
(roc.dt <- do.call(rbind, roc.dt.list))
roc.dt[, fn0 := fn-min(fn), by=.(set.fold)]
roc.dt[, min.fp.fn := ifelse(fp<fn0, fp, fn0)]
roc.dt[, width := max.thresh-min.thresh]
roc.dt[, area := ifelse(min.fp.fn==0, 0, min.fp.fn*width)]
(aum.dt <- roc.dt[, .(
  aum=sum(area)
), by=.(set.fold, pred.name)][order(aum)])
aum.dt[, log.aum := log10(aum+1)]
aum.wide <- dcast(aum.dt, set.fold ~ pred.name, value.var="log.aum")
aum.wide[, status := ifelse(
  initial==improved, "same", ifelse(
    initial>improved, "better", "worse"))]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    auc, set.name, color=pred.name),
    data=auc.improved)

auc.improved[, accuracy.percent := 100-error.percent]
auc.tall <- melt(auc.improved, measure.vars=c("accuracy.percent", "auc"))
auc.stats <- auc.tall[, .(
  mean=mean(value),
  sd=sd(value)
), by=.(set.name, variable, pred.name)]
auc.only <- auc.stats[variable=="auc" & pred.name=="improved"][order(mean)]
set.levs <- auc.only$set.name
auc.stats[, set.fac := factor(set.name, set.levs)]
auc.tall[, set.fac := factor(set.name, set.levs)]
auc.only.wide <- dcast(
  auc.stats[variable=="auc"], set.name ~ pred.name, value.var="mean")
auc.only.wide[, diff := improved - initial]
auc.only.wide[order(diff)]

roc.wide <- dcast(
  auc.improved,
  set.name + fold + set.fold ~ pred.name,
  value.var=c("auc", "accuracy.percent"))
roc.wide[, auc_status := ifelse(
  auc_initial==auc_improved, "same", ifelse(
    auc_initial<auc_improved, "better", "worse"))]
roc.wide[auc_initial>auc_improved]
roc.wide[, accuracy.percent_status := ifelse(
  accuracy.percent_initial==accuracy.percent_improved, "same", ifelse(
    accuracy.percent_initial<accuracy.percent_improved, "better", "worse"))]

err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
roc.dt[, seg.i := 1:.N, by=.(set.fold, pred.name)]
roc.dt[, mid.thresh := (min.thresh+max.thresh)/2]
roc.tall <- melt(
  roc.dt,
  measure.vars=names(err.sizes))
status.colors <- c(
  same="black",
  better="blue",
  worse="red")
tallrect.dt <- data.table(
  mid.thresh=seq(-10, 5, by=0.2))
roc.dots <- roc.dt[tallrect.dt, .(
  set.fold, pred.name, mid.thresh=i.mid.thresh, FPR, TPR
), on=.(
  min.thresh<mid.thresh, max.thresh>mid.thresh)]
animint(
  title="Minimizing area under min(FP,FN)",
  out.dir="figure-auc-improved-interactive",
  ## ggplot()+
  ##   ggtitle("Data sets ordered by mean improved AUC")+
  ##   guides(size="none", color="none")+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "lines"))+
  ##   theme_animint(width=800)+
  ##   facet_grid(. ~ variable, scales="free")+
  ##   geom_segment(aes(
  ##     mean+sd, set.fac,
  ##     color=pred.name, size=pred.name,
  ##     xend=mean-sd, yend=set.fac),
  ##     data=auc.stats)+
  ##   scale_size_manual(values=c(improved=2, initial=3))+
  ##   geom_point(aes(
  ##     mean, set.fac,
  ##     color=pred.name),
  ##     shape=21,
  ##     size=3,
  ##     fill="white",
  ##     data=auc.stats)+
  ##   geom_point(aes(
  ##     value, set.fac,
  ##     color=pred.name),
  ##     clickSelects="set.fold",
  ##     alpha=0.6,
  ##     size=4,
  ##     data=auc.tall)+
  ##   xlab("")+
  ##   ylab("Data set"),
  ggplot()+
    ggtitle("FP/FN curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(update_axes="y")+
    facet_grid(pred.name ~ .)+
    xlab("Prediction threshold")+
    ylab("Incorrectly predicted labels")+
    geom_vline(aes(
      xintercept=(min.thresh+max.thresh)/2,
      key=1),
      data=auc.improved,
      color="grey",
      showSelected="set.fold")+
    geom_line(aes(
      mid.thresh, value,
      key=variable, group=variable,
      color=variable, size=variable),
      data=roc.tall,
      showSelected="set.fold")+
    geom_polygon(aes(
      mid.thresh, min.fp.fn, key=1),
      color="grey",
      data=roc.dt,
      size=0,
      showSelected="set.fold")+
    make_tallrect(tallrect.dt, "mid.thresh")+
    ## geom_tallrect(aes(
    ##   xmin=min.thresh, xmax=max.thresh),
    ##   data=roc.dt,
    ##   showSelected="set.fold",
    ##   clickSelects="mid.thresh",
    ##   alpha=0.5)+
    scale_size_manual(values=err.sizes)+
    scale_color_manual(values=err.colors),
  selector.types=list(
    variable="multiple"),
  ggplot()+
    ggtitle("ROC curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_path(aes(
      FPR, TPR,
      key=pred.name,
      color=pred.name, group=pred.name),
      data=roc.dt,
      showSelected="set.fold")+
    geom_point(aes(
      FPR, TPR, color=pred.name, key=pred.name),
      fill="white",
      data=auc.improved,
      showSelected="set.fold")+
    geom_point(aes(
      FPR, TPR, color=pred.name, key=pred.name),
      data=roc.dots,
      showSelected=c("set.fold", "mid.thresh"),
      size=4,
      alpha=0.5)+
    geom_point(aes(
      FPR, TPR, color=pred.name, key=paste(pred.name, mid.thresh)),
      data=roc.dots,
      showSelected="set.fold",
      clickSelects="mid.thresh",
      size=4,
      alpha=0.5),
  ggplot()+
    ggtitle("Percent correctly predicted labels")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    scale_color_manual(values=status.colors)+
    guides(color="none")+
    coord_equal()+
    geom_point(aes(
      accuracy.percent_initial, accuracy.percent_improved,
      key=set.fold,
      color=accuracy.percent_status),
      clickSelects="set.fold",
      alpha=0.6,
      size=4,
      data=roc.wide),
  ggplot()+
    ggtitle("Log[Area under Min(FP,FN) + 1]")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    guides(color="none")+
    scale_color_manual("Status",values=status.colors)+
    geom_point(aes(
      initial, improved,
      key=set.fold,
      color=status),
      clickSelects="set.fold",
      size=4,
      alpha=0.6,
      data=aum.wide)+
    coord_equal(),
  ggplot()+
    ggtitle("Area under the ROC curve")+
    theme_bw()+
    theme_animint(width=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_abline(slope=1, intercept=0, color="grey")+
    scale_color_manual("Status",values=status.colors)+
    geom_point(aes(
      auc_initial, auc_improved,
      key=set.fold,
      color=auc_status),
      clickSelects="set.fold",
      size=4,
      alpha=0.6,
      data=roc.wide)+
    coord_equal(),
  duration=list(
    set.fold=500,
    mid.thresh=500),
  time=list(
    ms=500,
    variable="mid.thresh"),
  first=list(
    set.fold="H3K27me3_RL_cancer/2",
    mid.thresh=-5
    )
)
