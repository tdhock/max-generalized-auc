library(data.table)
library(animint2)
nb.comb <- readRDS("neuroblastomaProcessed.combinations.rds")
roc.dt <- nb.comb$auc[, data.table(
  roc[[1]]
), by=.(size, combo.i)]
perfect <- roc.dt[FPR==0 & TPR==1]
nb.comb$auc[!perfect, on=.(size, combo.i)]

one.combo <- roc.dt[combo.i==256]
ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=one.combo)+
  facet_wrap("size")
wrong.way <- roc.dt[, {
  diff.dt <- data.table(
    dtpr=diff(TPR),
    dfpr=diff(FPR))
  diff.dt[, .(
    TPR=sum(dtpr[dtpr<0]),
    FPR=sum(dfpr[dfpr<0])
    )]
}, by=.(size, combo.i)]
u.roc <- roc.dt[, {
  list(n.uniq=nrow(unique(data.table(FPR, TPR))))
}, by=.(size, combo.i)]
auc.stats <- nb.comb$auc[u.roc, .(
  size, combo.i, auc, aum=aub, n.finite, n.uniq,
  panel.key=paste("uniq", n.uniq, auc)
), on=.(size, combo.i)]

auc.stat.counts <- auc.stats[, .(
  count=.N
), by=.(size, n.uniq, panel.key, auc)]
ggplot()+
  facet_wrap("size")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_hline(yintercept=1, col="grey")+
  scale_fill_gradient(low="white", high="red")+
  scale_y_continuous(breaks=seq(0, 1.2, by=0.2))+
  scale_x_continuous(breaks=seq(0, 20, by=2))+
  geom_point(aes(
    n.uniq, auc, fill=count),
    shape=21,
    size=4,
    data=auc.stat.counts)

## ideas for manual histogram.
max.finite <- max(auc.stats$n.finite)
break.vec <- seq(0, 10, by=0.5)
edge.vec <- seq(
  0, max.finite, l=11)
rect.dt <- data.table(
  xmin=edge.vec[-length(edge.vec)],
  xmax=edge.vec[-1])
##link two plots, show details.
auc.stats[, round.auc := round(auc, 4)]
rfac <- 5
auc.stats[, round.aum := round(aum*rfac)/rfac]
auc.stats[, prop.finite := n.finite/max.finite]
panel.titles <- c(
  round.aum="Area under both",
  prop.finite="Prop. predictions in finite interval",
  round.auc="Area under ROC curve")
xlevs <- c("prop.finite", "round.aum")
ylevs <- c("prop.finite", "round.auc")
panel.dt <- data.table(
  expand.grid(
    xvar=xlevs,
    yvar=ylevs,
    stringsAsFactors=FALSE)
)[xvar!=yvar]
auc.panels <- panel.dt[, {
  auc.stats[, {
    getvar <- function(v)as.numeric(.SD[[v]])
    dt <- data.table(size, combo.i)
    for(xy in c("x", "y")){
      v <- get(paste0(xy, "var"))
      L <- list(
        val=v,
        orig=sub("round.", "", v))
      for(suffix in names(L)){
        to.col <- paste0(xy, suffix)
        from.col <- L[[suffix]]
        dt[[to.col]] <- as.numeric(.SD[[from.col]])
      }
    }
    dt
  }]
}, by=.(xvar, yvar)]
xfac <- function(val)factor(val, xlevs, panel.titles[xlevs])
yfac <- function(val)factor(val, ylevs, panel.titles[ylevs])
auc.panels[, xfac := xfac(xvar)]
auc.panels[, yfac := yfac(yvar)]
auc.panels[, key := paste(xval, yval)]
auc.panels[, panel.key := paste(xvar, yvar, key)]
count.dt <- auc.panels[, .(
  combos=.N
), by=.(size, xfac, yfac, xval, yval, key, panel.key)]
ggplot()+facet_wrap("size")+geom_bar(aes(combos), data=count.dt)

ggplot()+facet_wrap("size")+geom_bar(aes(log10(combos)), data=count.dt)

combos.for.panel.key <- count.dt[auc.panels, .(
  size, panel.key, combo.i
), on=.(
  size, xfac, yfac, xval, yval, key, panel.key)]
## use first.orig to avoid overplotting clickable points.
first.orig <- auc.panels[, .SD[1], by=.(size, xfac, yfac, xorig, yorig)]
auc.orig <- combos.for.panel.key[auc.panels, .(
  panel.key, xfac=i.xfac, yfac=i.yfac, size,
  combo.i,
  xorig,
  yorig
), allow.cartesian=TRUE, on=.(size, combo.i)]
## for each panel.key, find all original data points with combo.i values.
XPANEL <- function(val, ...){
  data.table(..., xfac=xfac(val))
}
YPANEL <- function(val, ...){
  data.table(..., yfac=yfac(val))
}
count.wide <- dcast(
  count.dt,
  xfac + yfac + xval + yval + key + panel.key ~ size,
  value.var="combos")
count.tall <- melt(
  count.wide,
  id.vars=c("xfac", "yfac", "xval", "yval", "key", "panel.key"),
  variable.name="size",
  value.name="combos")
count.tall[is.na(combos), combos := 0]
count.tall[, combos.chr := ifelse(combos==0, "none", "some")]
roc.color <- "deepskyblue"
small.size <- 4
small.alpha <- 0.7
med.size <- 6
big.size <- 7
err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
roc.dt[, thresh.i := 1:.N, by=.(size, combo.i)]
roc.dt[, min.fp.fn := ifelse(fp<fn, fp, fn)]
roc.dt[, width.new := ifelse(
  min.fp.fn>0, width.thresh, ifelse(
    width.thresh != 0.1, 0.1))]
roc.dt[, max.new := cumsum(width.new), by=.(size, combo.i)]
roc.dt[, min.new := c(0, max.new[-.N]), by=.(size, combo.i)]
roc.dt.tall <- melt(
  roc.dt,
  measure.vars=names(err.colors))
aum0 <- auc.stats[, .(
  aum0=sum(aum==0),
  auc1=sum(auc==1),
  auc.over1=sum(auc>1)
), by=.(size)]
aum0.tall <- melt(aum0, id.vars="size")
viz <- animint(
  title="Generalized ROC curve metrics",
  source="https://github.com/tdhock/max-generalized-auc/blob/master/figure-neuroblastomaProcessed-combinations-interactive.R",
  out.dir="figure-neuroblastomaProcessed-combinations-interactive",
  sizes=ggplot()+
    ggtitle("Select margin size")+
    theme_bw()+
    theme_animint(width=300, height=300)+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(variable ~ ., scales="free")+
    scale_x_continuous(
      "Margin size of prediction wrt inf. interval")+
    scale_y_continuous("Number of prediction combinations")+
    geom_point(aes(
      log10(size), value),
      data=aum0.tall)+
    geom_tallrect(aes(
      xmin=log10(size)-0.5,
      xmax=log10(size)+0.5),
      data=aum0,
      alpha=0.5,
      clickSelects="size"),
  scatter=ggplot()+
    ggtitle("ROC curve, AUC/AUM distribution, select prediction")+
    theme_bw()+
    theme_animint(width=500, height=500)+
    theme(panel.margin=grid::unit(0.5, "lines"))+
    guides(color="none")+
    facet_grid(yfac ~ xfac, scales="free")+
    scale_x_continuous("", breaks=break.vec)+
    scale_y_continuous("", breaks=break.vec)+
    scale_color_manual(values=c(none="white", some="black"))+
    scale_fill_gradient(low="white", high="red")+
    geom_path(aes(
      FPR, TPR, key=1),
      data=XPANEL("prop.finite", YPANEL("prop.finite", roc.dt)),
      color=roc.color,
      showSelected=c("size", "combo.i"))+
    geom_point(aes(
      FPR, TPR, key=1),
      data=XPANEL("prop.finite", YPANEL("prop.finite", roc.dt)),
      showSelected=c("size", "combo.i", "thresh.i"))+
    geom_hline(aes(
      yintercept=auc),
      data=YPANEL("round.auc", auc=1),
      color="grey50")+
    geom_vline(aes(
      xintercept=aum),
      data=XPANEL("round.aum", aum=0),
      color="grey50")+
    geom_point(aes(
      xval, yval,
      key=1),
      size=big.size,
      showSelected=c("size", "panel.key"),
      data=count.tall)+
    geom_point(aes(
      xval, yval,
      key=key,
      tooltip=sprintf(
        "%d combos with %s=%s, %s=%s for size %s",
        combos,
        xfac, paste(xval),
        yfac, paste(yval),
        paste(size)),
      color=combos.chr,
      fill=combos),
      size=med.size,
      clickSelects="panel.key",
      showSelected=c("size"),
      data=count.tall)+
    geom_point(aes(
      xorig, yorig, key=combo.i),
      showSelected=c("size", "panel.key"),
      clickSelects="combo.i",
      alpha=small.alpha,
      size=small.size,
      data=auc.orig)+
    geom_point(aes(
      xorig, yorig, key=1),
      color=roc.color,
      showSelected=c("size", "combo.i"),
      size=small.size,
      data=auc.orig)+
    geom_blank(aes(
      x, y),
      data=data.table(x=0, y=0)),
  err=ggplot()+
    ggtitle("Error curves, select threshold")+
    theme_bw()+
    theme_animint(
      width=300, height=300
    )+
    theme(panel.margin=grid::unit(0, "lines"))+
    scale_x_continuous("Distorted prediction threshold")+
    scale_y_continuous("Incorrectly predicted labels")+
    scale_color_manual(values=err.colors)+
    scale_size_manual(values=err.sizes)+
    geom_rect(aes(
      xmin=min.new, ymin=0,
      key=thresh.i,
      xmax=max.new, ymax=min.fp.fn),
      alpha=0.1,
      size=0,
      showSelected=c("size", "combo.i"),
      data=roc.dt)+
    geom_segment(aes(
      min.new, value,
      key=paste(variable, thresh.i),
      color=variable, size=variable,
      xend=max.new, yend=value),
      data=roc.dt.tall,
      showSelected=c("size", "combo.i"))+
    geom_rect(aes(
      xmin=min.new, ymin=-Inf,
      xmax=max.new, ymax=Inf,
      key=thresh.i),
      clickSelects="thresh.i",
      showSelected=c("size", "combo.i"),
      alpha=0.5,
      data=roc.dt),
    ## geom_point(aes(
    ##   thresh.i, value),
    ##   data=roc.dt.tall,
    ##   showSelected=c("size", "combo.i")),
  ## uniq=ggplot()+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "lines"))+
  ##   geom_hline(yintercept=1, col="grey")+
  ##   scale_fill_gradient(low="white", high="red")+
  ##   scale_y_continuous(breaks=seq(0, 1.2, by=0.2))+
  ##   scale_x_continuous(breaks=seq(0, 20, by=2))+
  ##   geom_point(aes(
  ##     n.uniq, auc,
  ##     key=1),
  ##     showSelected=c("size", "panel.key"),
  ##     size=big.size,
  ##     data=auc.stat.counts)+
  ##   geom_point(aes(
  ##     n.uniq, auc,
  ##     key=panel.key,
  ##     fill=count),
  ##     showSelected="size",
  ##     clickSelects="panel.key",
  ##     size=med.size,
  ##     data=auc.stat.counts)+
  ##   geom_point(aes(
  ##     n.uniq, auc,
  ##     key=combo.i),
  ##     showSelected=c("size", "panel.key"),
  ##     clickSelects="combo.i",
  ##     size=small.size,
  ##     alpha=small.alpha,
  ##     data=auc.stats)+
  ##   geom_point(aes(
  ##     n.uniq, auc,
  ##     key=1),
  ##     showSelected=c("size", "combo.i"),
  ##     color=roc.color,
  ##     size=small.size,
  ##     data=auc.stats),
  duration=list(size=1000, combo.i=1000, thresh.i=500),
  first=list(
    thresh.i=1,
    size=1e-4,
    combo.i=132,
    panel.key="prop.finite round.auc 0.625 1.1667"),
  time=list(variable="thresh.i", ms=500)
)
if(FALSE){
  animint2pages(viz, "2024-06-26-neuroblastomaProcessed-combinations")
}
