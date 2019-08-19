source("packages.R")

curveAlignment <- readRDS("curveAlignment.rds")

auc.dt.list <- list()
roc.dt.list <- list()
err.dt.list <- list()
roc.segs.list <- list()
roc.win.err.list <- list()
off.by <- 0.2
offset.prob <- curveAlignment$problems$prob.dir[2]
for(offset in seq(-5, 5, by=off.by)){
  print(offset)
  pred.dt <- data.table(
    curveAlignment$problems, pred.log.lambda=10+c(0, offset))
  pred.eval <- curveAlignment$evaluation[pred.dt, on=list(prob.dir)]
  pred.eval[, possible.fn := possible.tp]
  roc <- penaltyLearning::ROChange(
    pred.eval, pred.dt, "prob.dir")
  ## compute derivative of Area under min(FP, FN).
  thresh.dt <- pred.eval[order(-min.log.lambda), {
    fp.diff <- diff(fp)
    fp.change <- fp.diff != 0
    fn.diff <- diff(fn)
    fn.change <- fn.diff != 0
    fp.dt <- if(any(fp.change))data.table(
      log.lambda=min.log.lambda[c(fp.change, FALSE)],
      fp=as.numeric(fp.diff[fp.change]),
      fn=0)
    fn.dt <- if(any(fn.change))data.table(
      log.lambda=min.log.lambda[c(fn.change, FALSE)],
      fp=0,
      fn=as.numeric(fn.diff[fn.change]))
    ##browser(expr=sample.id=="McGill0322")
    rbind(fp.dt, fn.dt)
  }, by=.(prob.dir)]
  pred.with.thresh <- thresh.dt[pred.dt, on=.(prob.dir), nomatch=0L]
  pred.with.thresh[, thresh := log.lambda - pred.log.lambda]
  first.dt <- pred.eval[max.log.lambda==Inf]
  thresh.ord <- pred.with.thresh[order(-thresh), .(
    prob.dir=c(NA, prob.dir),
    min.thresh=c(-Inf, log.lambda),
    max.thresh=c(log.lambda, Inf),
    fp = cumsum(c(sum(first.dt$fp), fp)),
    fn = cumsum(c(sum(first.dt$fn), fn)),
    change=c(0, ifelse(fp==0, fn, fp))
    )]
  thresh.ord[, min.fp.fn := ifelse(fp<fn, fp, fn)]
  thresh.ord[, min.change := c(NA, diff(min.fp.fn))]
  prob.deriv <- thresh.ord[min.change==change, .(
    deriv=-sum(change)
  ), by=.(prob.dir)]
  offset.deriv <- prob.deriv[offset.prob, on=.(prob.dir)]
  ## save info for display.
  pred.eval[, min.thresh := min.log.lambda-pred.log.lambda]
  pred.eval[, max.thresh := max.log.lambda-pred.log.lambda]
  pred.eval[, piece := 1:.N]
  err.dt.list[[paste(offset)]] <- data.table(offset, pred.eval)
  roc$roc[, thresh := (min.thresh+max.thresh)/2]
  pred.some.cols <- pred.dt[, list(id=1, prob.dir, pred.log.lambda)]
  roc.off.id <- data.table(offset, id=1, roc$roc)
  roc.off <- roc.off.id[pred.some.cols, on=list(
    id), allow.cartesian=TRUE]
  roc.off[, log.lambda := thresh + pred.log.lambda]
  roc.segs.list[[paste(offset)]] <-
    curveAlignment$segments[roc.off, nomatch=0L, allow.cartesian=TRUE, on=list(
      prob.dir,
      min.log.lambda<=log.lambda,
      max.log.lambda>=log.lambda)]
  roc.win.err.list[[paste(offset)]] <-
    curveAlignment$errors[roc.off, nomatch=0L, on=list(
      prob.dir,
      min.log.lambda<=log.lambda,
      max.log.lambda>=log.lambda)]
  off.min <- roc$roc[errors==min(errors)]
  roc$roc[, min.fp.fn := ifelse(fp<fn, fp, fn)]
  roc$roc[, width.thresh := max.thresh-min.thresh]
  roc$roc[, min.change.after := c(diff(min.fp.fn), NA)]
  min.positive <- roc$roc[0<min.fp.fn]
  bad <- min.positive[width.thresh==Inf]
  if(nrow(bad)){
    print(bad)
    stop("infinite aub")
  }
  aub <- min.positive[, {
    sum(min.fp.fn*width.thresh)
  }]
  auc.dt.list[[paste(offset)]] <- with(roc, data.table(
    auc, aub, offset,
    aub.deriv=offset.deriv$deriv,
    min.errors=off.min$errors[1],
    n.min=nrow(off.min),
    thresholds[threshold=="min.error"]))
  roc.dt.list[[paste(offset)]] <- data.table(
    offset, roc$roc, piece=1:nrow(roc$roc), prob.dir="Total")
}
auc.dt <- do.call(rbind, auc.dt.list)
roc.dt <- do.call(rbind, roc.dt.list)
roc.segs <- do.call(rbind, roc.segs.list)
err.dt <- do.call(rbind, err.dt.list)
roc.win.err.dt <- do.call(rbind, roc.win.err.list)

ggplot()+
  geom_point(aes(
    offset, auc),
    data=auc.dt)

common.names <- intersect(names(roc.dt), names(err.dt))
both.dt <- rbind(
  err.dt[, common.names, with=FALSE],
  roc.dt[, common.names, with=FALSE])
err.dt.tall <- melt(
  both.dt,
  variable.name="error.type",
  measure.vars=c("fp", "fn", "errors"))
id2show <- function(seqID)gsub(
  "ATAC_JV_adipose/samples/AC1/|/problems/chrX-37148256-49242997", "",
  seqID)
roc.segs[, sample := id2show(prob.dir)]
curveAlignment$labels[, sample := id2show(prob.dir)]
roc.win.err.dt[, sample := id2show(prob.dir)]
curveAlignment$profiles[, sample := id2show(prob.dir)]
err.dt.tall[, sample := id2show(prob.dir)]
auc.dt[, thresh := (min.thresh+max.thresh)/2]
roc.dt[, Errors := ifelse(errors==min(errors), "Min", "More"), by=list(offset)]
auc.dt[, max.correct := as.numeric(labels-min.errors)]
auc.dt[, opt.models := as.numeric(n.min)]
auc.tall <- melt(
  auc.dt,
  measure.vars=c("aub", "aub.deriv", "auc", "max.correct", "opt.models"))
min.err <- roc.dt[Errors=="Min"]
min.err[, piece := 1:.N, by=list(offset)]
roc.size <- 5
roc.peaks <- roc.segs[status=="peak"]
text.y <- c(
  "offset"=200,
  "thresh"=175,
  "fp"=150,
  "fn"=125)
text.dt <- melt(
  roc.dt,
  id.vars=c("offset", "thresh"),
  measure.vars=names(text.y))
text.dt[, y := text.y[variable] ]
text.dt[, digits := ifelse(variable %in% c("fp", "fn"), 0, 1)]
text.dt[, value.num := round(value, digits)]
text.dt[, value.str := paste(value.num)]
err.dt.tall[, value.i := cumsum(
  c(FALSE, diff(value) != 0)
), by=list(sample, offset, error.type)]
err.dt.segs <- err.dt.tall[, list(
  min.thresh=min(min.thresh),
  max.thresh=max(max.thresh),
  value=value[1]
), by=list(sample, offset, error.type, value.i)]
min.tallrects <- data.table(
  err.dt.segs[error.type=="errors", {
    .SD[value==min(value)]
  }, by=list(sample, offset)],
  Errors="Min")
chunk_vars <- "offset"
chunk_vars <- character()
animint(
  title="Changepoint detection ROC curve alignment problem",
  ##first=list(offset=0.5),
  out.dir="figure-curveAlignment",
  duration=list(offset=250),
  time=list(variable="offset", ms=250),
  profiles=ggplot()+
    ylab("Number of aligned DNA sequence reads (coverage)")+
    ggtitle(
      "Noisy coverage data, labels, and predicted model")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=1300)+
    facet_grid(sample ~ window, scales="free", labeller=label_both)+
    ## geom_text(aes(
    ##   43447, 200, label=sprintf(
    ##     "offset=%.1f thresh=%.1f fp=%d fn=%d",
    ##     offset, thresh, fp, fn)),
    ##   hjust=0,
    ##   showSelected=c("offset", "thresh"),
    ##   data=data.table(
    ##     roc.dt,
    ##     window=1,
    ##     sample="MSC83"))+
    geom_text(aes(
      43447, y,
      key=paste(offset, thresh, variable),
      label=paste0(
        variable, "=", value.str)),
      hjust=0,
      showSelected=c("offset", "thresh"),
      data=data.table(
        text.dt,
        window=1,
        sample="MSC83"))+
    geom_tallrect(aes(
      xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
      data=curveAlignment$labels,
      alpha=0.5,
      color="grey")+
    scale_linetype_manual(
      "Error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      key=paste(chromStart, chromEnd),
      linetype=status),
      data=roc.win.err.dt,
      chunk_vars=chunk_vars,
      showSelected=c("offset", "thresh"),
      fill=NA,
      size=2,
      color="black")+
    scale_fill_manual(values=ann.colors)+
    geom_step(aes(
      chromStart/1e3, coverage),
      data=curveAlignment$profiles,
      color="grey50")+
    geom_segment(aes(
      segStart/1e3, mean,
      key=paste(segStart, segEnd),
      xend=segEnd/1e3, yend=mean),
      color="green",
      alpha=0.7,
      chunk_vars=chunk_vars,
      showSelected=c("offset", "thresh"),
      data=roc.segs)+
    geom_segment(aes(
      segStart/1e3, 0,
      key=paste(segStart, segEnd),
      xend=segEnd/1e3, yend=0),
      color="deepskyblue",
      showSelected=c("offset", "thresh"),
      chunk_vars=chunk_vars,
      size=3,
      alpha=0.7,
      data=roc.peaks)+
    geom_point(aes(
      segStart/1e3, 0,
      key=paste(segStart, segEnd)),
      color="deepskyblue",
      showSelected=c("offset", "thresh"),
      chunk_vars=chunk_vars,
      size=4,
      fill="white",
      alpha=0.7,
      data=roc.peaks)+
    scale_x_continuous(
      "Position on chrX (kb = kilo bases, reference genome hg19)",
      breaks=seq(4e4, 5e4, by=5)),
  metrics=ggplot()+
    ggtitle(
      "AUC, select offset")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(variable ~ ., scales="free")+
    geom_blank(aes(
      x, y),
      data=data.table(
        x=0, y=c(10.4, 8.6),
        variable="max.correct"))+
    xlab("Offset = Difference between predicted values of samples")+
    geom_point(aes(
      offset, value),
      fill=NA,
      data=auc.tall)+
    ylab("")+
    make_tallrect(auc.dt, "offset"),
  ## auc=ggplot()+
  ##   ggtitle(
  ##     "AUC, select offset")+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "lines"))+
  ##   geom_point(aes(
  ##     offset, auc),
  ##     fill=NA,
  ##     data=auc.dt)+
  ##   geom_tallrect(aes(
  ##     xmin=offset-off.by/2, xmax=offset+off.by/2),
  ##     clickSelects="offset",
  ##     alpha=0.5,
  ##     data=auc.dt),
  error=ggplot()+
    ggtitle(
      "Error curves, select threshold")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample ~ ., scales="free")+
    scale_color_manual(values=c(
      fp="red",
      fn="deepskyblue",
      errors="black"))+
    scale_size_manual(values=c(
      fp=5,
      fn=3,
      errors=1))+
    xlab("Prediction threshold")+
    scale_y_continuous(
      "Number of incorrectly predicted labels",
      breaks=seq(0, 20, by=2))+
    ## geom_vline(aes(
    ##   xintercept=xintercept),
    ##   color="grey",
    ##   data=data.table(xintercept=0))+
    ## geom_vline(aes(
    ##   xintercept=thresh, key=piece),
    ##   showSelected=c("offset", "Errors"),
    ##   color="grey50",
    ##   data=data.table(
    ##     min.err,
    ##     Errors="Min"))+
    ## geom_text(aes(
    ##   thresh+0.2, labels*0.9, key=1, label=paste0(
    ##     "Min Errors=", errors)),
    ##   showSelected=c("offset", "Errors"),
    ##   hjust=0,
    ##   color="grey50",
    ##   data=data.table(
    ##     auc.dt,
    ##     Errors="Min",
    ##     sample="Total"))+
    geom_tallrect(aes(
      xmin=min.thresh,
      xmax=max.thresh,
      key=min.thresh),
      showSelected=c("offset", "Errors"),
      color="grey50",
      alpha=0.5,
      data=min.tallrects)+
    geom_text(aes(
      min.thresh, labels*0.9, key=1, label=paste0(
        "Min Errors=", errors)),
      showSelected=c("offset", "Errors"),
      hjust=0,
      color="grey50",
      data=data.table(
        auc.dt,
        Errors="Min",
        sample="Total"))+
    geom_segment(aes(
      min.thresh, value,
      key=paste(piece, error.type),
      color=error.type,
      size=error.type,
      xend=max.thresh, yend=value),
      chunk_vars=chunk_vars,
      showSelected="offset",
      data=err.dt.tall)+
    geom_tallrect(aes(
      xmin=min.thresh, xmax=max.thresh,
      tooltip=sprintf(
        "%.1f<thresh<%.1f FP=%d FN=%d",
        min.thresh, max.thresh, fp, fn),
      key=paste(offset, thresh)),
      showSelected="offset",
      clickSelects="thresh",
      alpha=0.5,
      data=roc.dt),
  roc=ggplot()+
    ggtitle(
      "ROC curves, select threshold")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_path(aes(
      FPR, TPR, key=paste(offset, thresh)),
      showSelected="offset",
      data=roc.dt)+
    ## geom_point(aes(
    ##   FPR, TPR, key=paste(FPR, TPR)),
    ##   showSelected="offset",
    ##   data=min.err,
    ##   size=roc.size,
    ##   fill=NA)+
    geom_point(aes(
      FPR, TPR, fill=Errors,
      tooltip=sprintf(
        "%.1f<thresh<%.1f FP=%d FN=%d",
        min.thresh, max.thresh, fp, fn),
      key=paste(offset, thresh)),
      showSelected="offset",
      clickSelects="thresh",
      size=roc.size,
      alpha=0.7,
      data=roc.dt)+
    scale_fill_manual(values=c(
      Min="black",
      More="white"))+
    coord_equal()+
    geom_text(aes(
      0.75, 0.25, key=1, label=sprintf(
        "AUC=%.2f", auc)),
      showSelected="offset",
      data=auc.dt)+
    geom_abline(aes(
      slope=slope, intercept=intercept),
      color="grey",
      data=data.table(slope=1, intercept=0))
)

