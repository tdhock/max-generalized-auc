### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,...){
  local.lib <- file.path(getwd(), "library")
  dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
  .libPaths(c(local.lib, .libPaths()))
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      print(search())
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
options(repos=c(
  "http://www.bioconductor.org/packages/release/bioc",
  ##"http://r-forge.r-project.org",
  "http://cloud.r-project.org",
  "http://cran.r-project.org"))
works_with_R(
  "4.1.0",
  data.table="1.14.0",
  future="1.21.0",
  future.apply="1.7.0",
  RJSONIO="1.3.1.4",
  R.utils="2.10.1",
  "tdhock/penaltyLearning@4e14a0b0e022d919884277d68b8e47bd158459f3",
  jointseg="1.0.2",
  gridExtra="2.3",
  neuroblastoma="1.0",
  tikzDevice="0.12.3.1",
  microbenchmark="1.4.7",
  animint2="1.0")

curveAlignment <- readRDS("curveAlignment.rds")
AUC.dt.list <- list()
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
    stop("infinite AUM")
  }
  AUM <- min.positive[, {
    sum(min.fp.fn*width.thresh)
  }]
  AUC.dt.list[[paste(offset)]] <- with(roc, data.table(
    AUC=auc, AUM, offset,
    AUM.deriv=offset.deriv$deriv,
    min.errors=off.min$errors[1],
    n.min=nrow(off.min),
    thresholds[threshold=="min.error"]))
  roc.dt.list[[paste(offset)]] <- data.table(
    offset, roc$roc, piece=1:nrow(roc$roc), prob.dir="Total")
}
AUC.dt <- do.call(rbind, AUC.dt.list)
roc.dt <- do.call(rbind, roc.dt.list)
roc.segs <- do.call(rbind, roc.segs.list)
err.dt <- do.call(rbind, err.dt.list)
roc.win.err.dt <- do.call(rbind, roc.win.err.list)

ggplot()+
  geom_point(aes(
    offset, AUC),
    data=AUC.dt)

common.names <- intersect(names(roc.dt), names(err.dt))
both.dt <- rbind(
  err.dt[, common.names, with=FALSE],
  roc.dt[, common.names, with=FALSE])
both.dt[, `min(fp,fn)` := pmin(fp, fn)]
err.dt.tall <- melt(
  both.dt,
  variable.name="error.type",
  measure.vars=c("fp", "fn", "errors", "min(fp,fn)"))
id2show <- function(seqID)gsub(
  "ATAC_JV_adipose/samples/AC1/|/problems/chrX-37148256-49242997", "",
  seqID)
roc.segs[, sample := id2show(prob.dir)]
curveAlignment$labels[, sample := id2show(prob.dir)]
roc.win.err.dt[, sample := id2show(prob.dir)]
curveAlignment$profiles[, sample := id2show(prob.dir)]
err.dt.tall[, sample := id2show(prob.dir)]
AUC.dt[, thresh := (min.thresh+max.thresh)/2]
roc.dt[, Errors := ifelse(errors==min(errors), "Min", "More"), by=list(offset)]
AUC.dt[, max.correct := as.numeric(labels-min.errors)]
AUC.dt[, opt.models := as.numeric(n.min)]
AUC.tall <- melt(
  AUC.dt,
  measure.vars=c("AUM", "AUM.deriv", "AUC", "min.errors", "opt.models"))
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
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
err.dt.show <- err.dt.tall[
  error.type %in% c("fp","fn") |
    (error.type=="min(fp,fn)" & sample=="Total")]
area.show <- err.dt.show[error.type=="min(fp,fn)"]
AUM.text <- area.show[, .SD[value>0][1], by=offset][AUC.dt, .(
  offset, sample, min.thresh, AUM), on="offset"]
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
        x=0, y=c(3.4, 1.6),
        variable="min.errors"))+
    xlab("Offset = Difference between predicted values of samples")+
    geom_point(aes(
      offset, value),
      fill=NA,
      data=AUC.tall)+
    ylab("")+
    make_tallrect(AUC.dt, "offset"),
  error=ggplot()+
    ggtitle(
      "Error curves, select threshold")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample ~ ., scales="free")+
    scale_color_manual(values=c(
      fp="red",
      fn="deepskyblue",
      "min(fp,fn)"="black"))+
    scale_size_manual(values=c(
      fp=5,
      fn=3,
      "min(fp,fn)"=1))+
    xlab("Prediction threshold")+
    scale_y_continuous(
      "Number of incorrectly predicted labels",
      breaks=seq(0, 20, by=2))+
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
        "Min Error=", errors)),
      showSelected=c("offset", "Errors"),
      hjust=0,
      color="grey50",
      data=data.table(
        AUC.dt,
        Errors="Min",
        sample="Total"))+
    geom_rect(aes(
      xmin=min.thresh, ymin=0,
      key=piece,
      xmax=max.thresh, ymax=value),
      chunk_vars=chunk_vars,
      showSelected="offset",
      fill="black",
      data=area.show)+
    geom_segment(aes(
      min.thresh, value,
      key=paste(piece, error.type),
      color=error.type,
      size=error.type,
      xend=max.thresh, yend=value),
      chunk_vars=chunk_vars,
      showSelected="offset",
      data=err.dt.show)+
    geom_text(aes(
      min.thresh, 0.5, key=1,
      label=sprintf("AUM=%.1f", AUM)),
      showSelected="offset",
      hjust=1,
      data=AUM.text)+
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
        "AUC=%.2f", AUC)),
      showSelected="offset",
      data=AUC.dt)+
    geom_abline(aes(
      slope=slope, intercept=intercept),
      color="grey",
      data=data.table(slope=1, intercept=0))
)

