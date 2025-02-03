library(animint2)
library(data.table)

data(neuroblastomaProcessed, package="penaltyLearning")
data(neuroblastoma, package="neuroblastoma")
e <- function(label, profile.id, chromosome){
  data.table(label, profile.id=factor(profile.id), chromosome=factor(chromosome))
}
select.dt <- rbind(
  e("pos", 4, 2),
  e("neg", 513, 3))
nb.list <- lapply(neuroblastoma, data.table)
nb.some <- lapply(nb.list, "[", select.dt, on=.NATURAL)
max.segments <- max(neuroblastomaProcessed$errors$n.segments)
nb.segs <- nb.some$profiles[, {
  cum.vec <- cumsum(c(0, logratio))
  d <- diff(position)/2
  between <- position[-1]-d
  data.start.pos <- c(position[1]-d[1], between)
  data.end.pos <- c(between, position[.N]+d[.N-1])
  fit <- jointseg::Fpsn(logratio, max.segments)
  end.t <- t(fit$t.est)
  end.dt <- data.table(
    end=as.integer(end.t),
    segments=as.integer(col(end.t))
  )[!is.na(end)]
  end.dt[, start := c(0, end[-.N])+1, by=segments]
  end.dt[, mean := (cum.vec[end+1]-cum.vec[start])/(end-start+1)]
  end.dt[, `:=`(
    start.pos=data.start.pos[start],
    end.pos=data.end.pos[end]
  )]
}, by=label]

some.err <- neuroblastomaProcessed$errors[select.dt, .(
  segments=n.segments,
  fp, fn, possible.fp, possible.fn,
  min.log.lambda=-max.log.lambda,
  max.log.lambda=-min.log.lambda,
  errors, labels,
  label
), on=list(profile.id, chromosome)]
err.sizes <- c(
  "min(FP,FN)"=2,
  FP=6,
  FN=4)
err.colors <- c(
  correct="transparent",
  "min(FP,FN)"="black",
  FP="red",
  FN="deepskyblue")
some.err.tall <- melt(
  some.err,
  measure.vars=c("fp","fn"),
  variable.name="var.lower")
some.err.tall[, error.type := toupper(var.lower)]
leg <- "Error type"

dmin <- 4
dmax <- 6.5
some.err[, fp.diff := c(NA, diff(fp)), by=label]
some.err[, fn.diff := c(NA, diff(fn)), by=label]
some.diff <- some.err[fp.diff != 0 | fn.diff != 0, .(
  id=1, label, fp.diff, fn.diff, pred.log.lambda=min.log.lambda)]
some.diff[, fp.cum := cumsum(fp.diff), by=label]
some.diff[, fn.cum := rev(cumsum(rev(-fn.diff))), by=label]
dlist <- split(some.diff, some.diff[["label"]])
border.pred <- with(dlist, pos[ #orange dots
  neg,
  data.table(
    differentiable=FALSE,
    pos=pred.log.lambda,
    neg=i.pred.log.lambda),
  on="id",
  allow.cartesian=TRUE])
grid.pred <- data.table( #black dots
  differentiable=TRUE,
  pos=0,
  neg=seq(dmin, dmax, by=0.05))
both.pred <- rbind(border.pred, grid.pred)
both.pred[, pred.diff := neg-pos]
pred.tall <- melt(
  both.pred,
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda")[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[order(pred.diff)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "label")
  pos <- pred.log.lambda[label=="pos"]
  with(L, data.table(
    aum, auc,
    SM=roc[min.thresh < max.thresh, sum(min.fp.fn)],
    roc=list(roc[, `:=`(
      min.thresh=min.thresh+pos,
      max.thresh=max.thresh+pos
    )])
  ))
}, by=list(pred.diff, differentiable)]
metrics.wide[auc==max(auc)] #max auc => aum>0.
metrics.wide[14:15, roc ]

##compute slope and intercept of each of the 6 T_b(s) functions, plot
##them using geom_abline, and geom_point to represent the 9
##intersection points.
some.diff[, `:=`(
  slope=ifelse(label=="pos", 0, -1),
  intercept=pred.log.lambda-ifelse(label=="pos", 0, 6.5))]

##ignore rest.

show.roc.dt <- metrics.wide[, data.table(
  roc[[1]],
  AUC=auc, AUM=round(aum,3)
), by=pred.diff]
show.roc.tall <- melt(
  show.roc.dt,
  measure=c("fp","fn","min.fp.fn"),
  variable.name="lower.var")
show.roc.tall[, error.type := ifelse(
  lower.var=="min.fp.fn", "min(FP,FN)", toupper(lower.var))]

metrics.tall <- melt(
  metrics.wide,
  measure.vars=c("aum", "auc"),
  variable.name="var.lower"
)[order(-differentiable)]
metrics.tall[, variable := toupper(var.lower)]

show.roc.dt[, roc.point := rank(min.thresh), by=pred.diff]
thresh.offset <- 0.1
show.roc.dt[, text.constant := ifelse(
  min.thresh==-Inf, max.thresh-thresh.offset,
  ifelse(
    max.thresh==Inf,
    min.thresh+thresh.offset,
    (min.thresh+max.thresh)/2
  ))]
show.roc.dt[, text.roc.i := rank(roc.point), by=.(pred.diff, FPR, TPR)]
show.roc.dt[, text.FPR := text.roc.i*0.04+FPR]
text.size <- 15
text.color <- "blue"
both.pred.adj <- melt(both.pred[, .(
  differentiable,
  pos=0,
  neg=pred.diff,
  pred.diff
)],
measure.vars=select.dt$label,
variable.name = "label",
value.name="pred.log.lambda")
pred.tall.thresh <- both.pred.adj[
  show.roc.dt, on="pred.diff", allow.cartesian=TRUE]
pred.tall.thresh[, pred.plus.constant := pred.log.lambda+text.constant]
pred.tall.thresh.wide <- dcast(
  pred.tall.thresh,
  pred.diff + roc.point ~ label,
  value.var="pred.plus.constant"
)[, label := "neg"]
nb.models <- nb.segs[start==1, .(label, segments)]
nb.changes <- nb.segs[start>1]
err.list <- penaltyLearning::labelError(
  nb.models,
  nb.some$annotations,
  nb.changes,
  model.vars="segments",
  change.var="start.pos",
  problem.vars="label")
selected.dt <- pred.tall.thresh[
  some.err,
  data.table(label, pred.diff, roc.point, segments),
  nomatch=NULL,
  on=.(
    label,
    pred.plus.constant < max.log.lambda,
    pred.plus.constant > min.log.lambda
  )]
selected.segs <- nb.segs[selected.dt, on=.(
  label, segments), allow.cartesian=TRUE]
type.abbrev <- c(
  "false negative"="FN",
  "false positive"="FP",
  correct="correct")
selected.err <- err.list$label.errors[selected.dt, on=.(
  label, segments)][, error.type := type.abbrev[status] ]
viz <- animint(
  title="Simple non-monotonic ROC curve",
  video="https://vimeo.com/1053132517",
  out.dir="2021-11-12-aum-convexity",
  overview=ggplot()+
    ggtitle("Overview, select difference")+
    theme_bw()+
    theme(panel.margin=grid::unit(1, "lines"))+
    theme_animint(width=300, height=300)+
    facet_grid(variable ~ ., scales="free")+
    scale_fill_manual(values=c(
      "TRUE"="black",
      "FALSE"="orange"))+
    geom_point(aes(
      pred.diff, value, fill=differentiable),
      help="One dot per prediction difference which could be selected.",
      size=4,
      shape=21,
      data=metrics.tall)+
    make_tallrect(metrics.tall, "pred.diff")+ 
    xlab("Prediction difference, f(neg) - f(pos)")+
    coord_cartesian(xlim=c(dmin,dmax))+
    scale_y_continuous("", breaks=seq(0, 3, by=1)),
  data=ggplot()+
    ggtitle("Data, labels, predicted changepoint models")+
    theme_bw()+
    theme(legend.position="none")+
    theme_animint(width=600, height=300)+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6, fill=annotation),
      help="One rect per label. The negative label (orange/top) generates a false positive if any change-points are predicted inside. The positive label (violet/bottom) generates a false negative if no change-points are predicted inside.",
      alpha=0.5,
      data=nb.some$annotations)+
    scale_fill_manual(
      "label",
      values=c(
        breakpoint="violet",
        normal="orange"))+
    geom_point(aes(
      position/1e6, logratio),
      help="One dot per data point to segment.",
      color="grey50",
      data=nb.some$profiles)+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6,
      color=error.type),
      help="Rect border if predicted change-points generate a label error: grey=correct, red=false positive, blue=false negative.",
      data=selected.err,
      showSelected=c("pred.diff", "roc.point"),
      size=5,
      fill="transparent")+
    geom_segment(aes(
      start.pos/1e6, mean,
      xend=end.pos/1e6, yend=mean),
      help="Blue segments show predicted mean model.",
      data=selected.segs,
      size=3,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    geom_vline(aes(
      xintercept=start.pos/1e6),
      help="Blue vertical lines show predicted change-points.",
      data=selected.segs[start>1],
      size=2,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    scale_color_manual(leg,values=err.colors)+
    facet_grid(label ~ ., labeller=label_both)+
    scale_y_continuous(
      "DNA copy number (logratio)")+
    scale_x_continuous(
      "Position on chromosome"),
  obsErr=ggplot()+
    ggtitle("Example error functions")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme(legend.position="none")+
    theme_animint(width=300, height=300)+
    facet_grid(label ~ ., labeller=label_both)+
    geom_vline(aes(
      xintercept=pred.plus.constant),
      help="Black vertical line shows predicted value, f(x) = -log(penalty).",
      data=pred.tall.thresh,
      showSelected=c("pred.diff", "roc.point"))+
    geom_segment(aes(
      pos, -Inf,
      xend=neg, yend=-Inf),
      help="Black horizontal segment shows the difference in predicted values, between the two labeled data sequences (panels).",
      data=pred.tall.thresh.wide,
      showSelected=c("pred.diff", "roc.point"))+
    geom_text(aes(
      neg-0.1, -0.3,
      label=sprintf("pred.diff=%.2f", pred.diff)),
      help="Text shows the difference in predicted values, between the two labeled data sequences (panels).",
      hjust=1,
      data=pred.tall.thresh.wide,
      showSelected=c("pred.diff", "roc.point"))+
    geom_segment(aes(
      min.log.lambda, value,
      xend=max.log.lambda, yend=value,
      color=error.type, size=error.type),
      help="Blue and red segments show number of label errors, for each data sequence, as a function of predicted value f(x) = -log(penalty).",
      showSelected="error.type",
      data=some.err.tall)+
    scale_y_continuous(
      "Label errors",
      breaks=c(0,1),
      limits=c(-0.4, 1.4))+
    scale_color_manual(leg,values=err.colors)+
    scale_size_manual(leg,values=err.sizes)+
    scale_x_continuous(
      "Predicted value f(x)"),
  totals=ggplot()+
    ggtitle("Total error, select interval")+
    theme_bw()+
    theme(panel.grid.minor=element_blank())+
    theme_animint(width=300, height=300)+
    geom_rect(aes(
      xmin=min.thresh, xmax=max.thresh,
      ymin=0, ymax=min.fp.fn),
      fill="grey50",
      help="Grey rects represent AUM = Area Under Min of False Positive and False Negative rates.",
      color=NA,
      alpha=0.5,
      showSelected="pred.diff",
      show.roc.dt)+
    geom_segment(aes(
      min.thresh, value,
      xend=max.thresh, yend=value,
      color=error.type, size=error.type),
      showSelected="pred.diff",
      help="Blue/red/black segments show number of label errors, summed over all data sequences, as a function of the selected threshold, or constant c added to predicted values, f(x)+c = -log(penalty)+c.",
      data=show.roc.tall)+
    geom_vline(aes(
      xintercept=text.constant),
      showSelected=c("pred.diff", "roc.point"),
      color=text.color,
      help="Blue vertical line shows a constant that could be added to predicted values, which corresponds to the selected point on the ROC curve.",
      alpha=0.5,
      data=show.roc.dt)+
    geom_text(aes(
      text.constant, -0.25, label=roc.point),
      showSelected="pred.diff",
      help="One number for each interval of constant label error / point on the ROC curve.",
      size=text.size,
      color=text.color,
      data=show.roc.dt)+
    geom_text(aes(
      -1.5, 0.25, label=sprintf("AUM=%.2f", aum)),
      help="Text shows AUM = Area Under Min of False Positive and False Negative rates.",
      data=metrics.wide,
      showSelected="pred.diff")+
    geom_tallrect(aes(
      xmin=min.thresh, xmax=max.thresh),
      help="Light grey rect represents interval of constants that could be added to predicted values, to obtain the selected point on the ROC curve.",
      data=show.roc.dt,
      fill=text.color,
      clickSelects="roc.point",
      showSelected="pred.diff",
      color="transparent",
      alpha=0.1)+
    scale_y_continuous(
      "Label errors",
      breaks=c(0,1),
      limits=c(-0.4, 1.4))+
    scale_color_manual(leg,values=err.colors)+
    scale_size_manual(leg,values=err.sizes)+
    ## geom_blank(aes(
    ##   x, y),
    ##   help="Blank geom for enlarging 
    ##   data=data.table(x=0, y=c(-0.4,1.4)))+
    scale_x_continuous(
      "Constant added to pred. values"),
  roc=ggplot()+
    ggtitle("ROC curve, select point")+
    theme_bw()+
    theme(panel.grid.minor=element_blank())+
    theme_animint(width=300, height=300)+
    geom_path(aes(
      FPR, TPR),
      help="Black path represents ROC curve.",
      showSelected="pred.diff",
      data=show.roc.dt)+
    geom_text(aes(
      0.5, 0.5, label=paste0("AUC=", auc)),
      help="Text shows AUC = Area Under the ROC Curve.",
      data=metrics.wide,
      showSelected="pred.diff")+
    scale_x_continuous(
      "False Positive Rate",
      breaks=seq(0,1,by=0.5))+
    scale_y_continuous(
      "True Positive Rate",
      breaks=seq(0,1,by=0.5))+
    geom_point(aes(
      FPR, TPR),
      help="One blue dot per point on the ROC curve.",
      data=show.roc.dt,
      size=4,
      alpha=0.5,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    geom_text(aes(
      text.FPR, TPR+0.01, label=roc.point),
      help="One number per point on the ROC curve, lower numbers for points generated by smaller constants c added to predicted values, f(x)+c = -log(penalty)+c.",
      size=text.size,
      color=text.color,
      showSelected="pred.diff",
      clickSelects="roc.point",
      data=show.roc.dt),
  time=list(
    variable="pred.diff",
    ms=500),
  source="https://github.com/tdhock/max-generalized-auc/blob/master/figure-aum-convexity-interactive.R"
)

viz
if(FALSE){
  animint2pages(viz, "2025-02-03-aum-convexity")
}

