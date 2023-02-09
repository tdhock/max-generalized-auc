library(animint2)
library(data.table)

data(neuroblastomaProcessed, package="penaltyLearning")
data(neuroblastoma, package="neuroblastoma")
e <- function(label, profile.id, chromosome){
  data.table(
    label, 
    profile.id=factor(profile.id), 
    chromosome=factor(chromosome))
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
  profile.id, chromosome,
  segments=n.segments,
  fp, fn, possible.fp, possible.fn,
  min.log.lambda=-max.log.lambda,
  max.log.lambda=-min.log.lambda,
  min.lambda,
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

dmin <- 2
dmax <- 4
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
neg.seq <- seq(dmin, dmax, by=0.05)
grid.pred <- CJ(
  neg=neg.seq, pos=-neg.seq
)[, differentiable := TRUE][]
range(grid.dt[, neg-pos])

both.pred <- rbind(border.pred, grid.pred)
pred.tall <- melt(
  both.pred,
  id.vars=names(both.pred),
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda"
)[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[, {
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
}, keyby=list(pos, neg, differentiable)]
metrics.wide[auc==max(auc)] #max auc => aum>0.
metrics.wide[14:15, roc ]

one.pred.diff <- one.pred[["neg"]]-one.pred[["pos"]]
one.pred <- c(neg=3.5, pos=-3.5)
some.err[, example := label]
diff.dt <- aum::aum_diffs_penalty(some.err, names(one.pred))
ls.list <- aum::aum_line_search(diff.dt, pred.vec=one.pred, maxIterations = 10)

##compute slope and intercept of each of the 6 T_b(s) functions, plot
##them using geom_abline, and geom_point to represent the 9
##intersection points.
some.diff[, `:=`(
  slope=ifelse(label=="pos", 0, -1),
  intercept=pred.log.lambda-ifelse(label=="pos", 0, 6.5))]

metrics.tall <- melt(
  metrics.wide,
  measure.vars=c("aum", "auc"),
  variable.name="var.lower"
)[order(-differentiable)]
metrics.tall[, variable := toupper(var.lower)]
metrics.tall[, norm := (value-min(value))/(max(value)-min(value))]

ggplot()+
  geom_tile(aes(
    pos, neg, fill=norm),
    data=metrics.tall[differentiable==TRUE])+
  scale_fill_gradient(low="white", high="red")+
  facet_grid(. ~ variable)+
  coord_equal()

ls.points.tall <- melt(
  ls.list$line_search_result,
  id="step.size",
  measure=c("aum","auc"),
  variable.name="var.lower")
ls.points.tall[, variable := toupper(var.lower)]
metrics.tall[, pred.diff := neg-pos]
metrics.tall[, step.size := (one.pred.diff-pred.diff)/2]
extra <- 0.5
ls.segs.tall <- rbind(
  ls.list$line_search_result[, .(
    variable="AUC", 
    step.min=step.size, step.max=c(step.size[-1], step.size[.N]+extra), 
    value.min=auc.after, value.max=auc.after)],
  ls.list$line_search_result[, .(
    variable="AUM", 
    step.min=step.size, step.max=c(step.size[-1], step.size[.N]+extra), 
    value.min=aum, value.max=c(aum[-1], aum[.N]))])
ggplot()+
  ggtitle("Overview, select difference")+
  theme_bw()+
  theme(panel.margin=grid::unit(1, "lines"))+
  theme_animint(width=300, height=300)+
  facet_grid(variable ~ ., scales="free")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  geom_point(aes(
    step.size, value, fill=differentiable),
    size=3,
    shape=21,
    data=metrics.tall)+
  geom_point(aes(
    step.size, value),
    data=ls.points.tall,
    color="red")+
  geom_segment(aes(
    step.min, value.min,
    xend=step.max, yend=value.max),
    data=ls.segs.tall,
    color="red")+
  xlab("Prediction difference, f(neg) - f(pos)")+
  scale_y_continuous("", breaks=seq(0, 3, by=1))
