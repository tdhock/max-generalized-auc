library(animint2)
library(data.table)

data(neuroblastomaProcessed, package="penaltyLearning")
data(neuroblastoma, package="neuroblastoma")
e <- function(label, profile.id, chromosome){
  data.table(label, profile.id=factor(profile.id), chromosome=factor(chromosome))
}
select.dt <- rbind(
  e("positive", 4, 2),
  e("negative", 513, 3))
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
border.pred <- with(dlist, positive[ #orange dots
  negative,
  data.table(
    differentiable=FALSE,
    positive=pred.log.lambda,
    negative=i.pred.log.lambda),
  on="id",
  allow.cartesian=TRUE])
grid.pred <- data.table( #black dots
  differentiable=TRUE,
  positive=0,
  negative=seq(dmin, dmax, by=0.05))
both.pred <- rbind(border.pred, grid.pred)
both.pred[, pred.diff := negative-positive]
pred.tall <- melt(
  both.pred,
  measure.vars=select.dt$label,
  variable.name="label",
  value.name="pred.log.lambda")[select.dt, nomatch=0L, on="label"]
metrics.wide <- pred.tall[order(pred.diff)][, {
  L <- penaltyLearning::ROChange(some.err, .SD, "label")
  positive <- pred.log.lambda[label=="positive"]
  with(L, data.table(
    aum, auc,
    SM=roc[min.thresh < max.thresh, sum(min.fp.fn)],
    roc=list(roc[, `:=`(
      min.thresh=min.thresh+positive,
      max.thresh=max.thresh+positive
    )])
  ))
}, by=list(pred.diff, differentiable)]
metrics.wide[auc==max(auc)] #max auc => aum>0.
metrics.wide[14:15, roc ]
initial.pred.diff <- 6.5
metrics.wide[,step.size := initial.pred.diff - pred.diff]

#jadon
library(dplyr)
library(Rcpp)
aumLineSearchLib <- sourceCpp("RaumLineSearch.cpp")

##compute slope and intercept of each of the 6 T_b(s) functions, plot
##them using geom_abline, and geom_point to represent the 9
##intersection points.
some.diff[, `:=`(
  prediction.value.change=ifelse(label=="positive", 0, -1), # Delta-i s
  intercept=pred.log.lambda-ifelse(label=="positive", 0, initial.pred.diff))]
# pred.log.lambda == v_b
some.diff[, slope := -prediction.value.change]
line.data <- data.table(some.diff)

ggplot() +
  geom_point(aes(step.size, aum, color = differentiable), data = metrics.wide)

intersect.fn <- function(a, b){
  x <- (b[1] - a[1]) / (a[2] - b[2])
  y <- a[1] + a[2] * x
  return(xy=c(x, y))
}

# build the intersection points
intersection.col.names <- c(
  "x", "y", 
  "intercept.a", "slope.a", 
  "intercept.b", "slope.b",
  "fp.diff.a", "fn.diff.a",
  "fp.diff.b", "fn.diff.b"
)
intersection.points <- data.frame(matrix(ncol = length(intersection.col.names), nrow = 0))
colnames(intersection.points) <- intersection.col.names
for (a in 1:nrow(some.diff)) {
  for (b in a:nrow(some.diff)) {
    if (a != b) {
      line.a <- some.diff[a]
      line.b <- some.diff[b]
      intercept.a <- line.a$intercept
      intercept.b <- some.diff[b]$intercept
      slope.a <- some.diff[a]$slope
      slope.b <- some.diff[b]$slope
      found.intersection <- intersect.fn(
        c(intercept.a, slope.a),
        c(intercept.b, slope.b)
      )
      if (is.finite(found.intersection)) {
        inter <- c(
          x=found.intersection[1],
          y=found.intersection[2],
          intercept.a=intercept.a,
          slope.a=slope.a,
          intercept.b=intercept.b,
          slope.b=slope.b,
          fp.diff.a=line.a$fp.diff,
          fn.diff.a=line.a$fn.diff,
          fp.diff.b=line.b$fp.diff,
          fn.diff.b=line.b$fn.diff
        )
        intersection.points[nrow(intersection.points) + 1,] <- inter
      }
    }
  }
}
# sort intersections by x value
intersection.points <- distinct(intersection.points[order(intersection.points$x),])
intersection.points$id <- seq_along(intersection.points[,1])

aum.points <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(aum.points) <- c("x", "y")
aum.slope <- 0


for (i in 2:nrow(some.diff)) {
  slope.diff <- some.diff[i]$slope - some.diff[i-1]$slope
  aum.slope <- aum.slope + slope.diff * some.diff[i]$prediction.value.change
}

# `dataframe$column == DoubleValue` doesn't work?
doubleEquals <- function (table, column, value) {
  table[column >= (value - 0.00005) & column <= (value + 0.00005)]
}

fp.vec <- list()
fn.vec <- list()
min.vec <- list()
for (i in 2:nrow(some.diff)) {
  fp.delta.delta <- some.diff[i]$fp.diff - some.diff[i-1]$fp.diff
  fn.delta.delta <- some.diff[i]$fn.diff - some.diff[i-1]$fn.diff
  fp.vec <- append(fp.vec, fp.delta.delta)
  fn.vec <- append(fn.vec, fn.delta.delta)
  min.vec <- append(min.vec, min(fp.delta.delta, fn.delta.delta))
}

for (i in 1:nrow(intersection.points)) {
  point <- intersection.points[i,]
  line.a <- doubleEquals(line.data, line.data$intercept, point$intercept.a)
  line.b <- doubleEquals(line.data, line.data$intercept, point$intercept.b)
  print(line.b)
  print(paste("A fp.diff: ", line.a$fp.diff))
  print(paste("A fn.diff: ", line.a$fn.diff))
  print(paste("B fp.diff: ", line.b$fp.diff))
  print(paste("B fn.diff: ", line.b$fn.diff))
}

initial.aum <- metrics.wide[step.size == 0]$aum
aum.points[1,] <- c(0, initial.aum)
for (i in 1:nrow(intersection.points)) {
  point <- intersection.points[i,]
  new.aum <- ((point$x - aum.points[i,]$x) * aum.slope) + aum.points[i,]$y
  aum.points[nrow(aum.points) + 1,] <- c(point$x, new.aum)
  # the aum.slope should be updated here
  aum.slope <- 0
  
  # fake the slope update by using the next non-differentiable point
  # that lines up with the next intersection point
  if (i < nrow(intersection.points && i > 1)) {
    previous.point
    next.point <- intersection.points[i+1,]
    closest.point <- metrics.wide[step.size >= (next.point$x - 0.00005) & step.size <= (next.point$x + 0.00005)]
    next.aum <- closest.point$aum
    rise <- next.aum - new.aum
    run <- next.point$x - point$x
    aum.slope <- rise / run
  }
}

aum.segments <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(aum.segments) <- c("x", "y", "xend", "yend")
for (i in 1:nrow(aum.points)-1) {
  aum.segments[i,] <- c(
    x=aum.points[i,]$x,
    y=aum.points[i,]$y, 
    xend=aum.points[i+1,]$x, 
    yend=aum.points[i+1,]$y
  )
}

print(intersection.points)

add.kind <- function(df, kind){
  data.frame(df, kind=factor(kind, c("aum", "inter")))
}

aum.plot <- ggplot() +
  ggtitle("AUM Intersection Points") +
  # light grey verticle lines showing the intersection x values
  geom_vline(data = add.kind(intersection.points, "aum"), aes(xintercept = x), color="grey") +
  geom_vline(data = add.kind(intersection.points, "inter"), aes(xintercept = x), color="grey") +
  # normal step size vs aum graph
  geom_point(data = add.kind(metrics.wide, "aum"), aes(step.size, aum, color = differentiable)) +
  # threshold lines
  geom_abline(data = add.kind(some.diff, "inter"), aes(intercept = intercept, slope = slope)) +
  # intersection points on the threshold lines
  geom_point(data = add.kind(intersection.points, "inter"), aes(x = x, y = y), color="orange") +
  # computed aum line segments
  geom_segment(data = add.kind(aum.segments, "aum"), aes(x = x, y = y, xend = xend, yend = yend)) +
  facet_grid(kind ~ ., scales="free")
aum.plot

intersection.viz <- animint(
  title="AUM and intersections",
  aum=aum.plot
)
animint2gist(intersection.viz)

grid.arrange(aum.plot, intersections.plot)
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
  positive=0,
  negative=pred.diff,
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
)[, label := "negative"]
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
  overview=ggplot()+
    ggtitle("Overview, select difference")+
    theme(panel.margin=grid::unit(1, "lines"))+
    theme(text=element_text(size = 15))+
    theme_animint(width=300, height=300)+
    facet_grid(variable ~ ., scales="free")+
    scale_fill_manual(values=c(
      "TRUE"="black",
      "FALSE"="orange"))+
    geom_point(aes(
      pred.diff, value, fill=differentiable),
      size=4,
      shape=21,
      data=metrics.tall)+
    make_tallrect(metrics.tall, "pred.diff")+ 
    xlab("Prediction difference, f(negative) - f(positive)")+
    coord_cartesian(xlim=c(dmin,dmax))+
    scale_y_continuous("", breaks=seq(0, 3, by=1)),
  data=ggplot()+
    ggtitle("Data, labels, predicted changepoint models")+
    theme_bw()+
    theme(legend.position="none")+
    theme_animint(width=600, height=300)+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6, fill=annotation),
      alpha=0.5,
      data=nb.some$annotations)+
    scale_fill_manual(
      "label",
      values=c(
        breakpoint="violet",
        normal="orange"))+
    geom_point(aes(
      position/1e6, logratio),
      color="grey50",
      data=nb.some$profiles)+
    geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6,
      color=error.type),
      data=selected.err,
      showSelected=c("pred.diff", "roc.point"),
      size=5,
      fill="transparent")+
    geom_segment(aes(
      start.pos/1e6, mean,
      xend=end.pos/1e6, yend=mean),
      data=selected.segs,
      size=3,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    geom_vline(aes(
      xintercept=start.pos/1e6),
      data=selected.segs[start>1],
      size=2,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    scale_color_manual(leg,values=err.colors)+
    facet_grid(label ~ ., labeller=label_both)+
    scale_y_continuous(
      "DNA copy number (microarray logratio)")+
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
      data=pred.tall.thresh,
      showSelected=c("pred.diff", "roc.point"))+
    geom_segment(aes(
      min.log.lambda, value,
      xend=max.log.lambda, yend=value,
      color=error.type, size=error.type),
      showSelected="error.type",
      data=some.err.tall)+
    geom_segment(aes(
      positive, -Inf,
      xend=negative, yend=-Inf),
      data=pred.tall.thresh.wide,
      showSelected=c("pred.diff", "roc.point"))+
    geom_text(aes(
      negative-0.1, -0.3,
      label=sprintf("pred.diff=%.2f", pred.diff)),
      hjust=1,
      data=pred.tall.thresh.wide,
      showSelected=c("pred.diff", "roc.point"))+
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
      color=NA,
      alpha=0.5,
      showSelected="pred.diff",
      show.roc.dt)+
    geom_segment(aes(
      min.thresh, value,
      xend=max.thresh, yend=value,
      color=error.type, size=error.type),
      showSelected="pred.diff",
      data=show.roc.tall)+
    geom_vline(aes(
      xintercept=text.constant),
      showSelected=c("pred.diff", "roc.point"),
      color=text.color,
      alpha=0.5,
      data=show.roc.dt)+
    geom_text(aes(
      text.constant, -0.25, label=roc.point),
      showSelected="pred.diff",
      size=text.size,
      color=text.color,
      data=show.roc.dt)+
    geom_text(aes(
      -1.5, 0.25, label=sprintf("AUM=%.2f", aum)),
      data=metrics.wide,
      showSelected="pred.diff")+
    geom_tallrect(aes(
      xmin=min.thresh, xmax=max.thresh),
      data=show.roc.dt,
      fill=text.color,
      clickSelects="roc.point",
      showSelected="pred.diff",
      color="transparent",
      alpha=0.1)+
    scale_y_continuous(
      "Label errors",
      breaks=c(0,1))+
    scale_color_manual(leg,values=err.colors)+
    scale_size_manual(leg,values=err.sizes)+
    geom_blank(aes(
      x, y),
      data=data.table(x=0, y=c(-0.4,1.4)))+
    scale_x_continuous(
      "Constant added to predicted values"),
  roc=ggplot()+
    ggtitle("ROC curve, select point")+
    theme_bw()+
    theme(panel.grid.minor=element_blank())+
    theme_animint(width=300, height=300)+
    geom_path(aes(
      FPR, TPR),
      showSelected="pred.diff",
      data=show.roc.dt)+
    geom_text(aes(
      0.5, 0.5, label=paste0("AUC=", auc)),
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
      data=show.roc.dt,
      size=4,
      alpha=0.5,
      color=text.color,
      showSelected=c("pred.diff", "roc.point"))+
    geom_text(aes(
      text.FPR, TPR+0.01, label=roc.point),
      size=text.size,
      color=text.color,
      showSelected="pred.diff",
      clickSelects="roc.point",
      data=show.roc.dt),
  time=list(
    variable="pred.diff",
    ms=500)
)
animint2gist(viz)
