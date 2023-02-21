library(animint2)
library(data.table)

data(neuroblastomaProcessed, package="penaltyLearning")
data(neuroblastoma, package="neuroblastoma")
ex <- function(label, profile.id, chromosome){
  data.table(
    label, 
    profile.id=factor(profile.id), 
    chromosome=factor(chromosome))
}
select.dt <- rbind(
  ex("pos", 4, 2),
  ex("neg", 513, 3))
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
dmax <- 3.05
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
    pos=pred.log.lambda,
    neg=i.pred.log.lambda),
  on="id",
  allow.cartesian=TRUE]
)[, diff := neg-pos]
neg.seq <- seq(dmin, dmax, by=0.025)
grid.pred <- CJ(
  pos=-neg.seq, 
  neg=neg.seq
)[, diff := neg-pos]
range(grid.pred[, neg-pos])

grid.uniq.diff <- grid.pred[, .(
  pos=0,
  neg=unique(diff)
)][, diff := neg-pos]
both.pred <- rbind(
  data.table(differentiable=FALSE, border.pred), 
  data.table(differentiable=TRUE, grid.uniq.diff)
)
pred.tall <- melt(
  both.pred,
  id.vars=c("diff","differentiable"),
  measure.vars=c("neg","pos"),
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
}, keyby=list(diff, differentiable)]
myjoin <- function(d.val, dt.pred){
  metrics.wide[differentiable==d.val][dt.pred, on="diff"]
}
both.roc <- rbind(
  myjoin(TRUE, grid.pred),
  myjoin(FALSE, border.pred))

pred.list <- list(
  nine.cross=c(neg=3.05, pos=-2.95),
  flat=c(neg=2.75, pos=-3),
  after.two=c(neg=2.5, pos=-2.75),
  increasing=c(neg=2.5, pos=-2.5))
diff.grid.list <- list()
ls.points.list <- list()
ls.segs.list <- list()
abline.dt.list <- list()
vline.dt.list <- list()
pred.points.list <- list()
heat.step.list <- list()
for(pred.name in names(pred.list)){
  one.pred <- pred.list[[pred.name]]
  one.pred.diff <- one.pred[["neg"]]-one.pred[["pos"]]
  some.err[, example := label]
  diff.dt <- aum::aum_diffs_penalty(some.err, names(one.pred))
  ls.list <- aum::aum_line_search(diff.dt, pred.vec=one.pred, maxIterations = 10)
  ##compute slope and intercept of each of the 6 T_b(s) functions, plot
  ##them using geom_abline, and geom_point to represent the 9
  ##intersection points.
  some.diff[, `:=`(
    slope=ifelse(label=="pos", 0, -1),
    intercept=pred.log.lambda-ifelse(label=="pos", 0, 6.5))]
  denom <- sum(ls.list$gradient_pred*c(1,-1))
  ToStep <- function(d){
    ifelse(d==0, 0, if(denom==0)NA else d/denom)
  }
  metrics.tall <- melt(
    both.roc,
    measure.vars=c("aum", "auc"),
    variable.name="var.lower"
  )[order(-differentiable)][, pred.diff := neg-pos][
  , step.size := ToStep(one.pred.diff-pred.diff)
  ][!is.na(step.size)]
  metrics.tall[, variable := toupper(var.lower)]
  metrics.tall[
  , norm := (value-min(value))/(max(value)-min(value)), by=variable]
  max.step <- ls.list$line_search_result[, (3*step.size[.N]-step.size[.N-1])/2]
  heat.step.list[[pred.name]] <- 
    data.table(t(sapply(c(
      ls.list$line_search_result$step.size, max.step
    ), function(s){
      one.pred-ls.list$gradient*s
    })), pred.name)
  if(length(max.step)==0)max.step <- Inf
  ls.segs.list[[pred.name]] <- data.table(
    pred.name, rbind(
      ls.list$line_search_result[, .(
        variable="AUC", 
        step.min=step.size, step.max=c(step.size[-1], max.step), 
        value.min=auc.after, value.max=auc.after)],
      ls.list$line_search_result[, .(
        variable="AUM", 
        step.min=step.size, step.max=c(step.size[-1], max.step), 
        value.min=aum, value.max=c(
          aum[-1], 
          if(max.step==Inf)aum else (max.step-step.size[.N])*aum.slope.after[.N]+aum[.N]))]))
  ls.points.list[[pred.name]] <- melt(
    data.table(pred.name, ls.list$line_search_result),
    measure=c("aum","auc"),
    variable.name="var.lower"
  )[, variable := toupper(var.lower)]
  diff.grid.list[[pred.name]] <- unique(metrics.tall[, .(
    pred.name, pred.diff, step.size, differentiable, variable, value
  )])
  abline.dt.list[[pred.name]] <- data.table(
    pred.name, variable="threshold", search="exact", ls.list$line_search_input)
  vline.dt.list[[pred.name]] <- unique(metrics.tall[
    differentiable==FALSE & variable=="AUC", 
    .(pred.name, step.size)])
  pred.points.list[[pred.name]] <- data.table(
    pred.name, t(one.pred))
}
heat.step <- rbindlist(heat.step.list)
pred.points <- rbindlist(pred.points.list)
diff.grid <- rbindlist(diff.grid.list)
ls.points <- rbindlist(ls.points.list)
ls.segs <- rbindlist(ls.segs.list)
abline.dt <- rbindlist(abline.dt.list)
vline.dt <- rbindlist(vline.dt.list)
ggplot()+
  geom_tile(aes(
    pos, neg, fill=norm),
    data=metrics.tall[differentiable==TRUE])+
  geom_point(aes(
    pos, neg),
    color="red",
    data=heat.step[, .SD[-.N], by=pred.name])+
  geom_line(aes(
    pos, neg, group=pred.name),
    color="red",
    data=heat.step)+
  geom_point(aes(
    pos, neg),
    shape=21,
    fill="white",
    data=pred.points)+
  scale_fill_gradient(low="white", high="blue")+
  facet_grid(. ~ variable)+
  coord_equal()


ggplot()+
  ggtitle("Overview, select step size")+
  geom_vline(aes(
    xintercept=step.size),
    data=vline.dt,
    color="grey")+
  theme_bw()+
  theme(panel.margin=grid::unit(1, "lines"))+
  theme_animint(width=300, height=300)+
  geom_abline(aes(
    slope=slope, intercept=intercept, color=search),
    data=abline.dt)+
  geom_point(aes(
    0, intercept, color=search),
    data=abline.dt)+
  facet_grid(variable ~ pred.name, scales="free_y")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  scale_color_manual(values=c(
    exact="red",
    grid="black"))+
  geom_point(aes(
    step.size, value, fill=differentiable, color=search),
    size=3,
    shape=21,
    data=data.table(search="grid", diff.grid))+
  geom_point(aes(
    step.size, value, color=search),
    data=data.table(search="exact", ls.points))+
  geom_segment(aes(
    step.min, value.min,
    color=search,
    xend=step.max, yend=value.max),
    data=data.table(search="exact", ls.segs))+
  xlab("Step size")+
  ##scale_y_continuous("", breaks=seq(0, 3, by=1))
  scale_y_continuous("")

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(1, "lines"))+
  theme_animint(width=300, height=300)+
  facet_grid(variable ~ ., scales="free")+
  scale_fill_manual(values=c(
    "TRUE"="black",
    "FALSE"="orange"))+
  scale_color_manual(values=c(
    exact="red",
    grid="black"))+
  geom_point(aes(
    pred.diff, value, fill=differentiable, color=search),
    size=3,
    shape=21,
    data=data.table(search="grid", diff.grid))+
  scale_y_continuous("", breaks=seq(0, 3, by=1))

min.max.step <- diff.grid[, seq(min(step.size), max(step.size), l=51)]
slope.int.lines <- abline.dt[, data.table(
  step.size=min.max.step,
  threshold=intercept+slope*min.max.step
), 
by=.(pred.name, variable, search, intercept, slope)
][min(abline.dt$intercept) < threshold & threshold < max(abline.dt$intercept)]
animint(
  out.dir="figure-line-search-interactive",
  heat=ggplot()+
    ggtitle("Loss function, select predictions")+
    theme_bw()+
    theme_animint(width=600, height=400)+
    geom_tile(aes(
      pos, neg, fill=norm),
      data=metrics.tall[differentiable==TRUE])+
    geom_point(aes(
      pos, neg),
      showSelected="pred.name",
      color="red",
      data=heat.step[, .SD[-.N], by=pred.name])+
    geom_line(aes(
      pos, neg, group=pred.name),
      color="red",
      showSelected="pred.name",
      data=heat.step)+
    geom_point(aes(
      pos, neg),
      shape=21,
      fill="white",
      size=4,
      clickSelects="pred.name",
      data=pred.points)+
    scale_fill_gradient(low="white", high="blue")+
    facet_grid(. ~ variable)+
    coord_equal(),
  step=ggplot()+
    ggtitle("Line search for selected predictions")+
    geom_vline(aes(
      xintercept=step.size),
      showSelected="pred.name",
      data=vline.dt,
      color="grey")+
    theme_bw()+
    theme(panel.margin=grid::unit(1, "lines"))+
    theme_animint(width=400, height=400)+
    geom_line(aes(
      step.size, threshold, color=search, group=intercept),
      showSelected="pred.name",
      data=slope.int.lines)+
    facet_grid(variable ~ ., scales="free")+
    scale_fill_manual(values=c(
      "TRUE"="black",
      "FALSE"="orange"))+
    scale_color_manual(values=c(
      exact="red",
      grid="black"))+
    geom_point(aes(
      step.size, value, fill=differentiable, color=search),
      showSelected="pred.name",
      size=3,
      shape=21,
      data=data.table(search="grid", diff.grid))+
    geom_point(aes(
      step.size, value, color=search),
      showSelected="pred.name",
      data=data.table(search="exact", ls.points))+
    geom_segment(aes(
      step.min, value.min,
      color=search,
      xend=step.max, yend=value.max),
      showSelected="pred.name",
      data=data.table(search="exact", ls.segs))+
    xlab("Step size")+
    scale_y_continuous("")
)
