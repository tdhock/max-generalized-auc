source("packages.R")

## need to clone https://github.com/tdhock/feature-learning-benchmark
folds.dt <- fread("../feature-learning-benchmark/labeled_problems_folds.csv")
addMeta <- function(dt){
  dt[, set.name := sub("/.*", "", prob.dir)]
  dt[, problem := sub(".*/", "", prob.dir)]
  dt[folds.dt, on=list(set.name, problem)]
}
errors.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_errors.csv"))
possible.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_possible_errors.csv"))

diff.info <- errors.dt[, aum::aum_diffs_penalty(data.table(example=prob.dir, min.lambda=exp(min.log.penalty), fp, fn), unique(prob.dir)), by=.(set.name, fold)]
diff.counts <- diff.info[, .(
  breakpoints=.N,
  examples=length(unique(example))
), by=.(set.name, fold)]
diff.tall <- melt(diff.counts, measure.vars=c("breakpoints", "examples"))
diff.tall[, .(
  max=max(value),
  mean=mean(value),
  min=min(value),
  folds=.N
), by=variable]

fold.counts <- possible.dt[, list(
  examples=.N,
  labels=sum(labels)
), by=.(set.name, fold)]
fold.counts[, list(
  folds=.N,
  min.examples=min(examples),
  max.examples=max(examples)
)]

test.fold.info <- folds.dt[set.name=="H3K4me3_XJ_immune" & fold==4]
test.fold.errors <- errors.dt[test.fold.info, on=.(set.name, fold, problem)]
test.fold.errors[, min.log.lambda := min.log.penalty]
test.fold.errors[, max.log.lambda := max.log.penalty]
test.fold.errors[, seg.i := cumsum(
  c(1, diff(fp)!=0 | diff(fn) != 0)), by=.(prob.dir)]
possible.errors <- possible.dt[test.fold.errors, on=list(
  set.name, fold, prob.dir)]
possible.errors[, possible.fn := possible.tp]
test.fold.segs <- test.fold.errors[, .(
  min.log.lambda=min(min.log.lambda),
  max.log.lambda=max(max.log.lambda),
  fp=fp[1],
  fn=fn[1]
), by=.(prob.dir, seg.i)]
prob.dir.ord <- unique(test.fold.segs$prob.dir)
diff.fp.fn <- aum::aum_diffs_penalty(
  test.fold.segs[, `:=`(example=prob.dir, min.lambda = exp(min.log.lambda))],
  prob.dir.ord)
pred.vec <- initial.pred[prob.dir.ord, pred.log.lambda, on="prob.dir"]
aum.list <- aum::aum(diff.fp.fn, pred.vec)
descent.direction.vec <- -rowMeans(aum.list$derivative_mat)
direction.dt <- data.table(
  example=seq(0, length(pred.vec)-1),
  weight=pred.vec,
  direction=descent.direction.vec)
(diffs.with.direction <- direction.dt[diff.fp.fn, on="example"][, `:=`(
  slope = -direction,
  intercept=pred-weight
)][])
##TODO Jadon plug in your code.
ggplot()+
  geom_point(aes(
    0, intercept),
    data=diffs.with.direction)+
  geom_abline(aes(
    slope=slope, intercept=intercept),
    data=diffs.with.direction)+
  coord_cartesian(xlim=c(-5,5))

test.fold.segs[, mid.log.lambda := (max.log.lambda+min.log.lambda)/2]
test.fold.targets <- penaltyLearning::targetIntervals(
  test.fold.errors, "prob.dir")
test.fold.targets[, width := max.log.lambda-min.log.lambda]
initial.pred <- test.fold.targets[order(width==Inf, -width), data.table(
  prob.dir,
  pred.log.lambda=ifelse(
    max.log.lambda==Inf, min.log.lambda+1, ifelse(
      min.log.lambda==-Inf, max.log.lambda-1,
      (min.log.lambda+max.log.lambda)/2)
  )
)]
initial.pred[!is.finite(pred.log.lambda), pred.log.lambda := 0]

test.fold.breaks <- test.fold.errors[, .(breaks=.N-1), by=prob.dir]
test.fold.breaks[, .(
  total.breaks=sum(breaks),
  max.breaks=max(breaks),
  mean.breaks=mean(breaks),
  min.breaks=min(breaks),
  examples=.N
)]
diff.counts[test.fold.info[1], on=c("set.name", "fold")]

## initialization:
pred.dt <- data.table(initial.pred)
getROC <- function(p){
  L <- penaltyLearning::ROChange(possible.errors, p, "prob.dir")
  non.smooth <- L$aum.grad[lo != hi]
  if(nrow(non.smooth))print(non.smooth)
  L
}
step.number <- 1
step.size <- 1
roc.list <- getROC(pred.dt)
iterations.dt.list <- list()
improvement <- Inf
while(1e-6 < improvement){
  ## these depend on predictions:
  while({
    step.dt <- pred.dt[roc.list$aum.grad, .(
      prob.dir,
      pred.log.lambda = pred.log.lambda-step.size*lo
    ), on=.(prob.dir)]
    step.list <- getROC(step.dt)
    roc.list$aum < step.list$aum
  }){
    ## TODO Jadon replace step size halving with your algorithm,
    ## either go all the way to the end of the path, quadratic
    ## max_iterations=N*(N-1)/2 or just stop with log linear
    ## algorithm, max_iterations=N.
    step.size <- step.size/2
  }
  cat(sprintf(
    "step=%d size=%e aum=%f->%f auc=%f->%f\n",
    step.number,
    step.size,
    roc.list$aum,
    step.list$aum,
    roc.list$auc,
    step.list$auc))
  iterations.dt.list[[paste(step.number)]] <- data.table(
    step.number,
    aum=roc.list$aum,
    auc=roc.list$auc,
    min.errors=roc.list$thresholds[threshold=="min.error", errors])
  improvement <- roc.list$aum-step.list$aum
  pred.dt <- step.dt
  roc.list <- step.list
  step.number <- step.number + 1
  step.size <- step.size*2
}
mid.pred <- test.fold.segs[pred.dt, .(
  prob.dir,
  improved.pred=pred.log.lambda,
  mid.log.lambda), on=.(
    prob.dir,
    min.log.lambda < pred.log.lambda,
    max.log.lambda > pred.log.lambda)]
mid.pred[, pred.log.lambda := ifelse(
  is.finite(mid.log.lambda), improved.pred, mid.log.lambda)]

iterations.dt <- do.call(rbind, iterations.dt.list)

pred.list <- list(
  initial=initial.pred,
  ##mid=mid.pred,
  improved=pred.dt)
out.auc.list <- list()
out.roc.list <- list()
for(pred.name in names(pred.list)){
  pred <- pred.list[[pred.name]]
  L <- penaltyLearning::ROChange(possible.errors, pred, "prob.dir")
  print(L$auc)
  out.auc.list[[paste(pred.name)]] <- with(L, data.table(
    pred.name,
    thresholds[threshold=="min.error"],
    auc, aum))
  out.roc.list[[paste(pred.name)]] <- data.table(pred.name, L$roc)
}

out.list <- list(
  iterations=iterations.dt,
  roc=do.call(rbind, out.roc.list),
  auc=do.call(rbind, out.auc.list))

saveRDS(out.list, "figure-aum-optimized-data.rds")
