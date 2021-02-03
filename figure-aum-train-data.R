source("packages.R")

folds.dt <- fread("../feature-learning-benchmark/labeled_problems_folds.csv")
addMeta <- function(dt){
  dt[, set.name := sub("/.*", "", prob.dir)]
  dt[, problem := sub(".*/", "", prob.dir)]
  dt[folds.dt, on=list(set.name, problem)]
}
errors.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_errors.csv"))
possible.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_possible_errors.csv"))
features.dt <- fread("../feature-learning-benchmark/labeled_problems_features.csv")

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
  max.log.lambda=max(max.log.lambda)
), by=.(prob.dir, seg.i)]
test.fold.segs[, mid.log.lambda := (max.log.lambda+min.log.lambda)/2]
test.fold.targets <- penaltyLearning::targetIntervals(
  test.fold.errors, "prob.dir")
all.features.mat <- as.matrix(features.dt[, -1, with=FALSE])
test.fold.features.dt <- features.dt[test.fold.targets$prob.dir, on="prob.dir"]
all.features.mat <- as.matrix(test.fold.features.dt[, -1, with=FALSE])
sd.features <- apply(all.features.mat, 2, sd)
keep.feature <- is.finite(sd.features) & 0 < sd.features
finite.features.mat <- all.features.mat[, keep.feature]
test.fold.targets.mat <- test.fold.targets[, cbind(min.log.lambda, max.log.lambda)]
scaled.features.mat <- scale(finite.features.mat)
rownames(scaled.features.mat) <- test.fold.targets$prob.dir

## No need to set seed, unregularized learning algorithm is
## deterministic (resulting fit does not depend on random seed).
keep.obs <- apply(is.finite(test.fold.targets.mat), 1, any)
fit <- penaltyLearning::IntervalRegressionUnregularized(
  scaled.features.mat[keep.obs,], test.fold.targets.mat[keep.obs,])
computeAUM <- function(w, i){
  pred.pen.vec <- (scaled.features.mat %*% w) + i
  pred.dt <- data.table(
    prob.dir=rownames(pred.pen.vec),
    pred.log.lambda=as.numeric(pred.pen.vec))
  out.list <- penaltyLearning::ROChange(
    possible.errors, pred.dt, "prob.dir")
  out.list$intercept <- out.list$thresholds[
    threshold=="min.error", (min.thresh+max.thresh)/2]
  out.list$weight.vec <- w
  out.list
}
initial.roc <- computeAUM(coef(fit)[-1], 0)

this.roc <- initial.roc
neg.t.X.subtrain <- -t(scaled.features.mat)
iterations.dt.list <- list()
for(step.number in 1:20){
  print(iterations.dt.list[[paste(step.number)]] <- with(this.roc, data.table(
    step.number, aum, auc,
    min.errors=thresholds[threshold=="min.error", errors])))
  g.dt <- this.roc[["aum.grad"]]
  ## If aum.grad has some problems with no changes in error then
  ## they may be missing.
  g.vec <- rep(0, ncol(neg.t.X.subtrain))
  names(g.vec) <- colnames(neg.t.X.subtrain)
  g.vec[
    g.dt[["prob.dir"]]
  ] <- g.dt[["lo"]]
  direction.vec <- neg.t.X.subtrain %*% g.vec
  take.step <- function(s){
    this.roc$weight.vec + s*direction.vec
  }
  step.roc.list <- list()
  for(step.size in 10^seq(-5, 0, by=0.25)){
    step.roc.list[[paste(step.size)]] <- computeAUM(take.step(step.size), 0)
  }
  aum.vec <- sapply(step.roc.list, "[[", "aum")
  this.roc <- step.roc.list[[which.min(aum.vec)]]
}#iteration
iterations.dt <- do.call(rbind, iterations.dt.list)

roc.list <- list(initial=initial.roc, improved=this.roc)
roc.dt.list <- list()
aum.dt.list <- list()
for(pred.name in names(roc.list)){
  L <- roc.list[[pred.name]]
  roc.dt.list[[pred.name]] <- data.table(
    pred.name, L[["roc"]])
  aum.dt.list[[pred.name]] <- with(L, data.table(
    pred.name, aum, auc, thresholds[threshold=="min.error"]))
}

out.list <- list(
  iterations=iterations.dt,
  roc=do.call(rbind, roc.dt.list),
  auc=do.call(rbind, aum.dt.list))

saveRDS(out.list, "figure-aum-train-data.rds")
