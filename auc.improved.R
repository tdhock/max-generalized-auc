source("packages.R")

folds.dt <- fread("../feature-learning-benchmark/labeled_problems_folds.csv")
addMeta <- function(dt){
  dt[, set.name := sub("/.*", "", prob.dir)]
  dt[, problem := sub(".*/", "", prob.dir)]
  dt[folds.dt, on=list(set.name, problem)]
}
errors.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_errors.csv"))
possible.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_possible_errors.csv"))

auc.improved.list <- list()

## compute derivative of Area under min(FP, FN).
fold.possible <- unique(folds.dt[, .(set.name, fold)])
i.possible <- 1:nrow(fold.possible)
N.possible <- paste(i.possible, "improved")
i.todo <- i.possible[!N.possible %in% names(auc.improved.list)]
biggest.step <- 0.1
for(i in seq_along(i.todo)){
  test.fold.i <- i.todo[[i]]
  cat(sprintf("%4d / %4d test folds TODO=%d\n", i, length(i.todo), test.fold.i))
  test.fold.info <- fold.possible[test.fold.i]
  test.fold.errors <- errors.dt[test.fold.info, on=.(set.name, fold)]
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
  ## initialization:
  pred.dt <- data.table(initial.pred)
  getROC <- function(p){
    L <- penaltyLearning::ROChange(possible.errors, p, "prob.dir")
    non.smooth <- L$aum.subdiff[lower != upper]
    if(nrow(non.smooth))print(non.smooth)
    L
  }
  step.number <- 1
  step.size <- 1
  roc.list <- getROC(pred.dt)
  improvement <- Inf
  while(1e-6 < improvement){
    ## these depend on predictions:
    while({
      step.dt <- pred.dt[roc.list$aum.subdiff, .(
        prob.dir,
        pred.log.lambda = pred.log.lambda-step.size*lower
      ), on=.(prob.dir)]
      step.list <- getROC(step.dt)
      roc.list$aum < step.list$aum
    }){
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
  pred.list <- list(
    initial=initial.pred,
    ##mid=mid.pred,
    improved=pred.dt)
  for(pred.name in names(pred.list)){
    pred <- pred.list[[pred.name]]
    L <- penaltyLearning::ROChange(possible.errors, pred, "prob.dir")
    print(L$auc)
    auc.improved.list[[paste(test.fold.i, pred.name)]] <- with(L, data.table(
      test.fold.i,
      test.fold.info,
      pred.name,
      roc=list(list(roc)),
      thresholds[threshold=="min.error"],
      auc))
  }
}
(auc.improved <- do.call(rbind, auc.improved.list))

saveRDS(auc.improved, "auc.improved.rds")
