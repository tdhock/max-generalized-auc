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
    set.name, fold, prob.dir)][, possible.fn := possible.tp]
  possible.segs <- possible.errors[, .(
    min.log.lambda=min(min.log.lambda),
    max.log.lambda=max(max.log.lambda)
  ), by=.(
    prob.dir, seg.i, fp, fn, errors, possible.fp, possible.fn, labels
  )][, `:=`(
    min.lambda = exp(min.log.lambda),
    example=prob.dir
  )]
  ## Check for non-zero at end of err fun.
  possible.segs[min.log.lambda == -Inf & fn > 0]
  possible.segs[min.log.lambda == Inf & fp > 0]
  test.fold.targets <- penaltyLearning::targetIntervals(
    possible.segs, "prob.dir")
  prob.ord <- test.fold.targets$prob.dir
  aum.diffs <- aum::aum_diffs_penalty(possible.segs, prob.ord)
  min.err.pred.dt <- test.fold.targets[, data.table(
    prob.dir,
    pred.log.lambda=fcase(
      min.log.lambda>-Inf & max.log.lambda==Inf, min.log.lambda+1, 
      min.log.lambda==-Inf & max.log.lambda<Inf, max.log.lambda-1,
      min.log.lambda>-Inf & max.log.lambda<Inf, (min.log.lambda+max.log.lambda)/2,
      min.log.lambda==-Inf & max.log.lambda==Inf, 0)
  )]
  getROC <- function(p){
    L <- penaltyLearning::ROChange(possible.segs, p, "prob.dir")
    non.smooth <- L$aum.grad[lo != hi]
    if(nrow(non.smooth))print(non.smooth)
    L
  }
  getAUM <- function(pred.vec){
    L <- aum::aum(aum.diffs, pred.vec)
    L$grad <- with(L, ifelse(
      derivative_mat[,1] == derivative_mat[,2],
      (derivative_mat[,1]+derivative_mat[,2])/2,
      0))
    L
  }
  possible.segs[prob.dir==prob.ord[1], .(fp,fn,min.log.lambda)]
  aum.diffs[example==0]
  init.list <- list(
    min.error=-min.err.pred.dt$pred.log.lambda,
    zero=rep(0, nrow(min.err.pred.dt)))
  for(initialization in names(init.list)){
    current.pred <- initial.pred <- init.list[[initialization]]
    step.number <- 1
    step.size <- 1
    ##roc.list <- getROC(pred.dt)
    aum.list <- getAUM(current.pred)
    ##data.table(pred.dt, current.pred)[, pred.log.lambda-current.pred]
    improvement <- Inf
    while(1e-3 < improvement){
      ## these depend on predictions:
      while({
        ## step.dt <- pred.dt[roc.list$aum.grad, .(
        ##   prob.dir,
        ##   pred.log.lambda = pred.log.lambda-step.size*ifelse(
        ##     sign(lo)==sign(hi), (lo+hi)/2, 0)
        ## ), on=.(prob.dir)]
        ## step.list <- getROC(step.dt)
        step.pred <- current.pred - step.size*aum.list$grad
        step.list <- getAUM(step.pred)
        aum.list$aum < step.list$aum
      }){
        step.size <- step.size/2
      }
      improvement <- aum.list$aum-step.list$aum
      cat(sprintf(
        "step=%d size=%e aum=%f->%f diff=%f\n",
        step.number,
        step.size,
        aum.list$aum,
        step.list$aum,
        improvement))
      ## pred.dt <- step.dt
      current.pred <- step.pred
      aum.list <- step.list
      step.number <- step.number + 1
      step.size <- step.size*2
    }
    pred.list <- list(
      initial=initial.pred,
      improved=current.pred)
    for(pred.name in names(pred.list)){
      pred <- data.table(
        prob.dir=prob.ord, 
        pred.log.lambda=-pred.list[[pred.name]])
      L <- penaltyLearning::ROChange(possible.segs, pred, "prob.dir")
      print(L$aum)
      auc.improved.list[[paste(test.fold.i, improvement, pred.name)]] <- 
        with(L, data.table(
          test.fold.i,
          test.fold.info,
          initialization,
          pred.name,
          roc=list(roc),
          thresholds[threshold=="min.error"],
          auc, aum))
    }
  }
}
(auc.improved <- do.call(rbind, auc.improved.list))

saveRDS(auc.improved, "auc.improved.rds")
