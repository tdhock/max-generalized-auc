library(ggplot2)
library(data.table)
cache.name <- "figure-line-grid-search-interactive-cache.rds"
if(FALSE){
  unlink(file.path(testFold.vec, cache.name), recursive=TRUE)
}

##TODO add for loop over line search set, subtrain or validation?

## > mb[per.set, on=list(set)][order(labels)]
##     megabytes                      set labels
##  1:       554       H3K36me3_TDH_other    200
##  2:       377      H3K36me3_TDH_ENCODE    338
##  3:       375       H3K4me3_TDH_ENCODE    525
##  4:       592       H3K27me3_RL_cancer    570
##  5:       798         H3K27ac_TDH_some    627
##  6:       906      H3K36me3_TDH_immune    630
##  7:       296        H3K27me3_TDH_some    696
##  8:      2407          CTCF_TDH_ENCODE   1378
##  9:      3223           H3K4me1_TDH_BP   1584
## 10:      5871       H3K36me3_AM_immune   1743
## 11:      6407          ATAC_JV_adipose   3241
## 12:      3017       H3K4me3_PGP_immune   3780
## 13:      2902       H3K4me3_TDH_immune   3807
## 14:      5421 H3K27ac-H3K4me3_TDHAM_BP  15961
(testFold.vec <- Sys.glob("../neuroblastoma-data/data/*/cv/*/testFolds/*"))
testFold.path <- "../neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/4"
aum.type="rate"
OneBatch <- function(testFold.path, aum.type){
  cv.path <- dirname(dirname(testFold.path))
  folds.csv <- file.path(cv.path, "folds.csv")
  cv.type <- basename(cv.path)
  test.fold <- basename(testFold.path)
  data.dir <- dirname(dirname(cv.path))
  data.name <- basename(data.dir)
  data.list <- list()
  for(f in c("inputs", "outputs", "evaluation")){
    f.csv.xz <- file.path(data.dir, paste0(f, ".csv.xz"))
    if(file.exists(f.csv.xz)){
      system(paste("unxz", f.csv.xz))
    }
    f.csv <- file.path(data.dir, paste0(f, ".csv"))
    f.dt <- data.table::fread(f.csv)
    data.list[[f]] <- f.dt
  }
  ## replace positive fn at end with 0 to avoid AUM=Inf.
  data.list$evaluation[, min.fn := min(fn), by=sequenceID]
  data.list$evaluation[, `:=`(
    min.lambda = exp(min.log.lambda),
    example=sequenceID
  ), by=sequenceID]
  bad <- data.list$evaluation[min.log.lambda == -Inf & min.fn < fn]
  if(nrow(bad)){
    print(bad)
  }
  data.list$evaluation[min.log.lambda == -Inf & 0 < fn]
  ## code below not necessary since this does not happen in our real
  ## data sets, but it could theoretically in some data.
  data.list$aum.input <- data.table(data.list$evaluation)[, `:=`(
    possible.fn=possible.fn-min.fn,
    fn=fn-min.fn
  ), by=sequenceID]
  ## read folds.  
  folds.dt <- data.table::fread(folds.csv)
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set := rep(
    c("subtrain", "validation"), l=.N)]
  folds.dt[, table(fold, set)]
  X.all <- scale(data.list$inputs[, -1])#rm seqID.
  rownames(X.all) <- data.list$inputs$sequenceID
  X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
  set.vec <- folds.dt[rownames(X.finite), set, on="sequenceID"]
  seqs.list <- list()
  diffs.list <- list()
  aum.vec.list <- list()
  for(s in unique(folds.dt$set)){
    seqs.set <- folds.dt[s==set, sequenceID]
    seqs.list[[s]] <- seqs.set
    seqs.diff <- aum::aum_diffs_penalty(
      data.list$aum.input,
      seqs.set,
      denominator=aum.type)
    diffs.list[[s]] <- seqs.diff
    zero <- rep(0, length(seqs.set))
    aum.vec.list[[s]] <- rbind(
      "aum::aum"=aum::aum(seqs.diff, zero)$aum,
      "penaltyLearning::ROChange"=penaltyLearning::ROChange(
        data.list$evaluation,
        data.table(sequenceID=seqs.set, pred.log.lambda=zero),
        problem.vars="sequenceID")$aum
    )
  }
  print(do.call(data.frame, aum.vec.list))
  X.subtrain <- X.finite[set.vec=="subtrain",]
  neg.t.X.subtrain <- -t(X.subtrain)
  seqs.train <- with(seqs.list, c(subtrain, validation))
  y.train <- data.list[["outputs"]][
    seqs.train,
    cbind(min.log.lambda, max.log.lambda),
    on="sequenceID"]
  keep <- apply(is.finite(y.train), 1, any)
  X.train <- X.finite[seqs.train, ]
  init.fun.list <- list(
    IntervalRegressionCV=function(){
      fit <- penaltyLearning::IntervalRegressionCV(
        X.train[keep, ],
        y.train[keep, ])  
      fit[["param.mat"]]
    },
    zero=function(){
      N.param <- ncol(X.finite)+1
      rep(0, N.param)+rnorm(N.param)
    }
  )
  iteration.dt.list <- list()
  considered.dt.list <- list()
  obj.sign.list <- list(auc=-1, aum=1)
  for(seed in 1:4)for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    for(algo in c("grid","exact","hybrid"))for(objective in names(obj.sign.list)){
      obj.sign <- obj.sign.list[[objective]]
      weight.vec <- int.weights[-1]
      intercept <- int.weights[1]
      computeROC <- function(w, i, set){
        pred.pen.vec <- (X.finite %*% w) + i
        pred.dt <- data.table(
          sequenceID=rownames(pred.pen.vec),
          pred.log.lambda=-as.numeric(pred.pen.vec))
        is.set <- set.vec==set
        set.dt <- pred.dt[is.set]
        penaltyLearning::ROChange(
          data.list$evaluation, set.dt, "sequenceID")
      }
      prev.obj <- Inf*obj.sign
      new.obj <- -Inf*obj.sign
      step.number <- 0
      while(obj.sign*(new.obj-prev.obj) < 1e-6){
      ##while(step.number<2){
        step.number <- step.number+1
        summary.dt.list <- list()
        for(set in names(seqs.list)){
          set.PL <- computeROC(weight.vec, intercept, set)
          summary.dt.list[[set]] <- with(set.PL, data.table(
            set,
            thresholds[threshold=="predicted"],
            auc,
            aum))
        }
        prev.obj <- summary.dt.list$subtrain[[objective]]
        summary.dt <- do.call(rbind, summary.dt.list)
        iteration.dt.list[[paste(
          seed, init.name, algo, step.number, objective
        )]] <- data.table(
          seed, init.name, algo, step.number, objective, summary.dt)
        LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec)
        pred.vec <- X.subtrain %*% weight.vec
        aum.list <- aum::aum(diffs.list$subtrain, pred.vec)
        pred.grad.vec <- rowMeans(aum.list$derivative_mat)
        direction.vec <- neg.t.X.subtrain %*% pred.grad.vec
        step.grid <- 10^seq(-9, 0)
        take.step <- function(s){
          weight.vec+s*direction.vec
        }
        totals <- colSums(diffs.list$subtrain[, .(fp_diff, fn_diff)])
        grid.dt <- data.table(step.size=step.grid)[, {
          step.weight <- take.step(step.size)
          grid.aum <- aum::aum(diffs.list$subtrain, X.subtrain %*% step.weight)
          before.dt <- data.table(grid.aum$total_error, key="thresh")[, `:=`(
            TPR_before=1-fn_before/-totals[["fn_diff"]],
            FPR_before=fp_before/totals[["fp_diff"]])]
          auc <- before.dt[, .(
            FPR=c(FPR_before, 1),
            TPR=c(TPR_before, 1)
          )][, sum((FPR[-1]-FPR[-.N])*(TPR[-1]+TPR[-.N])/2)]
          data.table(auc, aum=grid.aum$aum)
        }, by=step.size]
        steps.considered <- rbind(
          if(algo!="grid")LS$line_search_result[, .(search="exact", step.size, auc, aum)],
          if(algo!="exact")grid.dt[, .(search="grid", step.size, auc, aum)]
        )[, step.prop := seq(1, .N)/.N, by=search][]
        considered.dt.list[[paste(
          seed, init.name, algo, objective, step.number
        )]] <- data.table(
          seed, init.name, algo, objective, step.number, steps.considered)
        best.step <- steps.considered[which.min(obj.sign*get(objective))]
        new.weight.vec <- take.step(best.step$step.size)
        after.roc <- computeROC(new.weight.vec, 0, "subtrain")
        new.obj <- after.roc[[objective]]
        cat(sprintf(
          "seed=%d init.name=%s algo=%s step=%d %s %f->%f\n",
          seed, init.name, algo, step.number, objective, prev.obj, new.obj))
        intercept <- after.roc$thresholds[
          threshold=="min.error", (max.thresh+min.thresh)/2]
        weight.vec <- new.weight.vec
      }#step.number
    }#algo/objective
  }#seed/init.name
  list(
    sets=data.table(
      do.call(rbind, iteration.dt.list),
      data.name, cv.type, test.fold),
    steps=data.table(
      rbindlist(considered.dt.list),
      data.name, cv.type, test.fold))
}

all.it.list <- list()
for(testFold.i in seq_along(testFold.vec)){
  fdir <- testFold.vec[testFold.i]
  cache.rds <- file.path(fdir, cache.name)
  all.it.list[[testFold.i]] <- if(file.exists(cache.rds)){
    readRDS(cache.rds)
  }else{
    cat(sprintf("%4d / %4d %s\n", testFold.i, length(testFold.vec), fdir))
    iteration.list <- OneFold(fdir)
    saveRDS(iteration.list, cache.rds)
    iteration.list
  }
}

cache.vec <- Sys.glob(file.path(
  "../neuroblastoma-data/data/*/cv/*/testFolds/*",
  cache.name))
for(cache.i in seq_along(cache.vec)){
  cache.rds <- cache.vec[[cache.i]]
  L <- readRDS(cache.rds)
  algo.cols <- c("seed","init.name","algo")
  step.cols <- c(algo.cols,"step.number")
  best.steps <- L$steps[, .SD[which.min(aum)], by=step.cols][,c(step.cols,"search"),with=FALSE]
  join.dt <- L$sets[set != "test"][best.steps, on=step.cols]
  min.dt <- join.dt[set=="validation", .SD[which.min(aum)], by=.(seed,init.name,set)]

  ggplot()+
    geom_line(aes(
      step.number, aum, color=algo),
      data=join.dt)+
    geom_point(aes(
      step.number, aum, color=algo),
      shape=1,
      data=join.dt[search=="exact"])+
    geom_point(aes(
      step.number, aum, color=algo),
      data=min.dt)+
    facet_wrap(~seed+ init.name + set,scales="free")
  
}