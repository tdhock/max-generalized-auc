library(data.table)
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
testFold.vec <- Sys.glob("../neuroblastoma-data/data/*/cv/*/testFolds/*")
testFold.path <- "../neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/4"
OneFold <- function(testFold.path){
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
      denominator="count")
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
  for(seed in 1:4)for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    for(algo in c("grid","exact","hybrid")){
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
      prev.aum <- Inf
      new.aum <- -Inf
      step.number <- 0
      ##while(new.aum < prev.aum){
      while(step.number<2){
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
        prev.aum <- summary.dt.list$subtrain$aum
        summary.dt <- do.call(rbind, summary.dt.list)
        iteration.dt.list[[paste(seed, init.name, algo, step.number)]] <- data.table(
          seed, init.name, algo, step.number, summary.dt)
        LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec)
        pred.vec <- X.subtrain %*% weight.vec
        aum.list <- aum::aum(diffs.list$subtrain, pred.vec)
        pred.grad.vec <- rowMeans(aum.list$derivative_mat)
        direction.vec <- neg.t.X.subtrain %*% pred.grad.vec
        step.grid <- 10^seq(-9, 0)
        take.step <- function(s){
          weight.vec+s*direction.vec
        }
        grid.dt <- data.table(step.size=step.grid, aum=sapply(step.grid, function(step){
          step.weight <- take.step(step)
          aum::aum(diffs.list$subtrain, X.subtrain %*% step.weight)$aum
        }))
        steps.considered <- rbind(
          if(algo!="grid")LS$line_search_result[, .(search="exact", step.size, aum)],
          if(algo!="exact")grid.dt[, .(search="grid", step.size, aum)]
        )[, step.prop := seq(1, .N)/.N, by=search][]
        considered.dt.list[[paste(seed, init.name, algo, step.number)]] <- data.table(
          seed, init.name, algo, step.number, steps.considered)
        best.step <- steps.considered[which.min(aum)]
        new.weight.vec <- take.step(best.step$step.size)
        after.roc <- computeROC(new.weight.vec, 0, "subtrain")
        new.aum <- after.roc$aum
        cat(sprintf(
          "seed=%d init.name=%s algo=%s step=%d aum %f->%f\n",
          seed, init.name, algo, step.number, prev.aum, new.aum))
        intercept <- after.roc$thresholds[
          threshold=="min.error", (max.thresh+min.thresh)/2]
        weight.vec <- new.weight.vec
      }#step.number
    }#algo
  }#seed/init.name
  list(
    sets=data.table(
      do.call(rbind, iteration.dt.list),
      data.name, cv.type, test.fold),
    steps=data.table(
      rbindlist(considered.dt.list),
      data.name, cv.type, test.fold))
}

cache.name <- "figure-line-grid-search-interactive-cache.rds"
if(FALSE){
  unlink(file.path(testFold.vec, cache.name), recursive=TRUE)
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

