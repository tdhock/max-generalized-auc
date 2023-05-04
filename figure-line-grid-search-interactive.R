library(ggplot2)
library(data.table)
  obj.sign.list <- list(auc=-1, aum=1)
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
testFold.path <- "../neuroblastoma-data/data/detailed/cv/R-3.6.0-chrom/testFolds/1"               
seed <- 1
init.name="zero"
aum.type="count"
OneBatch <- function(testFold.path, aum.type){
  library(data.table)
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
  data.list$evaluation[, `:=`(
    min.fn=min(fn),
    max.fp=max(fp),
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
    fn=fn-min.fn,
    possible.fp=max.fp
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
      data.list$evaluation,
      seqs.set,
      denominator=aum.type)
    diffs.list[[s]] <- seqs.diff
  }
  totals <- colSums(diffs.list$subtrain[, .(fp_diff, fn_diff)])
  X.subtrain <- X.finite[set.vec=="subtrain",]
  neg.t.X.subtrain <- -t(X.subtrain)
  seqs.train <- with(seqs.list, c(subtrain, validation))
  y.train <- data.list[["outputs"]][
    seqs.train,
    cbind(min.log.lambda, max.log.lambda),
    on="sequenceID"]
  keep <- apply(is.finite(y.train), 1, any)
  X.train <- X.finite[seqs.train, ]
  N.param <- ncol(X.finite)+1
  init.param <- structure(
    rep(0, N.param),
    names=c("(Intercept)",colnames(X.finite)))
  init.fun.list <- list(
    IntervalRegressionCV=function(){
      fit <- penaltyLearning::IntervalRegressionCV(
        X.train[keep, ],
        y.train[keep, ])
      some.param <- fit[["param.mat"]]
      init.param[names(some.param)] <- some.param
      init.param
    },
    zero=function(){
      init.param+rnorm(N.param)
    }
  )
  iteration.dt.list <- list()
  considered.dt.list <- list()
  for(seed in 1:4)for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    ##for(algo in c("grid","exact","hybrid"))
    for(algo in c("grid","hybrid"))
      ##for(objective in names(obj.sign.list)){
      for(objective in "aum"){
      computeROC <- function(w, i, set){
        pred.pen.vec <- (X.finite %*% w) + i
        pred.dt <- data.table(
          sequenceID=rownames(pred.pen.vec),
          pred.log.lambda=-as.numeric(pred.pen.vec))
        is.set <- set.vec==set
        set.dt <- pred.dt[is.set]
        L <- penaltyLearning::ROChange(
          data.list$evaluation, set.dt, "sequenceID")
        alist <- aum_auc(diffs.list[[set]], pred.pen.vec[ seqs.list[[set]], ])
        L$aum.diffs <- alist$aum
        L$auc.diffs <- alist$auc
        L
      }
      aum_auc <- function(diffs.dt, pred.vec){
        aum.list <- aum::aum(diffs.dt, pred.vec)
        before.dt <- data.table(aum.list$total_error, key="thresh")[, `:=`(
          TPR_before=1-fn_before/-totals[["fn_diff"]],
          FPR_before=fp_before/totals[["fp_diff"]])]
        aum.list$auc <- before.dt[, .(
          FPR=c(FPR_before, 1),
          TPR=c(TPR_before, 1)
        )][, sum((FPR[-1]-FPR[-.N])*(TPR[-1]+TPR[-.N])/2)]
        aum.list
      }
      obj.sign <- obj.sign.list[[objective]]
      weight.vec <- int.weights[-1]
      intercept <- int.weights[1]
      prev.obj <- Inf*obj.sign
      step.number <- 0
      while({
        summary.dt.list <- list()
        for(set in names(seqs.list)){
          set.PL <- computeROC(weight.vec, intercept, set)
          summary.dt.list[[set]] <- with(set.PL, data.table(
            set,
            thresholds[threshold=="predicted"],
            auc, aum, auc.diffs, aum.diffs))
        }
        summary.dt <- do.call(rbind, summary.dt.list)
        iteration.dt.list[[paste(
          seed, init.name, algo, step.number, objective
        )]] <- data.table(
          seed, init.name, algo, step.number, objective, summary.dt)
        new.obj <- summary.dt.list$subtrain[[paste0(objective,".diffs")]]
        improvement <- obj.sign*(prev.obj-new.obj)
        cat(sprintf(
          "seed=%d init=%s algo=%s step=%d %s %f->%f\n",
          seed, init.name, algo, step.number, objective, prev.obj, new.obj))
        1e-5 < improvement
      }){
        ##while(step.number<2){
        LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec)
        pred.vec <- X.subtrain %*% weight.vec
        aum.list <- aum::aum(diffs.list$subtrain, pred.vec)
        pred.grad.vec <- rowMeans(aum.list$derivative_mat)
        direction.vec <- neg.t.X.subtrain %*% pred.grad.vec
        step.grid <- 10^seq(-9, 0)
        take.step <- function(s){
          weight.vec+s*direction.vec
        }
        grid.dt <- data.table(step.size=step.grid)[, {
          step.weight <- take.step(step.size)
          grid.aum <- aum_auc(diffs.list$subtrain, X.subtrain %*% step.weight)
          with(grid.aum, data.table(auc, aum))
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
        weight.vec <- take.step(best.step$step.size)
        new.aum <- aum::aum(diffs.list$subtrain, X.subtrain %*% weight.vec)
        err.thresh <- data.table(
          new.aum$total_error,key="thresh"
        )[, err_before := fp_before+fn_before][, .(
          thresh=c(thresh[1]-1,thresh[-1]-diff(thresh)/2,thresh[.N]+1),
          err=c(err_before,sum(diffs.list$subtrain$fp_diff))
        )]
        intercept <- err.thresh[which.min(err), thresh]
        step.number <- step.number+1
        prev.obj <- new.obj
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

args.dt <- data.table::CJ(
  testFold.path=testFold.vec,
  aum.type=c("rate","count")
)

## Run on SLURM.
registry.dir <- "figure-line-grid-search-interactive-registry"
registry.dir <- "figure-line-grid-search-interactive-registry-fixed"
reg=batchtools::loadRegistry(registry.dir)
batchtools::getStatus(reg=reg)
batchtools::findExpired(reg=reg)
status.dt <- batchtools::getJobStatus(reg=reg)
status.dt[!is.na(error)]
status.dt[!is.na(done)]

#analyze.
job.id <- 9#join.dt[init.name=="zero" & seed==1 & objective=="aum" & search != "exact"]
done.ids <- status.dt[is.na(error), job.id]
for(done.i in seq_along(done.ids)){
  job.id <- done.ids[[done.i]]
  args.row <- args.dt[job.id]
  ls.dir <- file.path(args.row$testFold.path, "line_search", "sets")
  dir.create(ls.dir, showWarnings = FALSE, recursive = TRUE)
  ls.csv <- file.path(ls.dir, paste0(args.row$aum.type, ".csv"))
  if(!file.exists(ls.csv)){
    cat(sprintf("%4d / %4d %s\n", done.i, length(done.ids), ls.csv))
    res <- batchtools::loadResult(job.id)
    best.steps <- res$steps[
    , .SD[which.min(obj.sign.list[[objective]]*get(objective))], by=.(
      seed,init.name,algo,objective,step.number
    )][,.(seed,init.name,algo,objective,step.number=step.number+1,search)]
    join.dt <- best.steps[res$sets, on=.(
      seed,init.name,algo,objective,step.number
    )]
    join.dt[is.na(search), table(step.number)]
    fwrite(join.dt, ls.csv)
  }  
}

if(FALSE){
  unlink(registry.dir, recursive=TRUE)
}
reg <- batchtools::makeRegistry(registry.dir)
batchtools::batchMap(OneBatch, args=args.dt, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 24*60*60,#seconds
  memory = 32000,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=TRUE), reg=reg)

batchtools::getStatus(reg=reg)
status.dt <- batchtools::getJobStatus(reg=reg)
status.dt[!is.na(error)]
status.dt[!is.na(done)]

batchtools::testJob(4, reg=reg)
args.dt[21]

## seed=1 init.name=IntervalRegressionCV algo=exact step=33146 auc 0.955147->0.955147
##  *** caught bus error ***
## address 0x153dc1917f40, cause 'non-existent physical address'
## Traceback:
##  1: aum_sort_interface(error.diff.df, pred.vec)


##job.id=376 Error in penaltyLearning::ROChange(data.list$evaluation, data.table(sequenceID = seqs.set,  : \n  no positive labels => fix by excluding data?

##job.id=354 Error in X.finite %*% w : non-conformable arguments => fixed by processing output of IRCV, which returns w that may be smaller than number of features.

##job.id=4,12 Error in aum_sort_interface(error.diff.df, pred.vec) : \n  fp should be non-negative => fixed by updating check in aum C++ code.

##job.id=21 Error in while (obj.sign * (new.obj - prev.obj) < 1e-06) { : \n  missing value where TRUE/FALSE needed => fixed by using aum package which always gives finite aum.

## Run locally.
for(args.i in 1:nrow(args.dt)){
  args.row <- args.dt[args.i]
  cache.rds <- args.row[, file.path(testFold.path, paste0(aum.type, ".rds"))]
  all.it.list[[args.i]] <- if(file.exists(cache.rds)){
    readRDS(cache.rds)
  }else{
    cat(sprintf("%4d / %4d\n", args.i, length(args.dt)))
    print(args.row)
    iteration.list <- do.call(OneBatch, args.row)
    saveRDS(iteration.list, cache.rds)
  }
}

## analyze.
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

#analyze 2
type.csv.vec <- Sys.glob(file.path(testFold.vec, "line_search","sets", "*.csv"))
selected.dt.list <- list()
for(type.csv.i in seq_along(type.csv.vec)){
  type.csv <- type.csv.vec[[type.csv.i]]
  type.dt <- fread(type.csv)
  meta.dt <- type.dt[1, .(
    data.name, cv.type, test.fold,
    gradient=sub(".csv","",basename(type.csv)))]
  ## does max auc get better auc than min aum?
  valid.dt <- type.dt[
    set=="validation"
  ][, step.prop := step.number/max(step.number), by=.(seed,init.name,algo,objective)]
  compare.obj.dt <- valid.dt[
  , .SD[which.max(auc), .(step.number,step.prop,valid.auc=auc)], by=.(seed,init.name,algo,objective)]
  not.zero <- valid.dt[0 < step.number]
  search.counts <- dcast(
    compare.obj.dt[not.zero, on=.(seed,init.name,algo,objective,step.number>=step.number),nomatch=0L],
    seed+init.name+algo+objective~search,
    length)
  selected.dt.list[[type.csv]] <- data.table(
    meta.dt, search.counts[compare.obj.dt, on=.(seed,init.name,algo,objective)])
}
selected.dt <- rbindlist(selected.dt.list)
fwrite(selected.dt, "figure-line-grid-search-interactive-selected.csv")

ggplot()+
  facet_grid(init.name + objective ~ ., labeller=label_both, scales="free")+
  geom_point(aes(
    auc, algo),
    data=compare.obj.dt)+
  scale_x_continuous(
    "Best validation AUC")


selected.dt <- data.table::fread("figure-line-grid-search-interactive-selected.csv")
exact.used <- selected.dt[algo=="hybrid" & 0<exact]
exact.used[, table(objective, gradient)]
grid.only <- selected.dt[algo=="grid"]
exact.used[grid.only, grid.auc = i.valid.auc, on=.(
  data.name, cv.type, test.fold, gradient, seed, init.name, objective
)]
exact.used[, hist(valid.auc-grid.auc)]

ggplot()+
  geom_point(aes(
    valid.auc, grid.auc),
    shape=1,
    data=exact.used)+
  coord_equal()

selected.wide <- dcast(
  selected.dt[step.number>0],
  data.name + cv.type + test.fold + gradient + seed + init.name + objective ~ algo,
  value.var=c('step.number','valid.auc','exact')
)[!is.na(exact_hybrid)][, hybrid_used_exact := ifelse(exact_hybrid==0, "no", "yes")]
ggplot()+
  facet_grid(. ~ hybrid_used_exact)+
  geom_point(aes(
    valid.auc_grid, valid.auc_hybrid),
    data=selected.wide)
##Shouldn't grid=hybrid when exact=0?
one <- selected.wide[exact_hybrid==0 & valid.auc_hybrid != valid.auc_grid & objective=="aum"][1]
selected.dt[one, on=.(
  data.name, cv.type, test.fold, gradient, seed, init.name, objective
)]
##type.csv.i=1 [seed==1 & init.name=="zero" & objective=="auc"]
##type.csv="../neuroblastoma-data/data/ATAC_JV_adipose/cv/equal_labels/testFolds/1/line_search/sets/count.csv"
best.algos <- selected.dt[
, .SD[which.max(valid.auc)],
  by=.(data.name, cv.type, test.fold, gradient, seed, init.name, objective)]
best.algos[, .(count=.N), by=algo]
best.algos[, .(count=.N), keyby=.(algo, objective)]
best.algos[objective=="auc", .(count=.N), keyby=.(algo, gradient)]
