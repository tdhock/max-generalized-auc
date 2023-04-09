library(ggplot2)
library(data.table)
library(batchtools)
library(aum)
library(dplyr)
cache.name <- "figure-line-grid-search-interactive-cache.rds"
if(FALSE){
  unlink(file.path(testFold.vec, cache.name), recursive=TRUE)
}

##TODO add for loop over line search set, subtrain or validation?

# args.dt[49:52] H3K4me3_TDH_ENCODE

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

# get the size of all folds of every dataset
# store each pair of fold name & size in a data table
# combine them into a giant data table
# order by size
datasets.by.size <- rbindlist(lapply(Sys.glob("../neuroblastoma-data/data/*/cv/*/testFolds/*"), function(x) {
  folder.size <- sum(file.info(list.files(x, all.files=TRUE, recursive=TRUE, full.names=TRUE))$size)
  data.table(size=folder.size, name=x)
}))[order(size),]
# get the data name and test fold from the path and store it in the table
datasets.by.size <- rbindlist(lapply(datasets.by.size$name, function(testFold.path) {
  size <- datasets.by.size[name==testFold.path]$size
  cv.path <- dirname(dirname(testFold.path))
  folds.csv <- file.path(cv.path, "folds.csv")
  cv.type <- basename(cv.path)
  test.fold <- basename(testFold.path)
  data.dir <- dirname(dirname(cv.path))
  data.name <- basename(data.dir)
  data.table(size, name=testFold.path, data.name, test.fold)
}))

(testFold.vec <- Sys.glob("../neuroblastoma-data/data/*/cv/*/testFolds/*"))
testFold.path <- "../neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/3"
seed <- 1
init.name="IntervalRegressionCV"
aum.type="count"
OneBatch <- function(testFold.path, aum.type, init.name){
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
  obj.sign.list <- list(aum=1)#list(auc=-1, aum=1)
  seeds <- c(3) #4
  for(seed in seeds) { #for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    for(algo in c("grid","exact","exactq","hybridA","hybridB"))for(objective in names(obj.sign.list)){
      start.time <- microbenchmark::get_nanotime()
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
      elapsed.time <- 0
      max.iterations <- if (algo == "exactq") {
        nrow(diffs.list$subtrain) * (nrow(diffs.list$subtrain)) - 1 / 2
      } else {
        nrow(diffs.list$subtrain)
      }
      
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
        current.time <- microbenchmark::get_nanotime() - start.time
        iteration.dt.list[[paste(
          seed, init.name, algo, step.number, objective
        )]] <- data.table(
          seed, init.name, algo, step.number, objective, elapsed.time, time=current.time, summary.dt)
        new.obj <- summary.dt.list$subtrain[[paste0(objective,".diffs")]]
        improvement <- obj.sign*(prev.obj-new.obj)
        cat(sprintf(
          "seed=%d init=%s algo=%s step=%d %s %f->%f\n",
          seed, init.name, algo, step.number, objective, prev.obj, new.obj))
        1e-5 < improvement
      }){
        ##while(step.number<2){
        pred.vec <- X.subtrain %*% weight.vec
        aum.list <- aum::aum(diffs.list$subtrain, pred.vec)
        pred.grad.vec <- rowMeans(aum.list$derivative_mat)
        direction.vec <- neg.t.X.subtrain %*% pred.grad.vec
        take.step <- function(s){
          weight.vec+s*direction.vec
        }
        
        ptm <- proc.time() # timer
        grid.result <- NULL
        exact.result <- NULL
        if (algo == "grid") {
          step.grid <- 10^seq(-9, 0)
          grid.dt <- data.table(step.size=step.grid)[, {
            step.weight <- take.step(step.size)
            grid.aum <- aum_auc(diffs.list$subtrain, X.subtrain %*% step.weight)
            with(grid.aum, data.table(auc, aum))
          }, by=step.size]
          grid.result <- grid.dt[, .(search="grid", step.size, auc, aum)]
        } else if (algo == "exact" || algo == "exactq") {
          LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec, maxIterations=max.iterations)
          exact.result <- LS$line_search_result[, .(search="exact", step.size, auc, aum)]
        } else if (algo == "hybridA") {
          LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec, maxIterations=max.iterations)
          exact.result <- LS$line_search_result[, .(search="exact", step.size, auc, aum)]
          search.result <- data.table(LS$line_search_result)
          search.result[, kink := .I/.N]
          best.row <- search.result[which.min(aum)]
          
          if (best.row$kink == 1) {
            # if kink == 1, we have chosen the very last step size we looked at.
            # run a grid search where we're at to find a larger step.size
            step.grid <- best.row$step.size * 10^seq(1, 4)
            grid.dt <- data.table(step.size=step.grid)[, {
              step.weight <- take.step(step.size)
              grid.aum <- aum_auc(diffs.list$subtrain, X.subtrain %*% step.weight)
              with(grid.aum, data.table(auc, aum))
            }, by=step.size]
            grid.result <- grid.dt[, .(search="grid", step.size, auc, aum)]
            best.grid.row <- grid.result[which.min(aum)]
            # if we got a better value using a grid search,
            # adjust our max iterations to search for more values next time
            if (best.grid.row$aum < best.row$aum) {
              max.iterations <- max.iterations * 2
            }
          }
        } else if (algo == "hybridB") {
          LS=aum::aum_line_search(diffs.list$subtrain, X.subtrain, weight.vec, maxIterations=max.iterations)
          exact.result <- LS$line_search_result[, .(search="exact", step.size, auc, aum)]
          search.result <- data.table(LS$line_search_result)
          search.result[, kink := .I/.N]
          best.row <- search.result[which.min(aum)]
          
          if (best.row$kink == 1) {
            # if kink == 1, we have chosen the very last step size we looked at.
            # run a grid search where we're at to find a larger step.size
            step.grid <- best.row$step.size * 2
            grid.dt <- data.table(step.size=step.grid)[, {
              step.weight <- take.step(step.size)
              grid.aum <- aum_auc(diffs.list$subtrain, X.subtrain %*% step.weight)
              with(grid.aum, data.table(auc, aum))
            }, by=step.size]
            grid.result <- grid.dt[, .(search="grid", step.size, auc, aum)]
            best.grid.row <- grid.result[which.min(aum)]
          }
        }
        elapsed.time <- (proc.time() - ptm)[["elapsed"]] # timer end
        
        steps.considered <- rbind(
          exact.result,
          grid.result
        )[, step.prop := seq(1, .N)/.N, by=search][]
        #considered.dt.list[[paste(
        #  seed, init.name, algo, objective, step.number
        #)]] <- data.table(
        #  seed, init.name, algo, objective, step.number, steps.considered)
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
      testFold.path, data.name, cv.type, test.fold))
  #steps=data.table(
  #  rbindlist(considered.dt.list),
  #  data.name, cv.type, test.fold))
}

args.dt <- data.table::CJ(
  testFold.path=datasets.by.size$name[1:109],#testFold.vec,
  aum.type=c("rate"),#,"count")
  init.name=c("zero", "IntervalRegressionCV")
)

## Run on SLURM.
registry.dir <- "figure-line-grid-search-interactive-registry"
registry.dir <- "figure-line-grid-search-interactive-registry-6"#4 datasets w hybrid A,B (1.2mb)
registry.dir <- "figure-line-grid-search-interactive-registry-7"#23 datasets (1.6mb)
registry.dir <- "figure-line-grid-search-interactive-registry-8"#70 datasets (20mb)
registry.dir <- "figure-line-grid-search-interactive-registry-9"#109 datasets (28mb)
registry.dir <- "figure-line-grid-search-interactive-registry-10"#109 datasets w/ init.name=c("zero", "IntervalRegressionCV")

if (FALSE) {
  reg=batchtools::loadRegistry(registry.dir, writeable = TRUE)
  #batchtools::clearRegistry(reg)
  batchtools::getStatus(reg=reg)
  batchtools::findExpired(reg=reg)
  status.dt <- batchtools::getJobStatus(reg=reg)
  status.dt[!is.na(error)]
  status.dt[!is.na(done)]
}

#analyze.
if(FALSE){
  done.ids <- status.dt[is.na(error), job.id]
  for(done.i in done.ids){
    job.id <- done.ids[[done.i]]
    args.row <- args.dt[job.id]
    ls.dir <- file.path(args.row$testFold.path, "line_search", "sets")
    dir.create(ls.dir, showWarnings = FALSE, recursive = TRUE)
    ls.csv <- file.path(ls.dir, paste0(args.row$aum.type, ".csv"))
    if(!file.exists(ls.csv)){
      cat(sprintf("%4d / %4d %s\n", done.i, length(done.ids), ls.csv))
      res <- batchtools::loadResult(job.id)
      best.steps <- res$steps[
        , .SD[which.min(aum)], by=.(
          seed,init.name,algo,objective,step.number
        )][,.(seed,init.name,algo,objective,step.number=step.number+1,search)]
      join.dt <- best.steps[res$sets, on=.(
        seed,init.name,algo,objective,step.number
      )]
      join.dt[is.na(search), table(step.number)]
      fwrite(join.dt, ls.csv)
    }  
  }
}

# args.dt[5:8]   CTCF_TDH_ENCODE 2407
# args.dt[13:16] H3K27ac_TDH_some 798
# args.dt[17:20] H3K27me3_RL_cancer 592
# args.dt[21:24] H3K27me3_TDH_some 296
# args.dt[29:32] H3K36me3_TDH_ENCODE 377
# args.dt[37:40] H3K36me3_TDH_other 554
# args.dt[c(13:16,21:24,29:32,37:40,17:20)]

if(FALSE){
  unlink(registry.dir, recursive=TRUE)
}
# CREATE REGISTRY AND RUN JOBS
reg <- batchtools::makeRegistry(file.dir=registry.dir)
#parallel.job.count <- 8
reg$cluster.functions <- makeClusterFunctionsMulticore()
# sample random test folds
#batchtools::batchMap(OneBatch, args=args.dt[sample(1:nrow(args.dt), 5)], reg=reg)
batchtools::batchMap(OneBatch, args=args.dt, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table,chunk=job.table$job.id)
#chunks <- data.frame(job.table,chunk=(job.table$job.id%%parallel.job.count)+1)
#chunks <- data.frame(job.table, chunk=1)
options(batchtools.verbose=TRUE, batchtools.progress=TRUE)
runJobs <- function() {
  ptm <- proc.time()
  batchtools::submitJobs(chunks, resources=list(
    #walltime = 50*5,#24 * 60 * 60,#seconds
    memory = 5000,#megabytes per cpu
    #max.concurrent.jobs = parallel.job.count,
    ncpus=1,  #>1 for multicore/parallel jobs.
    ntasks=1, #>1 for MPI jobs.
    chunks.as.arrayjobs=FALSE), reg=reg)
  total.time <- (proc.time() - ptm)[["elapsed"]]
  total.time
}
cat(sprintf("Finished running jobs in %f minutes\n", runJobs() / 60))

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
all.it.list <-
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
  meta.dt <- type.dt[1, .(
    data.name, cv.type, test.fold,
    gradient=sub(".csv","",basename(type.csv)))]
  type.csv <- type.csv.vec[[type.csv.i]]
  type.dt <- fread(type.csv)
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


# jadon tests

# load the first test and graph auc by each algorithm
if(FALSE) {
  result.one <- loadResult(1, reg)
  ggplot(result.one$sets[objective=="aum"]) +
    geom_line(aes(x=step.number, y=aum, color=algo, linetype=as.factor(set))) +
    facet_grid(init.name ~ algo, scales="free") +
    scale_x_log10() +
    scale_y_log10()
}

# load all results and build one big data table
result.sets.list <- list()
status.dt <- batchtools::getJobStatus(reg=reg)
completed.jobs <- status.dt[is.na(error)]$job.id
for (result.id in completed.jobs) {
  # ensure this job is done
  if (!is.na(status.dt[job.id==result.id]$done)) {
    r <- batchtools::loadResult(result.id, reg)
    r$sets[,result.id:=result.id]
    result.sets.list[[result.id]] <- r$sets
    r <- NULL
  }
}
result.sets <- do.call(rbind, result.sets.list)
result.sets.list <- NULL

algo.time.by.dataset <- data.table(result.sets[init.name=="zero"][objective=="aum"][set=="validation"] %>%
                                     group_by(result.id, algo, testFold.path) %>%
                                     reframe(total.time = sum(elapsed.time)))
colnames(algo.time.by.dataset)[colnames(algo.time.by.dataset) == "testFold.path"] <- "name"
results.with.dataset.size <- merge(algo.time.by.dataset, datasets.by.size, on=.(name))


# graph the complete results
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#DA72B2", "#D55E00", "#CC79A7")

# plot elapsed time per step of gradient descent for each algo/dataset
ggplot(result.sets[init.name=="zero"][objective=="aum"][set=="validation"]) +
  geom_point(aes(x=step.number, y=elapsed.time, color=algo, linetype=as.factor(set))) +
  facet_grid(init.name ~ result.id, scale="free") +
  #scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values=cbPalette)
ggsave("elapsed.time.png", width=1920*3, height=1080*3, units="px")

# plot aum for each dataset
ggplot(result.sets[init.name=="zero"][objective=="aum"][set=="validation"]) +
  geom_line(aes(x=step.number, y=aum, color=algo, linetype=as.factor(set))) +
  facet_wrap(vars(result.id), scale="free") +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values=cbPalette)
ggsave("aum.png", width=1920*3, height=1080*3, units="px")


ggplot(result.sets[init.name=="zero"][objective=="aum"][set=="validation"]) +
  geom_line(aes(x=step.number, y=auc, color=algo, linetype=as.factor(set))) +
  facet_wrap(vars(result.id), scale="free") +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values=cbPalette)
ggsave("auc.png", width=1920*3, height=1080*3, units="px")

# total time by algo (bar chart)
result.sets[init.name=="zero"][objective=="aum"][set=="validation"] %>%
  group_by(result.id, algo) %>%
  summarize(total.time = sum(elapsed.time)) %>%
  ggplot() +
  geom_col(aes(x=algo, y=total.time, fill=algo)) +
  #scale_y_log10() +
  geom_text(aes(x=algo, y=total.time, label=round(total.time,digits=1)), position=position_dodge(width=0.9), vjust=-0.25) +
  facet_grid(. ~ result.id) +
  scale_colour_manual(values=cbPalette)
ggsave("total.time.png", width=1920*3, height=1080*3, units="px")

ggplot(result.sets[init.name=="zero"][objective=="aum"][set=="validation"]) +
  geom_line(aes(x=time, y=aum, color=algo, linetype=as.factor(set))) +
  facet_wrap(vars(result.id), scale="free") +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values=cbPalette)
ggsave("aum.over.time.png", width=1920*3, height=1080*3, units="px")

# boxplot total time by algo
result.sets[init.name=="zero"][objective=="aum"][set=="validation"] %>%
  group_by(result.id, algo) %>%
  reframe(total.time = sum(elapsed.time)) %>%
  ggplot() +
  #geom_violin(aes(x=algo, y=total.time, fill=algo)) +
  geom_boxplot(aes(x=algo, y=total.time, fill=algo)) +
  scale_y_log10() +
  #geom_text(aes(x=algo, y=total.time, label=round(total.time,digits=3)), position=position_dodge(width=0.9), vjust=-0.25) +
  #facet_grid(. ~ result.id) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
ggsave("total.time.across.datasets.png", width=1920*3, height=1080*3, units="px")

results.with.dataset.size %>%
  ggplot() +
  geom_line(aes(x=size, y=total.time, color=algo)) +
  scale_y_log10() +
  scale_x_log10() +
  scale_colour_manual(values=cbPalette)
ggsave("size.affects.time.png", width=1920*2, height=1080*2, units="px")

results.with.dataset.size %>%
  #group_by(algo, s=signif(size, 3)) %>%
  #reframe(t=mean(total.time)) %>%
  ggplot(aes(x=size, y=total.time, color=algo, fill=algo)) +
  geom_point() +
  geom_smooth(level=0.70,span=0.6) +
  scale_y_log10() +
  scale_x_log10() +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  xlab("Size of Dataset (MB)") +
  ylab("Total time (seconds)") +
  ggtitle("something")
ggsave("size.affects.time.png", width=1920*3, height=1080*3, units="px")
