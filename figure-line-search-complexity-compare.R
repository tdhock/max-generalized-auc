library(ggplot2)
library(data.table)
setDTthreads(1)
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
testFold.path <- "../neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/3"
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
  n.subtrain.diffs <- nrow(diffs.list$subtrain)
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
      init.param[rownames(some.param)] <- some.param
      init.param
    },
    near.zero=function(){
      init.param+rnorm(N.param)
    }
  )
  loss.dt.list <- list()
  time.dt.list <- list()
  iterations.dt.list <- list()
  Breaks <- nrow(diffs.list$subtrain)
  maxIterations.list <- list(
    grid=NULL,
    min.aum="min.aum",
    exactL=Breaks,
    exactQ=Breaks*(Breaks-1)/2)
  for(seed in 1:4)for(init.name in c("near.zero","IntervalRegressionCV")){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    for(maxIterations.name in names(maxIterations.list)){
      maxIterations <- maxIterations.list[[maxIterations.name]]
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
      obj.sign <- 1
      objective <- "aum"
      weight.vec <- int.weights[-1]
      intercept <- int.weights[1]
      prev.obj <- Inf*obj.sign
      step.number <- 0
      elapsed.seconds <- system.time({
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
          loss.dt.list[[paste(
            seed, init.name, step.number, maxIterations.name
          )]] <- data.table(
            seed, init.name, step.number, maxIterations.name, summary.dt)
          new.obj <- summary.dt.list$subtrain[[paste0(objective,".diffs")]]
          improvement <- obj.sign*(prev.obj-new.obj)
          cat(sprintf(
            "seed=%d init=%s step=%d %s %f->%f\n",
            seed, init.name, step.number, maxIterations.name, prev.obj, new.obj))
          ##step.number < 2 &&
          1e-5 < improvement
        }){
          pred.vec <- X.subtrain %*% weight.vec
          aum.list <- aum::aum(diffs.list$subtrain, pred.vec)
          pred.grad.vec <- rowMeans(aum.list$derivative_mat)
          direction.vec <- neg.t.X.subtrain %*% pred.grad.vec
          take.step <- function(s){
            weight.vec+s*direction.vec
          }
          candidate.dt <- if(identical(maxIterations.name, "grid")){
            step.grid <- 10^seq(-9, 0)
            iterations <- NA_integer_
            data.table(step.size=step.grid)[, {
              step.weight <- take.step(step.size)
              grid.aum <- aum_auc(diffs.list$subtrain, X.subtrain %*% step.weight)
              with(grid.aum, data.table(auc, aum))
            }, by=step.size]
          }else{
            LS=aum::aum_line_search(
              diffs.list$subtrain, X.subtrain, weight.vec, maxIterations=maxIterations)
            iterations <- LS[["line_search_result"]][
            , if(identical(maxIterations,"min.aum"))q.size else .N]
            LS$line_search_result
          }
          best.step.size <- candidate.dt[which.min(aum), step.size]
          weight.vec <- take.step(best.step.size)
          new.aum <- aum::aum(diffs.list$subtrain, X.subtrain %*% weight.vec)
          err.thresh <- data.table(
            new.aum$total_error,key="thresh"
          )[, err_before := fp_before+fn_before][, .(
            thresh=c(thresh[1]-1,thresh[-1]-diff(thresh)/2,thresh[.N]+1),
            err=c(err_before,sum(diffs.list$subtrain$fp_diff))
          )]
          intercept <- err.thresh[which.min(err), thresh]
          step.number <- step.number+1
          iterations.dt.list[[paste(
            seed, init.name, step.number, maxIterations.name
          )]] <- data.table(
            seed, init.name, step.number, maxIterations.name, iterations)
          prev.obj <- new.obj
        }#step.number
      })[["elapsed"]]
      time.dt.list[[paste(
        seed, init.name, maxIterations.name
      )]] <- data.table(
        seed, init.name, maxIterations.name, elapsed.seconds)
    }#maxIterations
  }#seed/init.name
  with_meta <- function(L){
    data.table(data.name, cv.type, test.fold, rbindlist(L))
  }
  list(
    sets=with_meta(loss.dt.list),
    iterations=with_meta(iterations.dt.list),
    time=with_meta(time.dt.list))
}
args.dt <- data.table::CJ(
  testFold.path=testFold.vec,
  aum.type=c("rate","count")
)

## Run on SLURM.
registry.dir <- "figure-line-search-complexity-compare-grid"
if(FALSE){
  unlink(registry.dir, recursive=TRUE)
}
reg <- batchtools::makeRegistry(registry.dir)
batchtools::batchMap(OneBatch, args=args.dt, reg=reg)
batchtools::testJob(4, reg=reg)
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

args.dt[21]

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
type.csv.vec <- Sys.glob(file.path(testFold.vec, "line_search","complexity", "*.csv"))
total.dt.list <- list()
for(type.csv.i in seq_along(type.csv.vec)){
  type.csv <- type.csv.vec[[type.csv.i]]
  aum.type <- sub(".csv","",basename(type.csv))
  type.dt <- fread(type.csv)
  type.total.dt <- type.dt[, .(
    aum.type, steps=.N, sum.iterations=sum(q.size), mean.iterations=mean(q.size)
  ), by=.(
    data.name, cv.type, test.fold, seed, init.name, maxIterations, n.subtrain.diffs
  )]
  total.dt.list[[type.csv]] <- type.total.dt
}
total.dt <- rbindlist(total.dt.list)
fwrite(total.dt, "figure-line-search-complexity.csv")
rfac <- 10
total.dt[, N:= 10^(round(log10(n.subtrain.diffs)*rfac)/rfac)]

L <- list(measurements=total.dt[maxIterations=="min.aum", data.table(
  iterations=mean(sum.iterations),
  min=min(sum.iterations),
  max=max(sum.iterations)
), by=.(expr.name=paste(aum.type, init.name), N)])
my_funs <- list(
  N=function(N)log10(N),
  "N \\log N"=function(N)log10(N) + log10(log(N)),
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec="iterations", fun.list=my_funs)
meas <- best[["measurements"]]
ref.dt <- best[["references"]]
ref.color <- "violet"
emp.color <- "black"
gg <- ggplot2::ggplot()+
  ggplot2::facet_grid(unit ~ expr.name, scales="free")+
  ggplot2::theme_bw()+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max),
    data=meas,
    fill=emp.color,
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical),
    size=2,
    color=emp.color,
    data=meas)+
  ggplot2::geom_line(ggplot2::aes(
    N, reference, group=fun.name),
    color=ref.color,
    size=1,
    data=ref.dt)+
  ggplot2::scale_y_log10("")+
  ggplot2::scale_x_log10()
if(requireNamespace("directlabels")){
  gg+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, label=fun.name),
      data=ref.dt,
      color=ref.color,
      method="bottom.polygons")
}else{
  gg
}

L <- list(measurements=total.dt[maxIterations=="min.aum", data.table(
  iterations=mean(sum.iterations),
  min=min(sum.iterations),
  max=max(sum.iterations)
), by=.(expr.name=paste(init.name), N)])
my_funs <- list(
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec="iterations", fun.list=my_funs)
meas <- best[["measurements"]]
ref.dt <- best[["references"]]
ref.color <- "violet"
emp.color <- "black"
gg <- ggplot2::ggplot()+
  ggplot2::facet_grid(unit ~ expr.name, scales="free")+
  ggplot2::theme_bw()+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max),
    data=meas,
    fill=emp.color,
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical),
    size=2,
    color=emp.color,
    data=meas)+
  ggplot2::geom_line(ggplot2::aes(
    N, reference, group=fun.name),
    color=ref.color,
    size=1,
    data=ref.dt)+
  ggplot2::scale_y_log10("")+
  ggplot2::scale_x_log10()
if(requireNamespace("directlabels")){
  gg+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, label=fun.name),
      data=ref.dt,
      color=ref.color,
      method="bottom.polygons")
}else{
  gg
}

L <- list(measurements=total.dt[maxIterations=="min.aum", data.table(
  sum.iterations=mean(sum.iterations),
  mean.iterations=mean(mean.iterations)
), by=.(expr.name=paste(aum.type, init.name), N)]
my_funs <- list(
  N=function(N)log10(N),
  "N \\log N"=function(N)log10(N) + log10(log(N)),
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(
  L, unit.col.vec=c("sum.iterations", "mean.iterations"), fun.list=my_funs)
meas <- best[["measurements"]]
ref.dt <- best[["references"]]
ref.color <- "violet"
emp.color <- "black"
gg <- ggplot2::ggplot()+
  ggplot2::facet_grid(unit ~ expr.name, scales="free")+
  ggplot2::theme_bw()+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical),
    size=2,
    color=emp.color,
    data=meas)+
  ggplot2::geom_line(ggplot2::aes(
    N, reference, group=fun.name),
    color=ref.color,
    size=1,
    data=ref.dt)+
  ggplot2::scale_y_log10("")+
  ggplot2::scale_x_log10()
if(requireNamespace("directlabels")){
  gg+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, label=fun.name),
      data=ref.dt,
      color=ref.color,
      method="bottom.polygons")
}else{
  gg
}

total.wide <- dcast(
  total.dt,
  N ~ .,
  value.var=c("sum.iterations", "steps"),
  fun.aggregate = list(median, min, max)
)[, expr.name := "line.search"]
L <- list(measurements=total.wide)
my_funs <- list(
  N=function(N)log10(N),
  "N \\log N"=function(N)log10(N) + log10(log(N)),
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec=c("sum.iterations_median", "steps_median"), fun.list=my_funs)
addUnit <- function(DT)DT[, Unit := sub("_median", "", unit)]
meas <- addUnit(best[["measurements"]])
ref.dt <- addUnit(best[["references"]])
ref.color <- "violet"
emp.color <- "black"
ribbon.dt <- nc::capture_melt_multiple(total.wide, Unit=".*", "_", column="min|max")
gg <- ggplot2::ggplot()+
  ggplot2::facet_grid(Unit ~ ., scales="free")+
  ggplot2::theme_bw()+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max),
    data=ribbon.dt,
    fill=emp.color,
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical),
    size=2,
    color=emp.color,
    data=meas)+
  ggplot2::geom_line(ggplot2::aes(
    N, reference, group=fun.name),
    color=ref.color,
    size=1,
    data=ref.dt)+
  ggplot2::scale_y_log10("")+
  ggplot2::scale_x_log10()
if(requireNamespace("directlabels")){
  gg+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, label=fun.name),
      data=ref.dt,
      color=ref.color,
      method="bottom.polygons")
}else{
  gg
}


value.var <- c("sum.iterations")
unit.col.vec <- paste0(value.var,"_median")
total.wide <- dcast(
  total.dt[maxIterations=="min.aum"],
  N + init.name ~ .,
  value.var=value.var,
  fun.aggregate = list(median, min, max, length)
)[, expr.name := init.name]
L <- list(measurements=total.wide)
my_funs <- list(
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec=unit.col.vec, fun.list=my_funs)
addUnit <- function(DT)DT[, Unit := sub("_median", "", unit)]
meas <- addUnit(best[["measurements"]])
ref.dt <- addUnit(best[["references"]])
ref.color <- "violet"
emp.color <- "black"
ribbon.dt <- nc::capture_melt_multiple(total.wide, Unit=".*", "_", column="min|max")
gg <- ggplot2::ggplot()+
  ggplot2::facet_grid(Unit ~ expr.name, scales="free")+
  ggplot2::theme_bw()+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max),
    data=ribbon.dt,
    fill=emp.color,
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical),
    size=2,
    color=emp.color,
    data=meas)+
  ggplot2::geom_line(ggplot2::aes(
    N, reference, group=fun.name),
    color=ref.color,
    size=1,
    data=ref.dt)+
  ggplot2::scale_y_log10("")+
  ggplot2::scale_x_log10()
dl <- if(requireNamespace("directlabels")){
  gg+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, label=fun.name),
      data=ref.dt,
      color=ref.color,
      method="bottom.polygons")
}else{
  gg
}
png('figure-line-search-complexity.png', width=8, height=4, units="in", res=200)
print(dl)
dev.off()


#analyze.
reg=batchtools::loadRegistry(registry.dir)

status.dt <- batchtools::getJobStatus(reg=reg)
status.dt[!is.na(error)]
status.dt[!is.na(done)]
done.ids <- status.dt[is.na(error) & !is.na(done), job.id]
data.list <- list()
for(done.i in seq_along(done.ids)){
  job.id <- done.ids[[done.i]]
  args.row <- args.dt[job.id]
  res <- batchtools::loadResult(job.id)
  for(data.type in names(res)){
    data.list[[data.type]][[done.i]] <- args.row[, .(aum.type, res[[data.type]])]
  }
}
for(data.type in names(data.list)){
  data.dt <- rbindlist(data.list[[data.type]])
  data.csv <- sprintf("figure-line-search-complexity-compare-%s.csv", data.type)
  data.table::fwrite(data.dt, data.csv)
}

data.csv.vec <- Sys.glob("figure-line-search-complexity-compare-*.csv")
dt.list <- list()
for(data.csv in data.csv.vec){
  data.type <- gsub(".csv|.*-", "", data.csv)
  dt.list[[data.type]] <- data.table::fread(data.csv)
}

q25 <- function(x)quantile(x,0.25)
q75 <- function(x)quantile(x,0.75)

comp.dt <- data.table::fread("figure-line-search-complexity.csv")
subtrain.sizes <- unique(comp.dt[, .(data.name, cv.type, test.fold, n.subtrain.diffs)])
rfac <- 5
iterations.tall <- dt.list[["iterations"]][, .(
  steps=.N,
  iterations=sum(iterations)
),
by=.(aum.type, data.name, cv.type, test.fold, seed, maxIterations.name)
][
  subtrain.sizes, on=.(data.name, cv.type, test.fold), nomatch=0L
][, B := as.integer(10^(round(log10(n.subtrain.diffs)*rfac)/rfac))][]
iterations.wide <- dcast(
  iterations.tall[!is.na(iterations)],
  B + maxIterations.name ~ .,
  list(median, min, max, q25, q75),
  value.var="iterations")
L <- list(measurements=iterations.wide[, data.table(
  iterations=iterations_median,
  N=B,
  expr.name=maxIterations.name)])
my_funs <- list(
  N=function(N)log10(N),
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec="iterations", fun.list=my_funs)
best$ref[, maxIterations.name := expr.name]
library(ggplot2)
ref.color <- "red"
gg <- ggplot()+
  facet_grid(. ~ maxIterations.name, labeller=label_both)+
  geom_ribbon(aes(
    B, ymin=iterations_q25, ymax=iterations_q75),
    alpha=0.5,
    data=iterations.wide)+
  geom_line(aes(
    B, iterations_median),
    data=iterations.wide)+
  geom_line(aes(
    N, reference, group=fun.name),
    color=ref.color,
    data=best$ref)+
  directlabels::geom_dl(aes(
    N, reference, group=fun.name, label=fun.name),
    color=ref.color,
    method="bottom.polygons",
    data=best$ref)+
  scale_x_log10(
    "B = breakpoints in subtrain set error functions")+
  scale_y_log10(
    "Total line search iterations\nuntil loss stops decreasing\nmedian and quartiles")
png("figure-line-search-complexity-compare-iterations-refs.png", width=9, height=3, units="in", res=200)
print(gg)
dev.off()

steps.wide <- dcast(
  iterations.tall,
  B + maxIterations.name ~ .,
  list(median, min, max, q25, q75),
  value.var="steps")
L <- list(measurements=steps.wide[, data.table(
  steps=steps_median,
  N=B,
  expr.name=maxIterations.name)])
my_funs <- list(
  N=function(N)log10(N),
  "\\log N"=function(N)log10(log(N)))
best <- atime::references_best(L, unit.col.vec="steps", fun.list=my_funs)
best$ref[, maxIterations.name := expr.name]
ref.color <- "red"
gg <- ggplot()+
  facet_grid(. ~ maxIterations.name, labeller=label_both)+
  geom_ribbon(aes(
    B, ymin=steps_q25, ymax=steps_q75),
    alpha=0.5,
    data=steps.wide)+
  geom_line(aes(
    B, steps_median),
    data=steps.wide)+
  geom_line(aes(
    N, reference, group=fun.name),
    color=ref.color,
    data=best$ref)+
  directlabels::geom_dl(aes(
    N, reference, group=fun.name, label=fun.name),
    color=ref.color,
    method="bottom.polygons",
    data=best$ref)+
  scale_x_log10(
    "B = breakpoints in subtrain set error functions")+
  scale_y_log10(
    "Gradient descent steps\nuntil loss stops decreasing\nmedian and quartiles")
png("figure-line-search-complexity-compare-steps-refs.png", width=9, height=3, units="in", res=200)
print(gg)
dev.off()


seconds.tall <- dt.list[["time"]][, seconds := elapsed.seconds][
  subtrain.sizes, on=.(data.name, cv.type, test.fold), nomatch=0L
][, B := as.integer(10^(round(log10(n.subtrain.diffs)*rfac)/rfac))][]
seconds.wide <- dcast(
  seconds.tall,
  B + maxIterations.name ~ .,
  list(median, min, max, q25, q75),
  value.var="seconds")
L <- list(measurements=seconds.wide[, data.table(
  seconds=seconds_median,
  N=B,
  expr.name=maxIterations.name)])
my_funs <- list(
  N=function(N)log10(N),
  "N^2"=function(N)2*log10(N))
best <- atime::references_best(L, unit.col.vec="seconds", fun.list=my_funs)
best$ref[, maxIterations.name := expr.name]
library(ggplot2)
ref.color <- "red"
gg <- ggplot()+
  facet_grid(. ~ maxIterations.name, labeller=label_both)+
  geom_ribbon(aes(
    B, ymin=seconds_q25, ymax=seconds_q75),
    alpha=0.5,
    data=seconds.wide)+
  geom_line(aes(
    B, seconds_median),
    data=seconds.wide)+
  geom_line(aes(
    N, reference, group=fun.name),
    color=ref.color,
    data=best$ref)+
  directlabels::geom_dl(aes(
    N, reference, group=fun.name, label=fun.name),
    color=ref.color,
    method="bottom.polygons",
    data=best$ref)+
  scale_x_log10(
    "B = breakpoints in subtrain set error functions")+
  scale_y_log10(
    "Total line search seconds\nuntil loss stops decreasing\nmedian and quartiles")
png("figure-line-search-complexity-compare-seconds-refs.png", width=9, height=3, units="in", res=200)
print(gg)
dev.off()

gg <- ggplot()+
  geom_line(aes(
    B, seconds_median, color=maxIterations.name),
    data=seconds.wide)+
  geom_ribbon(aes(
    B, ymin=seconds_q25, ymax=seconds_q75, fill=maxIterations.name),
    alpha=0.5,
    data=seconds.wide)+
  scale_x_log10(
    "B = breakpoints in subtrain set error functions",
    limits=c(NA,8000))+
  scale_y_log10(
    "Total line search seconds\nuntil loss stops decreasing\nmedian and quartiles")
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-line-search-complexity-compare-seconds.png", width=4.5, height=3, units="in", res=200)
print(dl)
dev.off()

max.valid.auc <- dt.list$sets[
  set=="validation",
  .SD[which.max(auc), .(auc)],
  by=.(aum.type, data.name, cv.type, test.fold, seed, init.name, maxIterations.name)]
max.valid.wide <- dcast(
  max.valid.auc,
  aum.type+data.name+cv.type+test.fold+seed+init.name ~ maxIterations.name,
  value.var="auc")
max.valid.wide[, hist(min.aum-exactL)]
max.valid.wide[, table(min.aum>exactL)]
max.valid.wide[, summary(min.aum-exactL)]
max.valid.wide[, summary(min.aum-exactQ)]
max.valid.wide[, summary(exactQ-exactL)]

ggplot()+
  geom_boxplot(aes(
    auc, maxIterations.name),
    data=max.valid.auc)

max.valid.folds <- dcast(
  max.valid.auc,
  aum.type+data.name+cv.type+test.fold+init.name ~ maxIterations.name,
  list(mean, length, min, max),
  value.var="auc")
max.valid.folds[auc_min_exactL < auc_max_exactL]
max.valid.folds[auc_min_exactQ < auc_max_exactQ]
max.valid.folds[auc_min_min.aum < auc_max_min.aum]
pdt <- max.valid.folds[init.name=="IntervalRegressionCV", {
  L <- t.test(auc_mean_grid, auc_mean_min.aum, alternative="two.sided", paired=TRUE)
  with(L, data.table(p.value, estimate))
}, by=.(aum.type, data.name, cv.type)][order(p.value)]
pdt[order(estimate)]
select.dt <- rbind(
  data.table(aum.type="rate", data.name="CTCF_TDH_ENCODE", cv.type="equal_labels"),
  data.table(aum.type="rate", data.name="systematic", cv.type="R-3.6.0-sequenceID"))
select.auc <- max.valid.auc[select.dt, on=.(aum.type, data.name, cv.type)]

ggplot()+
  geom_point(aes(
    auc, maxIterations.name),
    data=select.auc)+
  facet_grid(test.fold ~ aum.type+data.name+cv.type, labeller=label_both, scales="free")

select.dt <- rbind(
  ##data.table(aum.type="rate", data.name="CTCF_TDH_ENCODE", cv.type="equal_labels"),
  data.table(aum.type="count", data.name="H3K27ac-H3K4me3_TDHAM_BP", cv.type="equal_labels"),
  ##data.table(aum.type="rate", data.name="systematic", cv.type="chrom"),
  ##data.table(aum.type="rate", data.name="systematic", cv.type="sequenceID"),
  NULL)
select.auc <- max.valid.auc[select.dt, on=.(aum.type, data.name, cv.type)]
ggplot()+
  geom_point(aes(
    auc, maxIterations.name),
    shape=1,
    data=select.auc)+
  facet_grid(init.name ~ test.fold, labeller=label_both, scales="free")
