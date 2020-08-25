library(data.table)
library(ggplot2)

data.dir <- "../neuroblastoma-data/data/ATAC_JV_adipose"

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

folds.csv <- Sys.glob(file.path(data.dir, "cv", "*", "folds.csv"))[1]
folds.dt <- data.table::fread(folds.csv)
validation.fold <- 1
validation.ids <- folds.dt[fold==validation.fold, sequenceID]

X.all <- scale(data.list$inputs[, -1])
rownames(X.all) <- data.list$inputs$sequenceID
X.finite <- X.all[, apply(is.finite(X.all), 2, all)]

set.list <- list(
  validation=rownames(X.finite) %in% validation.ids)
set.list$train <- !set.list$validation
X.list <- lapply(set.list, function(i)X.finite[i, ])

y.train <- data.list[["outputs"]][
  !sequenceID %in% validation.ids,
  cbind(min.log.lambda, max.log.lambda)]
set.seed(1)
weight.vec <- rep(0, ncol(X.finite))
intercept <- 0

computeAUM <- function(w, i, is.set){
  pred.pen.vec <- (X.finite %*% w) + i
  pred.dt <- data.table(
    sequenceID=rownames(pred.pen.vec),
    pred.log.lambda=as.numeric(pred.pen.vec))
  set.dt <- pred.dt[is.set]
  penaltyLearning::ROChange(
    data.list$evaluation, set.dt, "sequenceID")
}

iteration.dt.list <- list()

for(iteration in 1:1000){
  if(! iteration %in% names(iteration.dt.list)){
    summary.dt.list <- list()
    set.roc.list <- list()
    for(set in names(set.list)){
      set.roc.list[[set]] <- computeAUM(weight.vec, intercept, set.list[[set]])
      summary.dt.list[[set]] <- with(set.roc.list[[set]], data.table(
        set,
        thresholds,
        aum))
    }
    summary.dt <- do.call(rbind, summary.dt.list)
    iteration.dt.list[[paste(iteration)]] <- data.table(
      iteration, summary.dt)
    print(iteration)
    g.dt <- set.roc.list[["train"]][["aum.grad"]]
    not.same <- g.dt[lo != hi]
    if(0 < nrow(not.same)){
      print(not.same)
      stop("not equal")
    }
    g.vec <- g.dt$lo
    direction.vec <- -t(X.list[["train"]]) %*% g.vec
    take.step <- function(s){
      weight.vec + s*direction.vec
    }
    set.aum.list <- list()
    for(step.size in 10^seq(-10, 0, by=0.5)){
      new.weight.vec <- take.step(step.size)
      for(set in "train"){
        set.roc <- computeAUM(new.weight.vec, 0, set.list[[set]])
        set.aum.list[[paste(step.size, set)]] <- data.table(
          step.size, set, aum=set.roc$aum,
          intercept=set.roc$thresholds[
            threshold=="min.error", (max.thresh+min.thresh)/2])
      }
    }
    set.aum <- do.call(rbind, set.aum.list)
    best.dt <- set.aum[, .SD[min(aum)==aum], by=set]
    ggplot()+
      geom_line(aes(
        step.size, aum),
        data=set.aum)+
      geom_point(aes(
        step.size, aum),
        data=best.dt)+
      geom_text(aes(
        step.size, aum, label=aum),
        vjust=0,
        data=best.dt)+
      scale_x_log10()+
      scale_y_log10()+
      facet_grid(set ~ ., scales="free")
    weight.vec <- take.step(best.dt[["step.size"]])
    intercept <- best.dt[["intercept"]]
  }
}

iteration.dt <- do.call(rbind, iteration.dt.list)
iteration.dt[set=="train", .(iteration, threshold, aum, errors)]
ggplot()+
  geom_line(aes(
    iteration, aum),
    data=iteration.dt[threshold=="predicted"])+
  facet_grid(set ~ ., scales="free")

ggplot()+
  geom_line(aes(
    iteration, errors),
    data=iteration.dt[threshold=="predicted" & iteration>1])+
  facet_grid(set ~ ., scales="free")
