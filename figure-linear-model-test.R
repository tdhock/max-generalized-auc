library(data.table)
library(ggplot2)

testFold.vec <- Sys.glob("../neuroblastoma-data/data/*/cv/*/testFolds/*")

OneFold <- function(testFold.path){
  cv.path <- dirname(dirname(testFold.path))
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
  folds.csv <- Sys.glob(file.path(data.dir, "cv", "*", "folds.csv"))[1]
  folds.dt <- data.table::fread(folds.csv)
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set := rep(
    c("subtrain", "validation"), l=.N)]
  seqs.train <- folds.dt[["sequenceID"]]
  X.all <- scale(data.list$inputs[, -1])
  rownames(X.all) <- data.list$inputs$sequenceID
  X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
  set.list <- list()
  for(s in unique(folds.dt$set)){
    set.list[[s]] <- rownames(X.finite) %in% folds.dt[s==set, sequenceID]
  }
  X.list <- lapply(set.list, function(i)X.finite[i, ])
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
  for(seed in 1:4)for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    weight.vec <- int.weights[-1]
    intercept <- int.weights[1]
    computeAUM <- function(w, i, is.set){
      pred.pen.vec <- (X.finite %*% w) + i
      pred.dt <- data.table(
        sequenceID=rownames(pred.pen.vec),
        pred.log.lambda=as.numeric(pred.pen.vec))
      set.dt <- pred.dt[is.set]
      penaltyLearning::ROChange(
        data.list$evaluation, set.dt, "sequenceID")
    }
    for(iteration in 1:50){
      summary.dt.list <- list()
      set.roc.list <- list()
      for(set in names(set.list)){
        set.roc.list[[set]] <- computeAUM(weight.vec, intercept, set.list[[set]])
        summary.dt.list[[set]] <- with(set.roc.list[[set]], data.table(
          set,
          thresholds[threshold=="predicted"],
          auc,
          aum))
      }
      summary.dt <- do.call(rbind, summary.dt.list)
      iteration.dt.list[[paste(seed, init.name, iteration)]] <- data.table(
        seed, init.name, iteration, summary.dt)
      cat(sprintf(
        "it=%d seed=%d init=%s\n",
        iteration, seed, init.name))
      g.dt <- set.roc.list[["subtrain"]][["aum.grad"]]
      g.vec <- g.dt$lo
      direction.vec <- -t(X.list[["subtrain"]]) %*% g.vec
      take.step <- function(s){
        weight.vec + s*direction.vec
      }
      set.aum.list <- list()
      for(step.size in 10^seq(-10, 0, by=0.5)){
        new.weight.vec <- take.step(step.size)
        for(set in "subtrain"){
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
    }#iteration
  }#seed/init.name
  data.table(
    do.call(rbind, iteration.dt.list),
    data.name, cv.type, test.fold)
}

all.it.list <- list()
for(testFold.i in seq_along(testFold.vec)){
  fdir <- testFold.vec[testFold.i]
  out.csv <- file.path(fdir, "linear-model-aum.csv")
  all.it.list[[testFold.i]] <- if(file.exists(out.csv)){
    data.table::fread(out.csv)
  }else{
    cat(sprintf("%4d / %4d %s\n", testFold.i, length(testFold.vec), fdir))
    iteration.dt <- OneFold(fdir)
    data.table::fwrite(iteration.dt, out.csv)
    iteration.dt
  }
}
all.it <- do.call(rbind, all.it.list)

subtrain.it <- all.it[set=="subtrain"]
subtrain.it[, diff := c(NA, diff(aum)), by=.(init.name, data.name, test.fold, seed)]
subtrain.it[, .(init.name, data.name, test.fold, iteration, aum, diff)]
subtrain.it[diff>1e-6]
gg <- ggplot()+
  ggtitle("check if train AUM decreases")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_line(aes(
    iteration, aum,
    group=paste(seed, init.name)),
    data=subtrain.it)+
  facet_grid(init.name + data.name + test.fold ~ ., scales="free", labeller=label_both)
png(
  "figure-linear-model-test-aum-train-decreases.png",
  width=4, height=35, res=100, units="in")
print(gg)
dev.off()

validation.it <- all.it[set=="validation"]
ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  scale_y_log10()+
  geom_line(aes(
    iteration, aum, color=init.name,
    group=paste(seed, init.name)),
    data=validation.it)+
  geom_point(aes(
    iteration, aum, color=init.name,
    group=paste(seed, init.name)),
    data=validation.it[
     ,
       .SD[which.min(aum)],
       by=.(data.name, test.fold, init.name, seed)])+
  facet_grid(data.name + test.fold ~ ., scales="free", labeller=label_both)

valid.best.ids <- all.it[
  set=="validation",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, test.fold, init.name, seed)]
test.best.ids <- all.it[
  set=="test",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, test.fold, init.name, seed)]

test.it1 <- all.it[set=="test" & iteration==1]
test.selected <- all.it[set=="test"][valid.best.ids, on=names(valid.best.ids)]
test.best <- all.it[set=="test"][test.best.ids, on=names(test.best.ids)]
test.show <- rbind(
  data.table(iterations="initial", test.it1),
  data.table(iterations="best", test.best),
  data.table(iterations="selected", test.selected))
test.show[, Iterations := factor(iterations, c("initial", "selected", "best"))]
gg <- ggplot()+
  ggtitle("Test AUM, selected=min valid aum, best=min test aum, max it=50")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    aum, Iterations, color=factor(test.fold)),
    shape=1,
    data=test.show)+
  facet_grid(init.name ~ data.name, scales="free", labeller=label_both)
png(
  "figure-linear-model-test-compare-init.png",
  width=8, height=6, res=100, units="in")
print(gg)
dev.off()

test.show[, neg.auc := -auc]
test.show.tall <- melt(
  test.show[init.name=="IntervalRegressionCV"],
  measure.vars=c("neg.auc", "error.percent", "aum"))
test.iCV <- dcast(
  test.show.tall,
  data.name + test.fold + variable + seed ~ iterations)
gg <- ggplot()+
  ggtitle("Test set metrics, init=IntervalRegressionCV, selected worse than initial")+
  geom_point(aes(
    initial, selected, color=factor(seed)),
    shape=1,
    data=test.iCV)+
  geom_abline()+
  facet_wrap(~ data.name + variable, scales="free", labeller=label_both, ncol=3)
png(
  "figure-linear-model-test-initial-selected.png",
  width=10, height=6, res=100, units="in")
print(gg)
dev.off()

gg <- ggplot()+
  ggtitle("Test set metrics, init=IntervalRegressionCV, best about the same as initial")+
  geom_point(aes(
    initial, best, color=factor(seed)),
    shape=1,
    data=test.iCV)+
  geom_abline()+
  facet_wrap(~ data.name + variable, scales="free", labeller=label_both, ncol=3)
png(
  "figure-linear-model-test-initial-best.png",
  width=10, height=6, res=100, units="in")
print(gg)
dev.off()




