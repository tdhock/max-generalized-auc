source("packages.R")

PairsDT <- function(output.vec){
  is.positive <- output.vec == 1
  data.table::data.table(expand.grid(
    positive=which(is.positive),
    negative=which(!is.positive)))
}
equal.class.weights <- function(output.vec){
  otab <- table(output.vec)
  as.numeric(1/otab[paste(output.vec)])
}
Logistic <- function(pred.vec, output.vec, obs.weights){
  list(
    gradient=-obs.weights*output.vec/(1+exp(output.vec*pred.vec)),
    loss=sum(obs.weights*log(1+exp(-output.vec*pred.vec))))
}
AUM <- function(pred.vec, diff.dt){
  L <- aum::aum(diff.dt, pred.vec)
  d <- L$derivative_mat
  non.diff <- abs(d[,1] - d[,2]) > 1e-6
  if(any(non.diff)){
    ## Some non-differentiable points that were actually observed!
    ## data=DNA fold=1 loss=aum.rate step=0.001000
    ##              [,1]         [,2]
    ## [1,] -0.001956947 -0.001175589
    ## data=DNA fold=1 loss=aum.rate step=1000.000000
    ##               [,1] [,2]
    ## [1,] -0.0006463963    0
    cat(sprintf("%d non-diff points\n", sum(non.diff)))
    print(d[non.diff, ])
  }
  ## ifelse( derivative_mat[,1] == 0 | derivative_mat[,2] == 0, 0, ??
  with(L, list(
    gradient=(derivative_mat[,1]+derivative_mat[,2])/2,
    loss=aum))
}
zip.X.list <- list()
zip.y.list <- list()
for(set in c("train", "test")){
  f <- sprintf("zip.%s.gz", set)
  if(!file.exists(f)){
    u <- paste0("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/", f)
    download.file(u, f)
  }
  zip.dt <- data.table::fread(f)
  y.vec <- zip.dt[[1]]
  is.01 <- y.vec %in% 0:1
  y01.dt <- data.table(label=y.vec[is.01])
  y01.dt[, cum := 1:.N, by=label]
  max.dt <- y01.dt[, .(max=max(cum)), by=label]
  keep <- y01.dt$cum <= min(max.dt[["max"]])
  zip.y.list[[set]] <- y01.dt[keep, label]
  zip.X.list[[set]] <- as.matrix(zip.dt[is.01, -1, with=FALSE][keep,])
}
(y.tab <- sapply(zip.y.list, table))

train.set.list <- list(
  full=list(X=zip.X.list[["train"]], y=zip.y.list[["train"]]))
prop.pos.vec <- some.props <- c(0.01, 0.05, 0.5)
##want p/(p + n) = 0.05 => 0.05*(p+n) = p => 0.05p + 0.05n = p => 0.05n = 0.95p => p = 0.05 / 0.95n
min.prop.pos <- min(prop.pos.vec)
min.n.pos <- as.integer(min.prop.pos/(1-min.prop.pos) * y.tab["0", "train"])
min.total <- min.n.pos + y.tab["0", "train"]
c(min.n.pos, y.tab["0", "train"])/min.total
N.obs <- 1000
train.y.dt <- data.table(label=zip.y.list[["train"]])
train.y.dt[, i := 1:.N]
test.y <- zip.y.list[["test"]]
test.X <- zip.X.list[["test"]]
result.dt.list <- list()
selected.dt.list <- list()
for(prop.pos in prop.pos.vec){
  prop.dt <- rbind(
    data.table(prop=prop.pos, label=1),
    data.table(prop=1-prop.pos, label=0))
  prop.dt[, class.N := as.integer(N.obs*prop) ]
  prop.dt[, weight := 1/class.N]
  for(seed in 1:10){
    cat(sprintf("prop=%f seed=%d\n", prop.pos, seed))
    set.seed(seed)
    index.dt <- prop.dt[train.y.dt, on="label"][, .(
      i=.SD[sample(1:.N), i[1:class.N] ]
    ), by=.(label, weight, class.N)]
    seed.i <- index.dt[["i"]]
    seed.y <- zip.y.list[["train"]][seed.i]
    seed.X <- zip.X.list[["train"]][seed.i,]
    weight.list <- list(
      identity=rep(1, length(seed.y)),
      balanced=index.dt[["weight"]])
    pred.list <- list()
    for(weight.name in names(weight.list)){
      weight.vec <- weight.list[[weight.name]]
      fit <- glmnet::cv.glmnet(seed.X, seed.y, weight.vec, family="binomial")
      seed.pred <- predict(fit, test.X)
      pred.list[[paste0("cv.glmnet.", weight.name)]] <- seed.pred
    }
    y.tilde <- ifelse(seed.y==0, -1, 1)
    is.validation <- rep(c(TRUE, FALSE), l=length(y.tilde))
    set.list <- list(subtrain=!is.validation, validation=is.validation)
    data.by.set <- list()
    for(set.name in names(set.list)){
      is.set <- set.list[[set.name]]
      data.by.set[[set.name]] <- list(X=seed.X[is.set,], y=y.tilde[is.set])
    }
    is.subtrain <- !is.validation
    y.subtrain <- y.tilde[is.subtrain]
    X.subtrain <- seed.X[is.subtrain,]
    diff.rate.dt <- aum::aum_diffs_binary(y.subtrain, denominator="rate")
    diff.count.dt <- aum::aum_diffs_binary(y.subtrain, denominator="count")
    pairs.dt <- PairsDT(y.subtrain)
    loss.list <- list(
      logistic=function(pred.vec){
        Logistic(pred.vec, y.subtrain, 1/length(pred.vec))
      },
      logistic.weighted=function(pred.vec){
          Logistic(pred.vec, y.subtrain, index.dt[is.subtrain, weight])
        },
      aum.count=function(pred.vec){
        AUM(pred.vec, diff.count.dt)
      },
      aum.rate=function(pred.vec){
        AUM(pred.vec, diff.rate.dt)
      },
      squared.hinge.all.pairs=function(pred.vec, margin=1){
        pairs.dt[, diff := pred.vec[positive]-pred.vec[negative]-margin]
        pairs.dt[, diff.clipped := ifelse(diff<0, diff, 0)]
        pairs.tall <- data.table::melt(
          pairs.dt,
          measure.vars=c("positive", "negative"),
          value.name="pred.i",
          variable.name="label")
        ## d/dx (x - y - m)^2 = x - y - m
        ## d/dy (x - y - m)^2 = -(x - y - m)
        pairs.tall[, grad.sign := ifelse(label=="positive", 1, -1)]
        N.pairs <- nrow(pairs.dt)
        grad.dt <- pairs.tall[, .(
          gradient=sum(grad.sign*diff.clipped)
        ), keyby=pred.i]
        list(
          gradient=grad.dt$gradient/N.pairs,
          loss=sum(pairs.dt$diff.clipped^2)/N.pairs)
      }
    )
    for(loss.name in names(loss.list)){
      loss.fun <- loss.list[[loss.name]]
      step.candidates <- 10^seq(-2, 2, by=0.5)
      selection.dt.list <- list()
      test.pred.list <- list()
      for(step.size in step.candidates){
        weight.vec <- rnorm(ncol(seed.X))
        for(iteration in 1:1000){
          pred.vec <- X.subtrain %*% weight.vec
          loss.info <- loss.fun(pred.vec)
          direction <- -t(X.subtrain) %*% loss.info[["gradient"]]
          weight.vec <- weight.vec + step.size * direction
          test.pred.list[[paste(step.size, iteration)]] <- test.X %*% weight.vec
          for(set.name in "validation"){
            Xy.list <- data.by.set[[set.name]]
            set.pred.vec <- Xy.list[["X"]] %*% weight.vec
            roc.df <- WeightedROC::WeightedROC(set.pred.vec, Xy.list[["y"]])
            auc <- WeightedROC::WeightedAUC(roc.df)
            selection.dt.list[[paste(step.size, iteration, set.name)]] <-
              data.table(step.size, iteration, set.name, auc)
          }#set.name
        }#iteration
      }#step.size
      selection.dt <- do.call(rbind, selection.dt.list)
      selected <- selection.dt[set.name=="validation"][which.max(auc)]
      selected.dt.list[[paste(prop.pos, seed, loss.name)]] <- data.table(
        prop.pos, seed, loss.name, selected)
      pred.list[[loss.name]] <- selected[
      , test.pred.list[[paste(step.size, iteration)]] ]
    }#loss.name
    for(model in names(pred.list)){
      seed.pred <- pred.list[[model]]
      roc.df <- WeightedROC::WeightedROC(seed.pred, test.y)
      seed.pred.class <- ifelse(0<seed.pred, 1, 0)
      accuracy <- mean(seed.pred.class == test.y)
      auc <- WeightedROC::WeightedAUC(roc.df)
      result.dt.list[[paste(prop.pos, seed, model)]] <- data.table(
        prop.pos, seed, model, accuracy, auc)
    }
  }#seed
}#prop.pos
(result.dt <- do.call(rbind, result.dt.list))
(selected.dt <- do.call(rbind, selected.dt.list))

saveRDS(list(result=result.dt, selected=selected.dt, N.obs=nrow(seed.X)), file="figure-unbalanced-grad-desc-data.rds")

