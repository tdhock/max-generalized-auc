OneSeed <- function(seed){
  library(data.table)
  data(package="mlbench")
  data(Sonar, package="mlbench")
  data(DNA, package="mlbench")
  data.list <- list(
    Sonar=list(
      input.mat=as.matrix(Sonar[,1:60]),
      output.vec=ifelse(Sonar$Class=="R", 1, -1)),
    DNA=list(
      input.mat=ifelse(as.matrix(DNA[,1:180])==0, 0, 1),
      output.vec=ifelse(DNA$Class=="n", -1, 1)))
  data.name <- "Sonar"
  input.output.list <- data.list[[data.name]]
  input.mat <- input.output.list[["input.mat"]]
  full.input.mat <- scale(input.mat)
  full.output.vec <- input.output.list[["output.vec"]]
  stopifnot(full.output.vec %in% c(-1, 1))
  unique.sets <- c("subtrain", "validation", "test")
  PairsDT <- function(output.vec){
    is.positive <- output.vec == 1
    data.table(expand.grid(
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
      cat(sprintf("%d non-diff points\n", sum(non.diff)))
      print(d[non.diff, ])
    }
    with(L, list(gradient=rowMeans(derivative_mat), loss=aum))
  }
  loss.list <- list(
    logistic=function(pred.vec, output.vec, ...){
      Logistic(pred.vec, output.vec, 1/length(pred.vec))
    },
    logistic.weighted=
      function(pred.vec, output.vec,
               obs.weights=subtrain.obs.weights, ...){
        Logistic(pred.vec, output.vec, obs.weights)
      },
    aum.count=function(pred.vec, diff.count.dt, ...){
      AUM(pred.vec, diff.count.dt)
    },
    aum.rate=function(pred.vec, diff.rate.dt, ...){
      AUM(pred.vec, diff.rate.dt)
    },
    squared.hinge.all.pairs=function(pred.vec, pairs.dt, margin=1, ...){
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
  seed.dt.list <- list()
  set.seed(seed)
  set.vec <- sample(rep(unique.sets, l=length(full.output.vec)))
  set.data.list <- list()
  for(set.name in unique.sets){
    is.set <- set.vec == set.name
    output.vec <- full.output.vec[is.set]
    set.data.list[[set.name]] <- list(
      output.vec=output.vec,
      obs.weights=equal.class.weights(output.vec),
      input.mat=full.input.mat[is.set,],
      diff.rate.dt=aum::aum_diffs_binary(output.vec, denominator="rate"),
      diff.count.dt=aum::aum_diffs_binary(output.vec, denominator="count"),
      pairs.dt=PairsDT(output.vec))
  }
  X.mat <- set.data.list$subtrain$input.mat
  for(loss.name in names(loss.list)){
    loss.grad.fun <- loss.list[[loss.name]]
    for(step.size in 10^seq(-2,1,by=0.5)){
      cat(sprintf("seed=%d loss=%s step.size=%f\n", seed, loss.name, step.size))
      set.seed(1)
      weight.vec <- last.w <- rnorm(ncol(X.mat))
      done <- FALSE
      iteration <- 0
      while(!done){
        iteration <- iteration+1
        loss.for.weight <- function(w, set.data=set.data.list$subtrain){
          pred <- set.data$input.mat %*% w
          set.data$pred.vec <- pred
          out <- do.call(loss.grad.fun, set.data)
          out$pred <- pred
          out
        }
        loss.before.step <- loss.for.weight(weight.vec)
        direction <- -t(X.mat) %*% loss.before.step[["gradient"]]
        loss.for.step <- function(step.size){
          new.weight <- weight.vec + step.size * direction
          out <- loss.for.weight(new.weight)
          out$new.weight <- new.weight
          out$step.size <- step.size
          out
        }
        loss.after.step <- loss.for.step(step.size)
        weight.vec <- loss.after.step[["new.weight"]]
        diff.w <- sum(abs(weight.vec-last.w))
        last.w <- weight.vec
        diverged <- !is.finite(diff.w)
        if(!diverged){
          for(set.name in names(set.data.list)){
            set.data <- set.data.list[[set.name]]
            set.loss <- loss.for.weight(weight.vec, set.data)
            roc.df <- WeightedROC::WeightedROC(
              set.loss[["pred"]],
              set.data[["output.vec"]])
            auc <- WeightedROC::WeightedAUC(roc.df)
            out.dt <- data.table(
              seed,
              loss.name,
              step.size,
              iteration,
              set.name,
              auc,
              loss.value=set.loss$loss)
            for(aum.type in c("count", "rate")){
              diff.name <- paste0("diff.", aum.type, ".dt")
              out.dt[[paste0("aum.", aum.type)]] <- if(
                all(is.finite(set.loss[["pred"]]))
              ){
                aum.list <- aum::aum(set.data[[diff.name]], set.loss[["pred"]])
                aum.list[["aum"]]
              }else{
                NA
              }
            }
            seed.dt.list[[paste(
              seed,
              loss.name,
              step.size,
              iteration,
              set.name
            )]] <- out.dt
          }#for(set.name
        }#if(!diverged
        if(2000 < iteration || diverged || diff.w < 1e-6){
          done <- TRUE
        }
      }#while(!done
    }#for(step.size
  }#for(loss.name
  do.call(rbind, seed.dt.list)
}

LAPPLY <- future.apply::future_lapply
LAPPLY <- lapply
future::plan("multisession")
out.loss.list <- LAPPLY(1:10, OneSeed)
out.loss <- do.call(rbind, out.loss.list)

data.table::fwrite(out.loss, "figure-sonar-comparisons-data.csv")
