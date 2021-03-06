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
N <- 10
set.seed(1)
rand.pred.vec <- rnorm(N)
subtrain.output.vec <- rep(c(-1, 1), l=N)
subtrain.diff.count.dt <- aum::aum_diffs_binary(subtrain.output.vec, denominator="count")
subtrain.diff.rate.dt <- aum::aum_diffs_binary(subtrain.output.vec, denominator="rate")
library(data.table)
PairsDT <- function(output.vec){
  is.positive <- output.vec == 1
  data.table(expand.grid(
    positive=which(is.positive),
    negative=which(!is.positive)))
}
subtrain.pairs.dt <- PairsDT(subtrain.output.vec)
margin <- 1
## Note: for efficiency subtrain labels are assumed to be pre-computed
## in the enclosing environment, once before the optimization starts.
AUM <- function(pred.vec, diff.dt){
  L <- aum::aum(diff.dt, pred.vec)
  d <- L$derivative_mat
  non.diff <- abs(d[,1] - d[,2]) > 1e-6
  if(any(non.diff)){
    cat(sprintf("%d non-diff points\n", sum(non.diff)))
    print(d[non.diff, ])
  }
  with(L, list(gradient=derivative_mat[,1], loss=aum))
}
loss.list <- list(
  logistic=function(pred.vec, output.vec=subtrain.output.vec, ...){
    N <- length(pred.vec)
    list(
      gradient=-output.vec/(1+exp(output.vec*pred.vec))/N,
      loss=sum(log(1+exp(-output.vec*pred.vec)))/N)
  },
  aum.count=function(pred.vec, diff.count.dt=subtrain.diff.count.dt, ...){
    AUM(pred.vec, diff.count.dt)
  },
  aum.rate=function(pred.vec, diff.rate.dt=subtrain.diff.rate.dt, ...){
    AUM(pred.vec, diff.rate.dt)
  },
  squared.hinge.all.pairs=function(pred.vec, pairs.dt=subtrain.pairs.dt, ...){
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
    list(gradient=grad.dt$gradient/N.pairs, loss=sum(pairs.dt$diff.clipped^2)/N.pairs)
  }
)
result.list <- list()
for(loss.name in names(loss.list)){
  fun <- loss.list[[loss.name]]
  result.list[[loss.name]] <- fun(rand.pred.vec)
}
str(result.list)
sapply(result.list, "[[", "gradient")
out.loss.list <- list()
for(data.name in names(data.list)){
  input.output.list <- data.list[[data.name]]
  input.mat <- input.output.list[["input.mat"]]
  full.input.mat <- scale(input.mat)
  full.output.vec <- input.output.list[["output.vec"]]
  stopifnot(full.output.vec %in% c(-1, 1))
  set.seed(1)
  n.folds <- 4
  unique.folds <- 1:n.folds
  fold.vec <- sample(rep(unique.folds, l=length(full.output.vec)))
  for(validation.fold in unique.folds){
    is.set.list <- list(
      validation=fold.vec == validation.fold,
      subtrain=fold.vec != validation.fold)
    set.data.list <- list()
    for(set.name in names(is.set.list)){
      is.set <- is.set.list[[set.name]]
      output.vec <- full.output.vec[is.set]
      set.data.list[[set.name]] <- list(
        output.vec=output.vec,
        input.mat=full.input.mat[is.set,],
        diff.rate.dt=aum::aum_diffs_binary(output.vec, denominator="rate"),
        diff.count.dt=aum::aum_diffs_binary(output.vec, denominator="count"),
        pairs.dt=PairsDT(output.vec))
    }
    X.mat <- set.data.list$subtrain$input.mat
    for(loss.name in names(loss.list)){
      loss.grad.fun <- loss.list[[loss.name]]
      step.candidates <- 10^seq(-2, 2, by=0.25)
      for(step.size in step.candidates){
        set.seed(1)
        weight.vec <- rnorm(ncol(X.mat))
        done <- FALSE
        iteration <- 0
        prev.set.loss.vec <- rep(1e10, 2)
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
          set.loss.vec <- numeric()
          for(set.name in names(set.data.list)){
            set.data <- set.data.list[[set.name]]
            set.loss <- loss.for.weight(weight.vec, set.data)
            set.loss.vec[[set.name]] <- set.loss[["loss"]]
            roc.df <- WeightedROC::WeightedROC(set.loss[["pred"]], set.data[["output.vec"]])
            auc <- WeightedROC::WeightedAUC(roc.df)
            out.dt <- data.table(
              data.name, validation.fold, loss.name, step.size, iteration, set.name,
              auc,
              loss.value=set.loss$loss)
            for(aum.type in c("count", "rate")){
              diff.name <- paste0("diff.", aum.type, ".dt")
              aum.list <- aum::aum(set.data[[diff.name]], set.loss[["pred"]])
              out.col <- paste0("aum.", aum.type)
              out.dt[[out.col]] <- aum.list[["aum"]]
            }
            out.loss.list[[paste(
              data.name, validation.fold, loss.name, step.size, iteration, set.name
            )]] <- out.dt
          }#for(set.name
          diff.set.loss.vec <- set.loss.vec - prev.set.loss.vec
          max.inc.iterations <- 10
          valid.increasing.iterations <- if(!is.finite(diff.set.loss.vec[["validation"]])){
            max.inc.iterations
          }else if(0 <= diff.set.loss.vec[["validation"]]){
            cat(sprintf(
              "data=%s fold=%d loss=%s step=%f it=%d non-dec-iterations=%d\n",
              data.name, validation.fold, loss.name, step.size, iteration,
              valid.increasing.iterations))
            valid.increasing.iterations+1
          }else{
            0
          }
          if(
            max.inc.iterations <= valid.increasing.iterations
            || loss.after.step$step.size == 0
            || 1000 < iteration
          ){
            done <- TRUE
          }
          prev.set.loss.vec <- set.loss.vec
        }#while(!done
      }#for(step.size
    }#for(loss.name
  }#for(validation.fold
}
out.loss <- do.call(rbind, out.loss.list)
data.table::fwrite(out.loss, "figure-DNA-Sonar-subtrain-valid-data.csv")
system("gzip figure-DNA-Sonar-subtrain-valid-data.csv")
