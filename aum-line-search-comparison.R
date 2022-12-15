library(ggplot2)
library(data.table)
library(aum)
library(Rcpp)

data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
nb.err <- with(neuroblastomaProcessed$errors, data.frame(
  example=paste0(profile.id, ".", chromosome),
  min.lambda,
  max.lambda,
  fp, fn))
signal.features <- neuroblastomaProcessed$feature.mat[,c("log2.n","log.hall")]
n.noise <- 20
set.seed(1)
noise.features <- matrix(
  rnorm(n.noise*nrow(signal.features)),
  nrow(signal.features), n.noise)
X.sc <- scale(cbind(signal.features, noise.features))
keep <- apply(is.finite(X.sc), 2, all)
X.keep <- X.sc[,keep]
subtrain.i <- sample(1:nrow(X.keep), nrow(X.keep)/2)
index.list <- list(
  subtrain=subtrain.i,
  validation=-subtrain.i)
diff.list <- lapply(index.list, function(set.i){
  aum::aum_diffs_penalty(nb.err, rownames(X.keep)[set.i])
})
loss.dt.list <- list()
for (factor in 10^seq(1,4)) {
  for (search.type in c("exact", "exactq", "grid", "halving")) {
    weight.vec <- rep(0, ncol(X.keep))
    improvement <- old.aum <- Inf
    iteration <- 0
    X.keep <- X.sc[,keep]
    halving.step.size <- 1
    while(improvement > 1e-4){
      iteration <- iteration+1
      valid.list <- aum::aum(diff.list$validation, X.keep[index.list$validation,] %*% weight.vec)
    
      subtrain.aum <- Inf
      best.row <- NULL
      search.gradient.weight <- NULL
      ptm <- proc.time() # timer
      if (search.type == "exact" || search.type == "exactq") {
        quadratic.iterations <- nrow(diff.list$subtrain) * (nrow(diff.list$subtrain) - 1) / 2
        max.iterations <- if (search.type == "exact") nrow(diff.list$subtrain) * factor else quadratic.iterations

        nb.weight.search <- aum::aum_line_search_grid(
          diff.list$subtrain,
          feature.mat=X.keep[index.list$subtrain,],
          weight.vec=weight.vec,
          maxIterations=max.iterations)
        search.result <- nb.weight.search$line_search_result
        exact.dt <- data.table(search.result)
        exact.dt[, kink := .I/.N]
        best.row <- exact.dt[which.min(aum)]
        subtrain.aum <- nb.weight.search$aum
        search.gradient.weight <- nb.weight.search$gradient_weight
      } else if (search.type == "grid") {
        grid.step.size <- seq(0, max(1 / (factor + 1)), l = factor)
        feature.mat <- X.keep[index.list$subtrain,]
        pred.vec <- feature.mat %*% weight.vec
        subtrain.list <- aum::aum(diff.list$subtrain, pred.vec)
        subtrain.gradient.pred <- rowMeans(subtrain.list$derivative_mat)
        subtrain.gradient.weight <- t(feature.mat) %*% subtrain.gradient.pred
        subtrain.gradient <- feature.mat %*% subtrain.gradient.weight
    
        step.mat <- matrix(grid.step.size, length(pred.vec), length(grid.step.size), 
                           byrow = TRUE)
        pred.mat <- as.numeric(pred.vec) - step.mat * as.numeric(subtrain.gradient)
        grid.aum <- data.table(step.size=grid.step.size, aum = apply(pred.mat, 2, 
                               function(pred) aum::aum(diff.list$subtrain, pred)$aum))
        grid.aum[, kink := .I/.N]
        best.row <- grid.aum[which.min(aum)]
        subtrain.aum <- subtrain.list$aum
        search.gradient.weight <- subtrain.gradient.weight
      } else if (search.type == "halving") {
        feature.mat <- X.keep[index.list$subtrain,]
        pred.vec <- feature.mat %*% weight.vec
        subtrain.list <- aum::aum(diff.list$subtrain, pred.vec)
        subtrain.gradient.pred <- rowMeans(subtrain.list$derivative_mat)
        subtrain.gradient.weight <- t(feature.mat) %*% subtrain.gradient.pred
        subtrain.gradient <- feature.mat %*% subtrain.gradient.weight

        halving.steps <- factor * 10
        halving.aum.dt.list <- list()
        while({
          step.pred.vec <- as.numeric(pred.vec) - halving.step.size * -as.numeric(subtrain.gradient)
          step.list <- aum::aum(diff.list$subtrain, step.pred.vec)
          halving.aum.dt.list[[paste(halving.steps)]] <- data.table(step.size=halving.step.size, aum=step.list$aum)
          valid.list$aum < step.list$aum && halving.step.size != 0 && halving.steps >= 0
        }){
          cat(sprintf("  halving.step.size=%f aum=%f->%f\n", halving.step.size, valid.list$aum, step.list$aum))
          halving.steps <- halving.steps - 1
          halving.step.size <- halving.step.size/2
        }
        halving.aum.dt <- do.call(rbind, halving.aum.dt.list)
        halving.aum.dt[, kink := .I/.N]
        best.row <- halving.aum.dt[which.min(aum)]

        halving.step.size <- halving.step.size*2
        subtrain.aum <- step.list$aum
        search.gradient.weight <- subtrain.gradient.weight
      }
      elapsed.time <- (proc.time() - ptm)[["elapsed"]]
      
      loss.dt.list[[paste0(search.type, ".", factor, ".", iteration)]] <- data.table(
        iteration,
        factor,
        search.type,
        time=elapsed.time,
        set=c("subtrain", "validation"),
        aum=c(subtrain.aum, valid.list$aum))
    
      if(interactive())cat(sprintf(
        "iteration=%4d search.type=%s factor=%d aum=%.4f step=%.4f kink=%f\n",
        iteration, search.type, factor, best.row$aum, best.row$step.size, best.row$kink))
      improvement <- old.aum-best.row$aum
      old.aum <- best.row$aum
      weight.vec <- weight.vec-best.row$step.size*search.gradient.weight
    }
  }
}
loss.dt <- do.call(rbind, loss.dt.list)
min.dt <- loss.dt[, .SD[which.min(aum)], by=list(search.type, set, factor)]
best.dt <- loss.dt[, .SD[which.min(aum)], by=list(search.type, set)]

ggplot() +
  geom_line(aes(
    iteration, aum, color=set),
    data=loss.dt) +
  geom_point(aes(
    iteration, aum, color=set),
    shape=1,
    data=min.dt) +
  geom_point(aes(iteration, aum), color="darkblue", data=best.dt) +
  geom_vline(aes(xintercept = iteration), color="grey", data=min.dt) +
  geom_hline(aes(yintercept = aum), color="grey", data=min.dt) +
  scale_y_log10() +
  scale_x_log10() +
  #facet_grid(search.type ~ factor)
  facet_grid(search.type ~ factor, scales="free")
  #facet_wrap("factor", scales="free")

# get the best graphs for every search type
best.valid.dt <- best.dt[set == "validation"][,.SD[which.min(aum)],by=list(search.type)]
loss.best.dt <- loss.dt[best.valid.dt, on=.(search.type, factor, set)]
ggplot() +
  geom_line(aes(
    iteration, aum, color=search.type),
    data=loss.best.dt) +
  scale_y_log10() +
  scale_x_log10()

total.time <- aggregate(loss.dt$time, by=list(search.type=loss.dt$search.type), FUN=sum)

ggplot(data = total.time, aes(y=x, x=search.type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(x)), vjust=1.6, color="white", size=3.5) +
  ylab("total time (seconds)") + xlab("line search")

