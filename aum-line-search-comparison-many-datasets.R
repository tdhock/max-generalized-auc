library(ggplot2)
library(data.table)
library(aum)
library(Rcpp)
library(dplyr)

folds.dt <- fread("../feature-learning-benchmark/labeled_problems_folds.csv")
addMeta <- function(dt){
  dt[, set.name := sub("/.*", "", prob.dir)]
  dt[, problem := sub(".*/", "", prob.dir)]
  dt[folds.dt, on=list(set.name, problem)]
}
errors.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_errors.csv"))
possible.dt <- addMeta(fread("../feature-learning-benchmark/labeled_problems_possible_errors.csv"))


auc.improved.list <- list()

fold.possible <- unique(folds.dt[, .(set.name, fold)])
i.possible <- 1:nrow(fold.possible)
N.possible <- paste(i.possible, "improved")
i.todo <- i.possible[!N.possible %in% names(auc.improved.list)]
biggest.step <- 0.1
completed.datasets.list <- list()
for(i in seq_along(i.todo)){
  test.fold.i <- i.todo[[i]]
  cat(sprintf("%4d / %4d test folds TODO=%d\n", i, length(i.todo), test.fold.i))
  test.fold.info <- fold.possible[test.fold.i]
  dataset <- paste0(test.fold.info$set.name, test.fold.info$fold)
  test.fold.errors <- errors.dt[test.fold.info, on=.(set.name, fold)]
  test.fold.errors[, min.log.lambda := min.log.penalty]
  test.fold.errors[, max.log.lambda := max.log.penalty]
  test.fold.errors[, seg.i := cumsum(
    c(1, diff(fp)!=0 | diff(fn) != 0)), by=.(prob.dir)]
  possible.errors <- possible.dt[test.fold.errors, on=list(
    set.name, fold, prob.dir)][, possible.fn := possible.tp]
  possible.segs <- possible.errors[, .(
    min.log.lambda=min(min.log.lambda),
    max.log.lambda=max(max.log.lambda)
  ), by=.(
    prob.dir, seg.i, fp, fn, errors, possible.fp, possible.fn, labels
  )][, `:=`(
    min.lambda = exp(min.log.lambda),
    example=prob.dir
  )]
  ## Check for non-zero at end of err fun.
  possible.segs[min.log.lambda == -Inf & fn > 0]
  possible.segs[min.log.lambda == Inf & fp > 0]
  test.fold.targets <- penaltyLearning::targetIntervals(
    possible.segs, "prob.dir")
  prob.ord <- test.fold.targets$prob.dir
  aum.diffs <- aum::aum_diffs_penalty(possible.segs, prob.ord)
  min.err.pred.dt <- test.fold.targets[, data.table(
    prob.dir,
    pred.log.lambda=fcase(
      min.log.lambda>-Inf & max.log.lambda==Inf, min.log.lambda+1, 
      min.log.lambda==-Inf & max.log.lambda<Inf, max.log.lambda-1,
      min.log.lambda>-Inf & max.log.lambda<Inf, (min.log.lambda+max.log.lambda)/2,
      min.log.lambda==-Inf & max.log.lambda==Inf, 0)
  )]
  possible.segs[prob.dir==prob.ord[1], .(fp,fn,min.log.lambda)]
  aum.diffs[example==0]
  init.list <- list(
    min.error=-min.err.pred.dt$pred.log.lambda,
    zero=rep(0, nrow(min.err.pred.dt)))
  
  current.pred <- initial.pred <- init.list$min.error
  
  maxIterations <- nrow(aum.diffs)*(nrow(aum.diffs)-1) / 2
  improvement <- 1
  initial.aum <- aum::aum(aum.diffs, initial.pred)$aum
  # while we're still improving the total AUM
  steps <- 0
  iterations.list <- list()
  while (improvement > 0.00001) {
    # line search for where we're at
    weight.search <- aum::aum_line_search(
      aum.diffs,
      pred.vec=current.pred
      #maxIterations=maxIterations
    )
    # plot(weight.search)
    line.search.result <- data.table(weight.search$line_search_result)
    gradient <- weight.search$gradient
    best.row <- line.search.result[which.min(aum)]

    most.intersections <- max(weight.search$line_search_result$intersections)
    if (most.intersections > 1) {
      cat(sprintf("      most.intersections=%d on step=%d\n", most.intersections, steps))
    }
    
    # take a step in the descent direction
    step.pred <- current.pred - best.row$step.size * gradient
    # calculate new aum and check the difference
    step.aum <- aum::aum(aum.diffs, step.pred)
    improvement <- weight.search$aum - step.aum$aum

    steps <- steps + 1
    current.pred <- step.pred
    iterations.list[[paste0(steps)]] <- data.table(
      step=steps,
      aum=best.row$aum,
      auc=best.row$auc,
      set.name=test.fold.info$set.name,
      fold=as.factor(test.fold.info$fold)
    )
  }
  cat(sprintf("      aum=%f->%f in %d steps\n", initial.aum, step.aum$aum, steps))
  completed.datasets.list[[paste0(dataset)]] <- do.call(rbind, iterations.list)
}
completed.datasets <- do.call(rbind, completed.datasets.list)

ggplot(completed.datasets) +
  geom_line(aes(x=step, y=aum, color=fold)) +
  facet_wrap(set.name ~., scales="free")


# OLD CODE from aum-line-search-comparison.R
data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
data(notConverging, package="penaltyLearning", envir=environment())
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
# quadratic amount of iterations
final.max.iterations <- nrow(diff.list$subtrain) * (nrow(diff.list$subtrain) - 1) / 2
loss.dt.list <- list()
best.dt.list <- list()
10^c(2, 2.5, 3, 3.5, 4, 4.5, 5, 6)
for (factor in 10^c(3, 4, 5, 6)) {
  last.iteration <- factor == 10^6
  for (search.type in c("exact", "hybrid03", "grid")) {
    if (search.type != "exact" && search.type != "hybrid" && last.iteration) break
    weight.vec <- rep(0, ncol(X.keep))
    improvement <- old.aum <- Inf
    iteration <- 0
    X.keep <- X.sc[,keep]
    while(improvement > 1e-4){
      iteration <- iteration+1
      valid.list <- aum::aum(diff.list$validation, X.keep[index.list$validation,] %*% weight.vec)
    
      subtrain.aum <- Inf
      best.row <- NULL
      search.gradient.weight <- NULL
      ptm <- proc.time() # timer
      factor.label <- factor
      if (search.type == "exact") {
        max.iterations <-  factor
        # if we're on the last factor iteration, use the quadratic amount instead
        if (last.iteration) { #special case this one
          max.iterations <- final.max.iterations
          factor.label <- max.iterations
        }

        nb.weight.search <- aum::aum_line_search(
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
        #grid.step.size <- 10^seq(-7, -1, l = factor)
        grid.step.size <- 10^seq(-7, -1, l = 10)
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
      } else if (startsWith(search.type, "hybrid")) {
        max.iterations <- factor
        # if we're on the last factor iteration, use the quadratic amount instead
        if (last.iteration) { #special case this one
          max.iterations <- final.max.iterations
          factor.label <- max.iterations
        }
        
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
        if (best.row$kink == 1) {
          # if kink == 1, we have chosen the very last step size we looked at.
          # run a grid search for where we're at to find a larger step.size
          
          # get the number of grid points we want to check
          grid.point.count <- as.numeric(gsub("hybrid", "", search.type))
          # start at 10 * the best step size we found
          grid.step.size <- best.row$step.size * 2
          feature.mat <- X.keep[index.list$subtrain,]
          pred.vec <- feature.mat %*% weight.vec
          subtrain.list <- aum::aum(diff.list$subtrain, pred.vec)
          subtrain.gradient.pred <- rowMeans(subtrain.list$derivative_mat)
          subtrain.gradient.weight <- t(feature.mat) %*% subtrain.gradient.pred
          subtrain.gradient <- feature.mat %*% subtrain.gradient.weight
          
          # new prediction vector stepping in the direction of the negative gradient
          points.searched <- 0
          repeat {
            points.searched <- points.searched + 1
            step.pred <- as.numeric(pred.vec) - grid.step.size * as.numeric(subtrain.gradient)
            step.aum <- aum(diff.list$subtrain, step.pred)$aum
            
            if (step.aum < best.row$aum) {
              best.row <- data.table(
                step.size=grid.step.size,
                aum=step.aum, 
                kink=1+(points.searched / 10)
              )
              subtrain.aum <- subtrain.list$aum
              search.gradient.weight <- subtrain.gradient.weight
            }
            
            if (step.aum > best.row$aum || points.searched >= grid.point.count) {
              break
            } else {
              grid.step.size <- grid.step.size * 10
            }
          }
          
          # OLD
          # step.mat <- matrix(grid.step.size, length(pred.vec), length(grid.step.size), 
          #                    byrow = TRUE)
          # pred.mat <- as.numeric(pred.vec) - step.mat * as.numeric(subtrain.gradient)
          # grid.aum <- data.table(step.size=grid.step.size, aum = apply(pred.mat, 2, 
          #                 function(pred) aum::aum(diff.list$subtrain, pred)$aum))
          # grid.aum[, kink := .I/.N]
          # best.grid.row <- grid.aum[which.min(aum)]
          
          # if (best.grid.row$aum < best.row$aum) {
          #   best.row <- best.grid.row
          #   subtrain.aum <- subtrain.list$aum
          #   search.gradient.weight <- subtrain.gradient.weight
          # }
        }
      } else if (search.type == "constant") {
        feature.mat <- X.keep[index.list$subtrain,]
        pred.vec <- feature.mat %*% weight.vec
        subtrain.list <- aum::aum(diff.list$subtrain, pred.vec)
        subtrain.gradient.pred <- rowMeans(subtrain.list$derivative_mat)
        subtrain.gradient.weight <- t(feature.mat) %*% subtrain.gradient.pred
        subtrain.gradient <- feature.mat %*% subtrain.gradient.weight
        
        constant.step.size <- 1 / (factor) #* 100)
        new.pred <- as.numeric(pred.vec) - constant.step.size * as.numeric(subtrain.gradient)
        constant.aum <- data.table(step.size=constant.step.size, aum=aum::aum(diff.list$subtrain, new.pred)$aum)
        constant.aum[, kink := .I/.N]
        best.row <- constant.aum
        subtrain.aum <- subtrain.list$aum
        search.gradient.weight <- subtrain.gradient.weight
      } else if (search.type == "halving") {
        feature.mat <- X.keep[index.list$subtrain,]
        pred.vec <- feature.mat %*% weight.vec
        subtrain.list <- aum::aum(diff.list$subtrain, pred.vec)
        subtrain.gradient.pred <- rowMeans(subtrain.list$derivative_mat)
        subtrain.gradient.weight <- t(feature.mat) %*% subtrain.gradient.pred
        subtrain.gradient <- feature.mat %*% subtrain.gradient.weight
        
        halving.step.size <- 1
        halving.steps <- factor * 10
        halving.aum.dt.list <- list()
        while({
          step.pred.vec <- as.numeric(pred.vec) - halving.step.size * as.numeric(subtrain.gradient)
          step.list <- aum::aum(diff.list$subtrain, step.pred.vec)
          halving.aum.dt.list[[paste(halving.steps)]] <- data.table(step.size=halving.step.size, aum=step.list$aum)
          # valid.list$aum < step.list$aum &&
          halving.step.size > 1 / (factor) && halving.steps >= 0
        }){
          cat(sprintf("  halving.step.size=%f aum=%f->%f\n", halving.step.size, valid.list$aum, step.list$aum))
          halving.steps <- halving.steps - 1
          halving.step.size <- halving.step.size/2
        }
        cat(sprintf("  halving.step.size=%f aum=%f->%f (last)\n", halving.step.size, valid.list$aum, step.list$aum))
        halving.aum.dt <- do.call(rbind, halving.aum.dt.list)
        halving.aum.dt[, kink := .I/.N]
        best.row <- halving.aum.dt[which.min(aum)]

        halving.step.size <- halving.step.size*2
        subtrain.aum <- step.list$aum
        search.gradient.weight <- subtrain.gradient.weight
      }
      elapsed.time <- (proc.time() - ptm)[["elapsed"]]
      
     # if ((search.type != "constant" && search.type != "halving") || old.aum-best.row$aum > 0) {
        loss.dt.list[[paste0(search.type, ".", factor, ".", iteration)]] <- data.table(
          iteration,
          factor=factor.label,
          search.type,
          time=elapsed.time,
          set=c("subtrain", "validation"),
          aum=c(subtrain.aum, valid.list$aum))
    # }
      
      # collect the best step sizes found
      best.dt.list[[paste0(search.type, ".", factor, ".", iteration)]] <- data.table(
        step=iteration,
        points=factor.label,
        search.type,
        aum=best.row$aum,
        step.size=best.row$step.size,
        kink=best.row$kink
      )
    
      if(interactive())cat(sprintf(
        "iteration=%4d search.type=%s factor=%.1f aum=%.4f->%.4f step=%.4f kink=%f\n",
        iteration, search.type, factor.label, old.aum, best.row$aum, best.row$step.size, best.row$kink))
      improvement <- old.aum-best.row$aum
      old.aum <- best.row$aum
      weight.vec <- weight.vec-best.row$step.size*search.gradient.weight
    }
  }
}
cat("done!")
loss.dt <- do.call(rbind, loss.dt.list)
best.dt <- do.call(rbind, best.dt.list)
min.dt <- loss.dt[, .SD[which.min(aum)], by=list(search.type, set, factor)]
min.best.dt <- loss.dt[, .SD[which.min(aum)], by=list(search.type, set)]

# big plot showing grid of search types and aum for subtrain & validation
ggplot() +
  geom_vline(aes(xintercept = iteration), color="darkgrey", data=min.dt) +
  geom_hline(aes(yintercept = aum), color="darkgrey", data=min.dt) +
  geom_line(aes(
    iteration, aum, color=set),
    size=1.3,
    data=loss.dt) +
  geom_point(aes(
    iteration, aum),
    color="lightblue",
    data=min.dt) +
  geom_point(aes(iteration, aum), color="darkblue", size=3, data=min.best.dt) +
  scale_fill_brewer(type = "seq", palette = "Pastel1") +
  scale_y_log10() +
  #facet_grid(search.type ~ factor)
  xlim(0, 1000) +
  scale_x_log10() +
  facet_grid(search.type ~ factor, scales="free")
  #facet_grid(factor ~ search.type, scales="free")
  #theme_light()
  #facet_wrap("factor", scales="free")

# get the best graphs for every search type
min.best.valid.dt <- min.best.dt[set == "validation"][,.SD[which.min(aum)],by=list(search.type)]
min.best.subtrain.dt <- min.best.dt[set == "subtrain"][,.SD[which.min(aum)],by=list(search.type)]
loss.best.valid.dt <- loss.dt[min.best.valid.dt, on=.(search.type, factor, set)]
loss.best.subtrain.dt <- loss.dt[min.best.subtrain.dt, on=.(search.type, factor, set)]
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

loss.dt %>%
  group_by(search.type, factor) %>%
  summarize(total.time = sum(time)) %>%
#  filter(search.type == "exact") %>%
  ggplot(aes(x=search.type, y=total.time, fill = search.type)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label=round(total.time)), vjust=1.6, color="white", size=3.5) +
  facet_grid(. ~ factor) +
  #scale_y_log10() +
  ylab("total time (ms)") + xlab("line search")

loss.dt %>%
  group_by(search.type, factor) %>%
  summarize(total.time = sum(time)) %>%
  ggplot(aes(color=search.type, x=factor, y=total.time)) +
  geom_line(size=2) +
  scale_y_log10() +
  scale_x_log10() +
  ylab("total time (seconds)") +
  xlab("points per step") +
  ggtitle("Total time for gradient descent")

loss.dt %>%
  group_by(factor, search.type) %>%
  summarize(average.time = mean(time)) %>%
  ggplot(aes(x=factor, y=average.time, color=search.type)) +
  geom_line(size=2) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("average time per gradient descent step (seconds)") +
  xlab("points per step (exact=I, grid=G)") +
  ggtitle("Average time per step of gradient descent")

best.dt %>%
  group_by(points, search.type, step) %>%
  ggplot() +
  geom_line(aes(x=step, y=step.size, color=search.type)) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(.~points, scales = "free") +
  ggtitle("Best step size chosen for each step of gradient descent")

best.dt %>%
  group_by(points, search.type, step) %>%
  ggplot() +
#  geom_line(aes(x=step, y=step.size, color=search.type)) +
  geom_line(aes(x=step, y=kink, color=search.type), linetype=1) +
#  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(.~points, scales = "free")

best.dt %>%
  group_by(points, search.type, step) %>%
  ggplot() +
  #geom_line(aes(x=step, y=step.size, color=search.type)) +
  geom_line(aes(x=step, y=kink, color=search.type), linetype=1) +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(search.type~points, scales = "free") +
  ggtitle("Step size chosen for hybrid search")

ggplot() +
  geom_line(aes(
    iteration, aum, color=set),
    data=loss.best.valid.dt[search.type == "exact"]) +
  geom_line(aes(
    iteration, aum, color=set),
    data=loss.best.subtrain.dt[search.type == "exact"]) +
  geom_point(aes(iteration, aum), color="darkblue", size=2, data=loss.best.valid.dt[search.type == "exact"][which.min(aum)]) +
  ggtitle("AUM Gradient Descent") +
  xlab("step") +
  ylab("AUM")


ggplot() +
  ggtitle("Runs with lowest overall time") +
  geom_line(aes(
    iteration, aum, color=set, linetype=search.type),
    size=1.3,
    data=loss.dt %>%
      filter((search.type=="exact" & (factor==1000)) |
               (search.type=="grid" & factor==100))
  ) +
  geom_point(aes(iteration, aum), color="darkblue", size=3, data=min.dt %>%
               filter((search.type=="exact" & (factor==1000)) |
                        (search.type=="grid" & factor==100))
  ) +
  scale_fill_brewer(type = "seq", palette = "Pastel1") +
  scale_y_log10() +
  xlim(0, 1000) +
  scale_x_log10()
