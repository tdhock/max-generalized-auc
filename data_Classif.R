library(ggplot2)
library(data.table)
(csv.vec <- Sys.glob(file.path("data_Classif","*.csv")))
result.list <- list()
for(csv.i in seq_along(csv.vec)){
  classif.csv <- csv.vec[[csv.i]]
  cat(sprintf("%d / %d %s\n", csv.i, length(csv.vec), classif.csv))
  data.name <- sub(".csv","",basename(classif.csv))
  raw.dt <- fread(classif.csv)
  set.seed(1)
  raw.y.tab <- table(raw.dt$y)
  neg.y <- if("0" %in% names(raw.y.tab)) 0 else raw.dt$y[1]
  classif.dt <- raw.dt[
    sample(.N)
  ][
  , y := ifelse(y==neg.y, 0, 1)
  ][]
  group.name <- names(classif.dt)[1]
  feature.names <- setdiff(names(classif.dt),c(group.name,"y"))
  y.tab <- table(classif.dt[["y"]])
  max.N <- nrow(classif.dt)
  min.N <- 10*max.N/min(y.tab)
  N.vec <- as.integer(10^seq(log10(min.N),log10(max.N),by=0.25))
  expr.list <- atime::atime_grid(list(
    maxIt=c("min.aum","linear","quadratic"),
    seed=1:4
  ), lm={
    set.seed(seed)
    fit <- aum:::aum_linear_model(
      feature.list, diff.list,
      maxIterations=maxit.list[[maxIt]],
      improvement.thresh=imp.thr)
    data.table(
      fit=list(fit),
      sum.iterations=sum(fit$search$iterations),
      steps=nrow(fit$search),
      mean.it.per.step=mean(fit$search$iterations),
      max.intersections,
      min.valid.AUM=fit$loss[set=="validation", min(aum)])
  })
  result.list[[data.name]] <- atime::atime(
    N=N.vec,
    setup={
      N.class <- round(y.tab/sum(y.tab)*N)
      N.dt <- classif.dt[, .SD[1:N.class[[paste(y)]]], by=y]
      y.vec <- as.integer(factor(N.dt[["y"]]))-1L
      fmat <- as.matrix(N.dt[,feature.names,with=FALSE])
      set.names <- c("subtrain","validation")
      set.vec <- rep(set.names, l=nrow(fmat))
      table(set.vec,y.vec)
      X.sc <- scale(fmat)
      keep <- attr(X.sc,"scaled:scale")!=0
      X.keep <- X.sc[,keep]
      feature.list <- list()
      diff.list <- list()
      out.dt <- data.table(set=set.vec, y=y.vec, X.keep)
      out.csv <- sprintf("data_Classif_scaled/%s_N=%d.csv", data.name, N)
      out.dir <- dirname(out.csv)
      dir.create(out.dir)
      fwrite(out.dt, out.csv)
      for(set.name in set.names){
        is.set <- set.vec==set.name
        feature.list[[set.name]] <- X.keep[is.set,]
        diff.list[[set.name]] <- aum::aum_diffs_binary(
          y.vec[is.set], denominator="rate")
      }
      imp.thr <- if(FALSE){
        diff.list[[1]][, .(
          diff=c(fp_diff,fn_diff)
        )][diff!=0, min(abs(diff))/100]
        print(imp.thr)
      }else 1e-3
      N.subtrain <- nrow(feature.list$subtrain)
      max.intersections <- N.subtrain*(N.subtrain-1)/2
      maxit.list <- list(
        min.aum="min.aum",
        ##max.auc="max.auc",#??
        linear=N.subtrain,
        quadratic=max.intersections)
    },
    verbose=TRUE,
    result=TRUE,
    times=1,
    seconds.limit=10,
    expr.list=expr.list)
}
if(FALSE){
  plot(result.list$FishSonar)
  plot(ref.list$N)
}
save(result.list, file="data_Classif.RData")

