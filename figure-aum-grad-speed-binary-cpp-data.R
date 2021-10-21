source("packages.R")

N.vec <- as.integer(10^seq(1, 6, by=0.5))
max.N <- max(N.vec)
all.labels.vec <- rep(c(-1,1), l=max.N)
all.diffs.dt <- aum::aum_diffs_binary(all.labels.vec)
set.seed(1)
all.pred.vec <- rnorm(max.N)
timing.dt.list <- list()

do.sub <- function(...){
  mcall <- match.call()
  L <- as.list(mcall[-1])
  for(arg.name in names(L)){
    maybe.lang <- L[[arg.name]]
    if(is.language(maybe.lang)){
      L[[arg.name]] <- substitute(
        result.list[[NAME]] <- EXPR,
        list(NAME=arg.name, EXPR=maybe.lang))
    }
  }
  L
}

for(N in N.vec){
  print(N)
  N.pred.vec <- all.pred.vec[1:N]
  N.diffs.dt <- all.diffs.dt[1:N]
  N.labels.vec <- sort(all.labels.vec[1:N])
  order.list <- list(sorted=sort(N.pred.vec), unsorted=N.pred.vec)
  for(prediction.order in names(order.list)){
    order.pred.vec <- order.list[[prediction.order]]
    result.list <- list()
    m.args <- c(do.sub(`Logistic`={
      aum:::logistic_grad(order.pred.vec, N.labels.vec)
    }, AUM={
      aum::aum(N.diffs.dt, order.pred.vec)
    }),
    if(N < 1e4)do.sub(`Squared Hinge All Pairs`={
      is.positive <- N.labels.vec == 1
      pairs.dt <- data.table(expand.grid(
        positive=which(is.positive),
        negative=which(!is.positive)))
      margin <- 1
      pairs.dt[, diff := order.pred.vec[positive]-order.pred.vec[negative]-margin]
      pairs.dt[, diff.clipped := ifelse(diff<0, diff, 0)]
      pairs.tall <- data.table::melt(
        pairs.dt,
        measure.vars=c("positive", "negative"),
        value.name="pred.i",
        variable.name="label")
      pairs.tall[, grad.sign := ifelse(label=="positive", 1, -1)]
      grad.dt <- pairs.tall[, .(
        gradient=sum(grad.sign*diff.clipped)
      ), keyby=pred.i]
      grad.dt[["gradient"]]
    }),
    times=10)
    timing.df <- do.call(microbenchmark::microbenchmark, m.args)
    timing.dt.list[[paste(N, prediction.order)]] <- with(timing.df, data.table(
      N, prediction.order, seconds=time/1e9, algorithm=expr))
  }
}
(timing.dt <- do.call(rbind, timing.dt.list))

data.table::fwrite(timing.dt, "figure-aum-grad-speed-binary-cpp-data.csv")

