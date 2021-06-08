source("packages.R")

N.vec <- as.integer(10^seq(1, 6, by=0.5))
max.N <- max(N.vec)
all.labels.vec <- rep(c(-1,1), l=max.N)
all.diffs.dt <- aum::aum_diffs_binary(all.labels.vec)
set.seed(1)
all.pred.vec <- rnorm(max.N)
timing.dt.list <- list()
for(N in N.vec){
  print(N)
  N.pred.vec <- all.pred.vec[1:N]
  N.diffs.dt <- all.diffs.dt[1:N]
  N.labels.vec <- sort(all.labels.vec[1:N])
  order.list <- list(sorted=sort(N.pred.vec), unsorted=N.pred.vec)
  for(prediction.order in names(order.list)){
    order.pred.vec <- order.list[[prediction.order]]
    timing.df <- microbenchmark::microbenchmark(logistic.grad={
      aum:::logistic_grad(order.pred.vec, N.labels.vec)
    }, sort={
      aum:::do_qsort(order.pred.vec)
    }, aum={
      aum::aum(N.diffs.dt, order.pred.vec)
    }, times=10)
    timing.dt.list[[paste(N, prediction.order)]] <- with(timing.df, data.table(
      N, prediction.order, seconds=time/1e9, algorithm=expr))
  }
}
(timing.dt <- do.call(rbind, timing.dt.list))

timing.stats <- timing.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(N, prediction.order, algorithm)]
gg <- ggplot()+
  facet_grid(. ~ prediction.order, labeller=label_both)+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=algorithm),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N, median, color=algorithm),
    data=timing.stats)+
  scale_x_log10(
    "Number of predictions",
    limits=c(10, max.N*10))+
  scale_y_log10(
    "Computation time in seconds,
median line, min/max band, 10 timings")
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-aum-grad-speed-binary-cpp-algos.png", width=10, height=3, res=200, units="in")
print(dl)
dev.off()

gg <- ggplot()+
  facet_grid(. ~ algorithm, labeller=label_both)+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=prediction.order),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N, median, color=prediction.order),
    data=timing.stats)+
  scale_x_log10(
    "Number of predictions",
    limits=c(10, max.N*10))+
  scale_y_log10(
    "Computation time in seconds,
median line, min/max band, 10 timings")
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-aum-grad-speed-binary-cpp.png", width=10, height=3, res=200, units="in")
print(dl)
dev.off()
