library(data.table)
library(ggplot2)
error.dt <- data.table::fread(
  "../feature-learning-benchmark/labeled_problems_errors.csv")
error.dt[, min.lambda := exp(min.log.penalty)]
levs <- unique(error.dt$prob.dir)
error.dt[, efac := factor(prob.dir, levs)]
error.dt[, example := as.integer(efac)-1L]
error.dt[, run := cumsum(c(1, diff(errors)!=0)), by=example]
diff.dt <- aum::aum_diffs_penalty(error.dt)
run.dt <- error.dt[, .(
  max=max(max.log.penalty),
  min=min(min.log.penalty)
), by=.(example,errors,run)]
min.dt <- run.dt[, .SD[errors==min(errors)], by=example]
big.dt <- min.dt[, .SD[size==max(size)][1], by=example]
bad <- big.dt[, .(count=.N), by=example][count>1]
run.dt[bad, on="example"]
some.ids <- big.dt[order(size)][c(1,10,100,500,1000), example]
some.vline <- big.dt[example %in% some.ids]
some.diffs <- diff.dt[example %in% some.ids]
some.err <- aum::aum_errors(some.diffs)
some.tall <- melt(
  some.err,
  measure.vars = c("fp","fn"),
  variable.name="error_type",
  value.name="label_errors")
err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
ggplot()+
  geom_vline(aes(
    xintercept=log.penalty),
    data=some.vline)+
  geom_segment(aes(
    min.pred, label_errors,
    color=error_type, size=error_type,
    xend=max.pred, yend=label_errors),
    data=some.tall)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)+
  facet_grid(example ~ ., scales="free")

diff.tall <- nc::capture_melt_single(
  diff.dt,
  error_type="fp|fn",
  "_diff",
  value.name="diff")[diff != 0]
diff.tall[, `:=`(
  next.diff = c(diff[-1],NA), 
  next.pred = c(pred[-1],NA)
), by=.(example, error_type)]
small.pred.diff <- diff.tall[
  !is.na(next.diff) & next.diff>0 & diff<0 & next.pred-pred<1]
small.pred.diff[-diff < next.diff]
small.pred.diff[-diff > next.diff]
small.pred.diff[error_type=="fn"][order(next.diff+diff)]
small.pred.diff[error_type=="fn"][order(next.diff-diff)]
some.ids <- c(291,778,4808,4483)
some.diffs <- diff.dt[example %in% some.ids]
some.err <- aum::aum_errors(some.diffs)

some.tall <- melt(
  some.err,
  measure.vars = c("fp","fn"),
  variable.name="error_type",
  value.name="label_errors")
err.sizes <- c(
  fp=3,
  fn=2,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
ggplot()+
  geom_segment(aes(
    min.pred, label_errors,
    color=error_type, size=error_type,
    xend=max.pred, yend=label_errors),
    data=some.tall)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)+
  facet_grid(example ~ ., scales="free")
## adaptive margin size?
