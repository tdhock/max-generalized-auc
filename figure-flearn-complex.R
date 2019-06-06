source("packages.R")

data.list <- list()
for(data.type in c("possible_errors", "errors")){
  f <- paste0(
    "../feature-learning-benchmark/labeled_problems_", data.type, ".csv")
  data.list[[data.type]] <- fread(f)
}


prob.dir.vec <- data.list$errors[min.log.penalty==-Inf & 0<fn, prob.dir]
some.err <- data.list$errors[prob.dir.vec, on=list(prob.dir)]
some.err[, list(min.fn=min(fn)), by=list(prob.dir)]
err.sizes <- c(
  fp=2,
  fn=1.5,
  errors=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  errors="black")
some.err.tall <- melt(
  some.err,
  measure.vars=names(err.colors))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(prob.dir ~ ., scales="free")+
  geom_segment(aes(
    min.log.penalty, value,
    xend=max.log.penalty, yend=value,
    color=variable, size=variable),
    data=some.err.tall)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)+
  scale_y_continuous(limits=c(-2, NA),breaks=seq(0, 20, by=4))
