source("packages.R")

data(neuroblastoma, package="neuroblastoma")
nb.dt <- data.table(neuroblastoma$profiles)
p4.2 <- nb.dt[profile.id==4 & chromosome==2]
p4.2[, i := 1:.N]
label.dt <- data.table(problem=1, min=20, max=80, annotation="1breakpoint", label="one change")
max.segments <- 20
fit <- jointseg::Fpsn(p4.2$logratio, max.segments)
segs.dt.list <- list()
models.dt.list <- list()
for(segments in 1:max.segments){
  end <- fit$t.est[segments, 1:segments]
  start <- c(0L, end[-length(end)])+1L
  start.end <- data.table(start, end)
  segs.dt.list[[paste(segments)]] <- data.table(
    segments,
    p4.2[start.end, .(start, end, mean=mean(logratio)), on=.(i >= start, i <= end), by=.EACHI])
  models.dt.list[[paste(segments)]] <- data.table(
    segments,
    loss=fit$J.est[segments])
}
##models.dt <- do.call(rbind, models.dt.list)
models.dt <- data.table(problem=1, segments=1:max.segments, loss=fit$J.est)
segs.dt <- do.call(rbind, segs.dt.list)[, .(problem=1, segments, start, end, mean)]
change.dt <- segs.dt[1 < start][, change := start-0.5][]
selected.dt <- penaltyLearning::modelSelection(models.dt, "loss", "segments")
err.list <- penaltyLearning::labelError(
  selected.dt, label.dt, change.dt,
  problem.vars="problem", change.var="change", model.vars="segments")

err.tall <- data.table::melt(err.list$model.errors, measure=c("fp", "fn"))
err.sizes <- c(
  fp=3,
  fn=2)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  correct="black")
err.colors[c("false positive", "false negative")] <- err.colors[c("fp", "fn")]
lab <- "Error type"
some.segs <- c(2:4, 15)
gg <- ggplot()+
  facet_grid(panel ~ ., scales="free")+
  theme_bw()+
  ##theme(panel.spacing=grid::unit(0, "lines"))+
  ylab("")+
  geom_vline(aes(
    xintercept=(min.log.lambda+max.log.lambda)/2),
    color="grey50",
    data=err.list$model.errors[segments %in% some.segs])+
  geom_segment(aes(
    min.log.lambda, value,
    color=variable, size=variable,
    xend=max.log.lambda, yend=value),
    data=data.table(panel="Label errors", err.tall))+
  geom_blank(aes(
    x, y),
    data=data.table(x=0, y=c(-0.4,1.4)))+
  geom_segment(aes(
    min.log.lambda, segments,
    xend=max.log.lambda, yend=segments),
    size=1,
    data=data.table(panel="Segments", err.list$model.errors))+
  scale_color_manual(lab, values=err.colors)+
  scale_size_manual(lab, values=err.sizes)+
  scale_x_continuous("Predicted log(penalty)")
png("figure-fn-not-monotonic-error.png", 5, 3, units="in", res=200)
print(gg)
dev.off()

some.segs.dt <- data.table(segments=some.segs)
show.labels <- err.list$label.errors[some.segs.dt, on="segments"]
show.segs <- segs.dt[show.labels, on="segments"]
show.change <- change.dt[show.labels, on="segments"]
model.color <- "violet"
gg <- ggplot()+
  theme_bw()+
  geom_rect(aes(
    xmin=min, xmax=max,
    ymin=-Inf, ymax=Inf,
    fill=label),
    alpha=0.5,
    data=label.dt)+
  geom_rect(aes(
    xmin=min, xmax=max,
    ymin=-Inf, ymax=Inf,
    color=status),
    fill=NA,
    size=2,
    data=show.labels)+
  scale_color_manual(values=err.colors)+
  ## geom_rect(aes(
  ##   xmin=min, xmax=max,
  ##   ymin=-Inf, ymax=Inf,
  ##   linetype=status),
  ##   fill=NA,
  ##   color="black",
  ##   size=2,
  ##   data=show.labels)+
  scale_fill_manual(values=c("one change"="grey50"))+
  scale_linetype_manual(
    "Error type",
    values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))+
  geom_point(aes(
    i, logratio),
    data=p4.2)+
  geom_segment(aes(
    start-0.5, mean,
    xend=end+0.5, yend=mean),
    color=model.color,
    size=1,
    data=show.segs)+
  geom_vline(aes(
    xintercept=change),
    data=show.change,
    size=1,
    color=model.color)+
  facet_grid(segments ~ ., labeller=label_both)+
  scale_x_continuous("Data sequence index")+
  scale_y_continuous("Data value")+
  geom_label(aes(
    (min+max)/2, -0.5,
    label=status,
    color=status),
    data=show.labels)+
  theme(legend.position="none")
png("figure-fn-not-monotonic.png", 6, 4, units="in", res=200)
print(gg)
dev.off()
