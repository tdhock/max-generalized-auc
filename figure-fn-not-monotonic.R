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

err.list$model.errors[, diff.fp := c(diff(fp), NA)]
err.list$model.errors[, diff.fn := c(diff(fn), NA)]
diff.dt <- err.list$model.errors[diff.fp!=0 | diff.fn!=0, .(
  pred.log.lambda=max.log.lambda,
  diff.fp,
  diff.fn
)]
err.tall <- data.table::melt(err.list$model.errors, measure=c("fp", "fn"))
err.sizes <- c(
  fp=3,
  fn=2)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  correct="grey50")
err.colors[c("false positive", "false negative")] <- err.colors[c("fp", "fn")]
lab <- "Error type"
some.segs <- c(2:4, 15)
vline.dt <- err.list$model.errors[segments %in% some.segs]
vline.dt[, x := -(min.log.lambda+max.log.lambda)/2]
gg <- ggplot()+
  facet_grid(panel ~ ., scales="free")+
  theme_bw()+
  ##theme(panel.spacing=grid::unit(0, "lines"))+
  ylab("")+
  geom_vline(aes(
    xintercept=x),
    color="grey50",
    data=vline.dt)+
  geom_segment(aes(
    -min.log.lambda, value,
    color=variable, size=variable,
    xend=-max.log.lambda, yend=value),
    data=data.table(panel="Label errors", err.tall))+
  geom_blank(aes(
    x, y),
    data=data.table(x=0, y=c(-0.4,1.4)))+
  geom_segment(aes(
    -min.log.lambda, segments,
    xend=-max.log.lambda, yend=segments),
    size=1,
    data=data.table(panel="Segments", err.list$model.errors))+
  scale_color_manual(lab, values=err.colors)+
  scale_size_manual(lab, values=err.sizes)+
  scale_x_continuous(
    "Predicted value, f(x) = -log(penalty)",
    limits=c(-2, 4))
png("figure-fn-not-monotonic-error.png", 3.3, 2, units="in", res=200)
print(gg)
dev.off()

model.color <- "violet"
standAlone <- TRUE
suffix <- if(standAlone)"-standAlone" else ""
no.ext <- paste0("figure-fn-not-monotonic-error", suffix)
f.tex <- paste0(no.ext, ".tex")
tikz(f.tex, width=3, height=3, standAlone = standAlone)
left.lines <- 4
other.lines <- 1
ax.label.offset <- 1.5
par(mar=c(0,left.lines,other.lines,other.lines), cex.axis=1.5)
layout(cbind(c(rep(1,1),rep(2,2))))
xrange <- c(-2, 4)
xsize <- diff(xrange)
diff.dt[, x := c(0, -1, -1.5, -2)]
diff.dt[, y := -.I/2]
diff.dt[, pos := c(4,4,4,4)]
plot(xrange, range(err.tall[["segments"]]), type="n", yaxt="n", xaxt="n",ylab="",xlab="")
vline.dt[, abline(v=x, col=model.color)]
vline.dt[, text(x, 19, segments, pos=2, offset=0.2, col=model.color)]
err.list$model.errors[, mysegs(-max.log.lambda, -min.log.lambda, segments, lwd=3)]
axis(2,c(1,10,20),las=1)
mtext("Segments", 2, left.lines-ax.label.offset)
bottom.lines <- 4
par(mar=c(bottom.lines,left.lines,0,other.lines))
plot(xrange, c(min(diff.dt[["y"]]), 1.4), type="n", yaxt="n", xaxt="n",xlab="",ylab="")
axis(2,c(0,1),las=1)
mysegs <- function(x0, x1, y, ...)segments(
  ifelse(x0==-Inf, xrange[1]-xsize, x0), y,
  ifelse(x1==Inf, xrange[2]+xsize, x1), y,
  lend=1,
  ...)
err.tall[, mysegs(
  -max.log.lambda, -min.log.lambda, value,
  lwd=err.sizes[paste(variable)]*4,
  col=err.colors[paste(variable)])]
diff.dt[, text(x, y, sprintf(
  "$(v=%.3f,\\Delta\\text{FP}=%d,\\Delta\\text{FN}=%d)$",
  -pred.log.lambda, diff.fp, diff.fn),
  cex=1, pos=pos, offset=-0.5)]
leg.dt <- data.table(
  variable=c("fp","fn"),
  x=c(2.5,-0.5))
vline.dt[, segments(x, 0, x, 2, col=model.color)]
leg.dt[, text(x, 0.9, toupper(variable), col=err.colors[paste(variable)], cex=1.5)]
diff.dt[, segments(-pred.log.lambda, y+0.3, -pred.log.lambda, 0)]
mtext("Label errors", 2, left.lines-ax.label.offset)
axis(1)
mtext("Predicted value, $f(\\mathbf x_i) = -\\log \\hat \\lambda_i$", 1,bottom.lines-ax.label.offset)
dev.off()
if(standAlone)system(paste("pdflatex", no.ext))

some.segs.dt <- data.table(segments=some.segs)
show.labels <- err.list$label.errors[some.segs.dt, on="segments"]
show.segs <- segs.dt[show.labels, on="segments"]
show.change <- change.dt[show.labels, on="segments"]
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
    size=3,
    data=show.labels)+
  theme(legend.position="none")
png("figure-fn-not-monotonic.png", 5, 3.8, units="in", res=200)
print(gg)
dev.off()
