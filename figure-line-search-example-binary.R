source("packages.R")
label.vec <- c(0,0,1,1)
pred.vec <- c(-8,-2,0,-4)
(bin.diffs <- aum::aum_diffs_binary(label.vec))
bin.line.search <- aum::aum_line_search(bin.diffs, pred.vec=pred.vec)
if(requireNamespace("ggplot2"))plot(bin.line.search)

roc.df <- WeightedROC::WeightedROC(pred.vec, label.vec)
ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=roc.df)+
  geom_point(aes(
    FPR, TPR),
    data=roc.df)

step.dt <- data.table(step.size=seq(0,2))
pred.dt <- step.dt[, .(
  pred=pred.vec-step.size*bin.line.search$gradient
), by=step.size]
roc.dt <- pred.dt[
, WeightedROC::WeightedROC(pred, label.vec), by=step.size
][
  order(step.size,-threshold)
][, `:=`(
  FNR = 1-TPR,
  const = -threshold,
  letter = LETTERS[1:.N]
)][
, `min(FPR,FNR)` := pmin(FPR,FNR)
][
, max.const := c(const[-1],Inf)
, by=step.size
][
, letter.const := ifelse(const== -Inf, max.const-1, ifelse(
  max.const==Inf, const+1, (const+max.const)/2))
][]
roc.segs <- roc.dt[, .(
  next.FPR=FPR[-1], FPR=FPR[-.N],
  next.TPR=TPR[-1], TPR=TPR[-.N],
  seg.i=seq(1,.N-1)
), by=step.size
][
, updated := ifelse(seg.i %in% c(1,.N), "no","yes")
, by=step.size][]
poly.dt <- roc.segs[updated=="yes", .(
  FPR=c(FPR,FPR,next.FPR,next.FPR),
  TPR=c(0,TPR,next.TPR,0)
), by=.(step.size,seg.i)]
gg <- ggplot()+
  theme_bw()+
  geom_polygon(aes(
    FPR, TPR, group=seg.i),
    fill="grey",
    data=poly.dt)+
  ## geom_path(aes(
  ##   FPR, TPR),
  ##   data=roc.dt)+
  geom_segment(aes(
    FPR, TPR,
    xend=next.FPR, yend=next.TPR,
    color=updated),
    size=1,
    data=roc.segs)+
  scale_color_manual(
    values=c(no="black",yes="orange"))+
  geom_point(aes(
    FPR, TPR),
    data=roc.dt)+
  geom_text(aes(
    FPR, TPR, label=letter),
    hjust=1.2,
    vjust=-0.2,
    data=roc.dt)+
  coord_equal()+
  facet_grid(. ~ step.size, labeller=label_both)+
  scale_x_continuous(
    "False Positive Rate",
    limits=c(-0.2,1.2),
    breaks=seq(0,1,by=0.5))+
  scale_y_continuous(
    "True Positive Rate",
    limits=c(-0.2,1.2),
    breaks=seq(0,1,by=0.5))
png("figure-line-search-example-binary-roc.png", 4.5, 1.8, units="in", res=400)
print(gg)
dev.off()

err.sizes <- c(
  "min(FPR,FNR)"=1,
  FPR=3,
  FNR=2)
err.colors <- c(
  "min(FPR,FNR)"="black",
  FPR="red",
  FNR="deepskyblue")
roc.long <- melt(
  roc.dt,
  measure=c("FNR","FPR","min(FPR,FNR)")
)
seg <- function(x,xend,y,step.size){
  data.table(x,xend,y,step.size)
}
seg.dt <- rbind(
  seg(2,3,0.9,0),
  seg(3,4,0.9,1),
  seg(4,5,0.9,2),
  seg(4,3,0.8,0),
  seg(3,2,0.8,1),
  seg(2,1,0.8,2))
gg <- ggplot()+
  theme_bw()+
  scale_size_manual(
    "Error type",
    values=err.sizes)+
  scale_color_manual(
    "Error type",
    values=err.colors)+
  geom_rect(aes(
    xmin=const, xmax=max.const,
    ymin=0, ymax=`min(FPR,FNR)`),
    fill="grey80",
    data=roc.dt)+
  geom_vline(aes(
    xintercept=const),
    color="grey",
    size=1,
    data=roc.dt[is.finite(const)])+
  geom_segment(aes(
    const, value,
    xend=max.const, yend=value,
    color=variable,
    size=variable),
    data=roc.long)+
  facet_grid(. ~ step.size, labeller=label_both)+
  geom_text(aes(
    letter.const,
    0.25,
    label=letter),
    size=3,
    data=roc.dt)+
  scale_x_continuous(
    "Constant added to predicted values",
    breaks=seq(-10,10,by=2),
    limits=c(-1.4,9.4))+
  scale_y_continuous("")+
  geom_segment(aes(
    x, y,
    xend=xend, yend=y),
    size=0.5,
    color="grey50",
    arrow=grid::arrow(length=grid::unit(0.1,"cm"),type="open"),
    data=seg.dt)
png("figure-line-search-example-binary-error.png", 4.5, 1.5, units="in", res=400)
print(gg)
dev.off()

aum.df <- data.frame(panel="AUM", bin.line.search$line_search_result)
last <- bin.line.search$line_search_result[.N]
step.after <- last$step*1.05
step.diff <- step.after-last$step
last.seg <- data.table(
  panel="AUM",
  last,
  step.after, 
  aum.after=last[, step.diff*aum.slope.after+aum])
auc.segs <- bin.line.search$line_search_result[, data.table(
  panel="AUC", 
  min.step.size=step.size,
  max.step.size=c(step.size[-1],Inf),
  auc=auc.after)]
all.points <- suppressWarnings(melt(
  bin.line.search$line_search_result[, AUC := auc],
  measure.vars=c("AUC"),
  variable.name="panel"))
Constant <- "Constant added\nto pred. values"
abline.df <- data.frame(panel=Constant, bin.line.search$line_search_input)
gg <- ggplot2::ggplot()+
  ggplot2::theme_bw()+
  ggplot2::geom_vline(ggplot2::aes(
    xintercept=step.size),
    color="grey",
    size=1,
    data=step.dt)+
  ggplot2::geom_line(ggplot2::aes(
    step.size, aum),
    size=1,
    data=aum.df)+
  ggplot2::geom_segment(ggplot2::aes(
    min.step.size, auc,
    xend=max.step.size, yend=auc),
    size=1,
    data=auc.segs)+
  ggplot2::geom_segment(ggplot2::aes(
    step.size, aum,
    xend=step.after, yend=aum.after),
    linetype="dotted",
    size=1,
    data=last.seg)+
  ggplot2::geom_point(ggplot2::aes(
    step.size, value),
    shape=1,
    data=all.points)+
  ggplot2::facet_grid(panel ~ ., scales="free")+
  ggplot2::geom_abline(ggplot2::aes(
    slope=slope, intercept=intercept),
    size=1,
    data=abline.df)+
  geom_text(aes(
    step.size, letter.const, label=letter),
    size=3,
    data=data.table(panel=Constant,roc.dt))+
  geom_blank(aes(
    0, y),
    data=data.table(y=c(-1.5,11),panel=Constant))+
  ggplot2::scale_y_continuous("")+
  ggplot2::scale_x_continuous(
    "Step size",
    breaks=seq(0,10))
png("figure-line-search-example-binary.png", 3, 3.5, units="in", res=400)
print(gg)
dev.off()


if(FALSE){
  N <- 100
  set.seed(2)
  (bin.diffs <- aum::aum_diffs_binary(round(runif(N))))
  bin.line.search <- aum::aum_line_search(
    bin.diffs, pred.vec=rnorm(N), maxIterations=N*N)
  if(requireNamespace("ggplot2"))plot(bin.line.search)
  bin.line.search$line_search_result$aum.slope.after
  ## check slope not increasing => convex.
}
