source("packages.R")
model.ord <- c("best","good","ok","bad")

d <- function(min.log.lambda, fp, fn){
  data.table(min.log.lambda, fp, fn)
}
profile <- function(..., possible.fp, possible.fn, errors, labels){
  dt <- do.call(rbind, list(...))
  if(missing(possible.fp))possible.fp <- max(dt$fp)
  if(missing(possible.fn))possible.fn <- max(dt$fn)
  errors <- dt[, fp+fn]
  if(missing(labels))labels <- max(errors)
  dt[, data.table(
    min.log.lambda,
    max.log.lambda=c(min.log.lambda[-1], Inf),
    fp, fn, errors, possible.fp, possible.fn, labels)]
}
profile.list <- list(
  best=profile(
    d(-Inf, 0, 10),
    d(1, 0, 4),
    d(3,0,2),
    d(4, 0, 0),
    d(5,3,0),
    d(7, 7, 0),
    d(8, 10, 0)),
  good=profile(
    d(-Inf, 0, 10),
    d(1,1,8),
    d(2, 1, 4),
    d(3, 1, 1),
    d(4,3,1),
    d(5, 8, 1),
    d(8, 10, 0)),
  ok=profile(
    d(-Inf, 0, 10),
    d(1, 0, 8),
    d(2, 4, 5),
    d(3, 4, 3),
    d(5, 7, 3),
    d(7, 7, 0),
    d(8, 10, 0)),
  bad=profile(
    d(-Inf, 0, 10),
    d(3,3,7),
    d(4,4,6),
    d(5,5,5),
    d(8, 10, 0)))  
pred.dt <- data.table(problem=1, pred.log.lambda=0)
roc.dt.list <- list()
auc.dt.list <- list()
for(profile.i in seq_along(profile.list)){
  p <- data.table(profile.list[[profile.i]], problem=1)
  roc.list <- penaltyLearning::ROChange(p, pred.dt, problem.vars="problem")
  model <- names(profile.list)[[profile.i]]
  roc.dt.list[[profile.i]] <- data.table(model, roc.list$roc)
  auc.dt.list[[profile.i]] <- with(roc.list, data.table(model, auc, aum))
}
roc.dt <- do.call(rbind, roc.dt.list)[
, aum := min.fp.fn*(max.thresh-min.thresh)
][
, AUM := sum(ifelse(is.finite(aum), aum, 0)), by=model
][
, `:=`(FP=fp, FN=fn, `Min(FP,FN)`=min.fp.fn, Model = factor(model, model.ord))
][
, FNR := 1-TPR
][
, `Min(FPR,FNR)` := pmin(FPR,FNR)
][
, letter := c(LETTERS,letters)[1:.N]
##, by=model
][]

auc.dt <- do.call(rbind, auc.dt.list)
fp.fn.dt <- data.table::melt(roc.dt, measure.vars=c("FP", "FN", "Min(FP,FN)"))
err.sizes <- c(
  FP=3,
  FN=2,
  "Min(FP,FN)"=1)
err.colors <- c(
  FP="red",
  FN="deepskyblue",
  "Min(FP,FN)"="black")
rate.names <- function(x)structure(x, names=gsub(
  "(?<=[PN])", "R", names(x), perl=TRUE))
rate.sizes <- rate.names(err.sizes)
rate.colors <- rate.names(err.colors)
poly.dt <- roc.dt[, {
  right <- .SD[-.N]
  left <- .SD[-1]
  zero <- rep(0, nrow(left))
  i <- 1:nrow(left)
  m <- left$model
  area <- ifelse(left$FPR < right$FPR, "negative", "positive")
  data.table(
    FPR=c(left$FPR, right$FPR, right$FPR, left$FPR),
    TPR=c(left$TPR, right$TPR, zero,  zero),
    area=rep(area, 4),
    seg=rep(i, 4))
}, by=model]
ggplot()+
  theme_bw()+
  scale_fill_manual(values=c(positive="black", negative="red"))+
  geom_polygon(aes(
    FPR, TPR, group=paste(seg, model), fill=area),
    alpha=0.2,
    data=poly.dt)+
  geom_path(aes(
    FPR, TPR),
    data=roc.dt)+
  geom_point(aes(
    FPR, TPR),
    data=roc.dt)+
  facet_grid(. ~ model, labeller=label_both)+
  coord_equal()+
  geom_text(aes(
    0.5, 0.5, label=sprintf("auc=%.2f", auc)),
    data=auc.dt)

auc.dt[, Model := factor(model, model.ord)]
auc.dt[, AUC := round(auc, 2)]
roc.join <- auc.dt[roc.dt, on="model"]
poly.join <- auc.dt[poly.dt, on="model"]
ggplot()+
  theme_bw()+
  scale_fill_manual(values=c(positive="black", negative="red"))+
  geom_polygon(aes(
    FPR, TPR, group=paste(seg, model), fill=area),
    alpha=0.2,
    data=poly.join)+
  geom_path(aes(
    FPR, TPR),
    data=roc.join)+
  facet_grid(. ~ AUC + model, labeller=label_both)+
  coord_equal()+
  scale_x_continuous(
    "False Positive Rate",
    breaks = seq(0, 1, by=0.5))+
  scale_y_continuous(
    "True Positive Rate",
    breaks = seq(0, 1, by=0.5))

some <- function(dt)dt[model %in% model.ord]
gg <- ggplot()+
  theme_bw()+
  scale_fill_manual(values=c(positive="black", negative="red"))+
  geom_polygon(aes(
    FPR, TPR, group=paste(seg, model)),
    alpha=0.2,
    data=some(poly.join))+
  geom_path(aes(
    FPR, TPR),
    data=some(roc.join))+
  facet_grid(.~Model+AUC, labeller=label_both)+
  coord_equal()+
  scale_x_continuous(
    "False Positive Rate",
    labels=c("0","0.5","1"),
    breaks = seq(0, 1, by=0.5))+
  scale_y_continuous(
    "True Positive Rate",
    breaks = seq(0, 1, by=0.5),
    labels=c("0","0.5","1"))
png("figure-more-than-one-new-binary.png", width=5, height=2, units="in", res=200)
print(gg)
dev.off()

leg <- "Error type"
gg <- ggplot()+
  facet_grid(. ~ Model + AUM, labeller=label_both)+
  theme_bw()+
  geom_rect(aes(
    xmin=min.thresh, xmax=max.thresh,
    ymin=0, ymax=value),
    color="grey",
    fill="grey",
    data=some(fp.fn.dt[variable=="Min(FP,FN)"]))+
  geom_segment(aes(
    min.thresh, value,
    color=variable, size=variable,
    xend=max.thresh, yend=value),
    data=some(fp.fn.dt))+
  scale_color_manual(leg, values=err.colors)+
  scale_size_manual(leg, values=err.sizes)+
  scale_x_continuous(
    "Constant added to predictions",
    breaks=seq(0, 10, by=2))+
  scale_y_continuous("Label errors", breaks=seq(0, 10, by=2))
png(
  "figure-more-than-one-new-binary-aum.png", 
  width=6, height=2, units="in", res=200)
print(gg)
dev.off()

roc.rate <- data.table(roc.dt)[
, aum := `Min(FPR,FNR)`*(max.thresh-min.thresh)
][
, AUM := sum(ifelse(is.finite(aum), aum, 0)), by=model
]
fpr.fnr.dt <- data.table::melt(
  roc.rate,
  measure.vars=c("FPR", "FNR", "Min(FPR,FNR)"))
ar <- function(model, x, y, xend, yend){
  data.table(Model=factor(model, model.ord), x, y, xend, yend)
}
arrows.only <- rbind(
  ar("best", 4.5, 0.8, 4.5, 0),
  ar("good", 6, 0.5, 6, 0.05),
  ar("ok", 4, 0.9, 4, 0.1),
  ar("bad", 6, 0.76, 6, 0.1))
arrows.auc <- fpr.fnr.dt[arrows.only, on="Model", mult="first"]
arrow.color <- "grey30"
gg <- ggplot()+
  facet_grid(. ~ Model, labeller=label_both)+
  theme_bw()+
  geom_vline(aes(
    xintercept=min.thresh),
    color="grey50",
    data=some(fpr.fnr.dt))+
  geom_rect(aes(
    xmin=min.thresh, xmax=max.thresh,
    ymin=0, ymax=value),
    color="grey",
    fill="grey",
    data=some(fpr.fnr.dt[variable=="Min(FPR,FNR)"]))+
  geom_segment(aes(
    min.thresh, value,
    color=variable, size=variable,
    xend=max.thresh, yend=value),
    data=some(fpr.fnr.dt))+
  geom_text(aes(
    ifelse(
      min.thresh==-Inf, max.thresh-1, ifelse(
        max.thresh==Inf, min.thresh+1, (min.thresh+max.thresh)/2)),
    -0.2,
    label=letter),
    vjust=0,
    size=3,
    data=some(fpr.fnr.dt))+
  geom_segment(aes(
    x, y, xend=xend, yend=yend),
    data=arrows.auc,
    color=arrow.color,
    size=1,
    arrow=grid::arrow(length=grid::unit(0.5,"lines"), type="closed"))+
  geom_label(aes(
    x, y,
    label=sprintf("AUM=%.1f", AUM)),
    vjust=0,
    color=arrow.color,
    size=3,
    data=arrows.auc)+
  scale_color_manual(leg, values=rate.colors)+
  scale_size_manual(leg, values=rate.sizes)+
  scale_x_continuous(
    "Constant added to predictions",
    limits=c(0,9),
    breaks=seq(0, 10, by=2))+
  scale_y_continuous("Label error rate", breaks=seq(0, 1, by=0.5))
print(gg)
png(
  "figure-more-than-one-new-binary-aum-rate.png", 
  width=6.6, height=2, units="in", res=200)
print(gg)
dev.off()

roc.join[, min.FPR.FNR := pmin(FPR,1-TPR)]
roc.join[, `sum(min)` := sum(min.FPR.FNR), by=Model]
rate.grid <- seq(0, 1, by=0.1)
grid.dt <- data.table(expand.grid(
  FPR=rate.grid, TPR=rate.grid
))[
, FNR := 1-TPR
][
, `Min(FPR,FNR)` := pmin(FPR,FNR)
][]
text.off <- 0.03
roc.lim <- c(-0.1,1.1)
gg <- ggplot()+
  theme_bw()+
  scale_fill_gradient(
    "Min(FPR,FNR)",
    low="white",
    high="purple")+
  geom_tile(aes(
    FPR, TPR, fill=`Min(FPR,FNR)`),
    data=grid.dt)+
  geom_path(aes(
    FPR, TPR),
    data=some(roc.join))+
  ## ggrepel::geom_text_repel(aes(
  ##   FPR, TPR, label=min.FPR.FNR),
  ##   size=3,
  ##   data=some(roc.join))+
  geom_text(aes(
    FPR-text.off,
    TPR+text.off,
    label=letter),
    hjust=1, vjust=0,
    size=2.5,
    data=some(roc.join))+
  geom_text(aes(
    FPR+text.off,
    TPR-text.off,
    label=min.FPR.FNR),
    hjust=0, vjust=1,
    size=2.5,
    data=some(roc.join))+
  ## geom_text(aes(
  ##   FPR+text.off,
  ##   TPR-text.off,
  ##   label=ifelse(
  ##     FPR==1,
  ##     paste0(min.FPR.FNR, "\n", letter),
  ##     paste(min.FPR.FNR, letter))),
  ##   hjust=0, vjust=1,
  ##   size=2.5,
  ##   data=some(roc.join))+
  geom_point(aes(
    FPR, TPR),
    data=some(roc.join))+
  facet_grid(.~`sum(min)`+AUC, labeller=label_both)+
  coord_equal()+
  scale_x_continuous(
    "False Positive Rate (FPR)",
    labels=c("0","0.5","1"),
    limits=roc.lim,
    breaks = seq(0, 1, by=0.5))+
  scale_y_continuous(
    "True Positive Rate\n(TPR = 1-FNR)",
    breaks = seq(0, 1, by=0.5),
    limits=roc.lim,
    labels=c("0","0.5","1"))
print(gg)
png(
  "figure-more-than-one-new-binary-heat.png", 
  width=6.6, height=2.2, units="in", res=200)
print(gg)
dev.off()

