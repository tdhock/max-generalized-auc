source("packages.R")

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
    d(1.5,0,2),
    d(2, 0, 0),
    d(2.5,3,0),
    d(3, 7, 0),
    d(4, 10, 0)),
  bad=profile(
    d(-Inf, 0, 10),
    d(1,1,9),
    d(2,2,8),
    d(3,3,7),
    d(4,4,6),
    d(5,5,5),
    d(6, 10, 0)),
  less=profile(
    d(-Inf, 0, 10),
    d(2, 1, 1),
    d(4, 10, 0)),
  good=profile(
    d(-Inf, 0, 10),
    d(0.5,1,8),
    d(1, 1, 4),
    d(2, 1, 1),
    d(2.5,3,1),
    d(3, 8, 1),
    d(4, 10, 0)),
  ok=profile(
    d(-Inf, 0, 10),
    d(2, 0, 8),
    d(5, 4, 5),
    d(7, 4, 3),
    d(8, 7, 3),
    d(9, 7, 0),
    d(10, 10, 0)),
  more=profile(
    d(-Inf, 0, 10),
    d(2, 8/3, 8/3),
    d(5, 10, 8/3),
    d(7, 10, 25/3),
    d(8, 5/3, 25/3),
    d(9, 5/3, 8/3),
    d(10, 10, 0)))
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
roc.dt <- do.call(rbind, roc.dt.list)
auc.dt <- do.call(rbind, auc.dt.list)
fp.fn.dt <- data.table::melt(roc.dt, measure.vars=c("fp", "fn", "min.fp.fn"))
err.sizes <- c(
  fp=3,
  fn=2,
  min.fp.fn=1)
err.colors <- c(
  fp="red",
  fn="deepskyblue",
  min.fp.fn="black")
ggplot()+
  facet_grid(model ~ ., labeller=label_both)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_rect(aes(
    xmin=min.thresh, xmax=max.thresh,
    ymin=0, ymax=value),
    color="grey",
    fill="grey",
    data=fp.fn.dt[variable=="min.fp.fn"])+
  geom_segment(aes(
    min.thresh, value,
    color=variable, size=variable,
    xend=max.thresh, yend=value),
    data=fp.fn.dt)+
  scale_color_manual(values=err.colors)+
  scale_size_manual(values=err.sizes)+
  scale_x_continuous("Constant added to predictions")

ggplot()+
  geom_path(aes(
    FPR, TPR, color=model, size=model, group=model),
    data=roc.dt)

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

model.ord <- c("best","good","ok","bad")
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
png("figure-more-than-one-binary.png", width=5, height=2, units="in", res=200)
print(gg)
dev.off()

roc.join[, min.FPR.FNR := pmin(FPR,1-TPR)]
roc.join[, `sum(min)` := sum(min.FPR.FNR), by=Model]
gg <- ggplot()+
  theme_bw()+
  scale_fill_gradient2(
    "min(FPR,FNR)",
    midpoint=0.25,
    low="blue",
    mid="white",
    high="red")+
  geom_path(aes(
    FPR, TPR),
    data=some(roc.join))+
  geom_point(aes(
    FPR, TPR, fill=min.FPR.FNR),
    shape=21,
    data=some(roc.join))+
  facet_grid(.~Model+`sum(min)`, labeller=label_both)+
  coord_equal()+
  scale_x_continuous(
    "False Positive Rate",
    labels=c("0","0.5","1"),
    breaks = seq(0, 1, by=0.5))+
  scale_y_continuous(
    "True Positive Rate",
    breaks = seq(0, 1, by=0.5),
    labels=c("0","0.5","1"))
print(gg)
png(
  "figure-more-than-one-binary-dots.png", 
  width=6, height=2, units="in", res=200)
print(gg)
dev.off()

err.sizes <- c(
  FP=3,
  FN=2,
  "min(FP,FN)"=1)
err.colors <- c(
  FP="red",
  FN="deepskyblue",
  "min(FP,FN)"="black")
binary.list <- list(
  "1"=profile(
    d(-Inf, 0, 1),
    d(0, 0, 0)),
  "0"=profile(
    d(-Inf, 0, 0),
    d(0, 1, 0)))
leg <- "Error type"
label.tall.dt.list <- list()
for(label in names(binary.list)){
  label.dt <- binary.list[[label]]
  label.tall <- melt(label.dt, measure=c("fp","fn"))
  label.tall[, Variable := toupper(variable)]
  label.tall.dt.list[[label]] <- data.table(label, label.tall)
  ggplot()+
    scale_color_manual(leg, values=err.colors)+
    scale_size_manual(leg, values=err.sizes)+
    scale_y_continuous("Label Errors")+
    scale_x_continuous(
      "Predicted score")+
    geom_segment(aes(
      min.log.lambda, value,
      color=Variable,
      size=Variable,
      xend=max.log.lambda, yend=value),
      data=label.tall)
}
label.tall.dt <- do.call(rbind, label.tall.dt.list)
lab.info <- rbind(
  data.table(label="0", hjust=1.1),
  data.table(label="1", hjust=-0.1))
label.non.zero <- label.tall.dt[value>0][lab.info, on="label"]
label.non.zero[, x := ifelse(label=="0", min.log.lambda, max.log.lambda)]
gg <- ggplot()+
  scale_color_manual(leg, values=err.colors)+
  scale_size_manual(leg, values=err.sizes)+
  scale_y_continuous(
    "Label Errors",
    breaks=0:1)+
  scale_x_continuous(
    "Predicted score f(x)",
    limits=c(-2,2))+
  geom_segment(aes(
    min.log.lambda, value,
    color=Variable,
    size=Variable,
    xend=max.log.lambda, yend=value),
    data=label.tall.dt)+
  geom_text(aes(
    x, 1, hjust=hjust, label=Variable, color=Variable),
    vjust=1, 
    data=label.non.zero)+
  facet_grid(label ~ ., labeller=label_both)
png(
  "figure-more-than-one-binary-errors.png",
  width=3, height=2.5, units="in", res=200)
print(gg)
dev.off()

for(m in names(profile.list)){
  p.roc <- roc.dt[model==m]
  p.auc <- auc.dt[model==m]
  p.poly <- poly.dt[model==m]
  p.fp.fn <- fp.fn.dt[model==m]
  p.fp.fn[, Variable := ifelse(
    variable=="min.fp.fn", "min(FP,FN)", toupper(paste(variable)))]
  p.rect <- p.fp.fn[Variable=="min(FP,FN)"]
  p.best <- p.roc[errors==min(errors)]
  best.color <- "green"
  p.roc[, q := .I]
  p.roc[, hjust := 0]
  p.roc[, vjust := 1.2]
  p.roc[q==6, `:=`(hjust=1, vjust=-0.4)]
  g <- ggplot()+
    theme_bw()+
    scale_fill_manual(values=c(positive="black", negative="red"))+
    geom_polygon(aes(
      FPR, TPR, group=paste(seg, model), fill=area),
      alpha=0.2,
      data=p.poly)+
    geom_path(aes(
      FPR, TPR),
      data=p.roc)+
    geom_text(aes(
      FPR+0.01, TPR, label=paste0("q=",q), hjust=hjust, vjust=vjust),
      size=3,
      data=p.roc)+
    geom_point(aes(
      FPR, TPR),
      data=p.roc)+
    coord_equal(xlim=c(0,1.1), ylim=c(-0.05, 1))+
    scale_y_continuous("True Positive Rate", breaks=c(0,0.5,1))+
    scale_x_continuous("False Positive Rate", breaks=c(0,0.5,1))+
    geom_text(aes(
      0.5, 0.5, label=sprintf("AUC=%.2f", auc)),
      data=p.auc)+
    theme(legend.position="none")
  ##if(all(p.poly$area=="positive"))g <- g+theme(legend.position="none")
  g.aum <- ggplot()+
    theme_bw()+
    theme(
      panel.grid.minor=element_blank(),
      panel.spacing=grid::unit(0, "lines"))+
    geom_rect(aes(
      xmin=min.thresh, xmax=max.thresh,
      ymin=0, ymax=value),
      color="grey",
      alpha=0.75,
      fill="grey",
      data=p.rect)+
    geom_text(aes(
      6, 4, label=sprintf("AUM=%.0f", aum)),
      data=p.auc)+
    geom_segment(aes(
      min.thresh, value,
      color=Variable, size=Variable,
      xend=max.thresh, yend=value),
      data=p.fp.fn)+
    geom_text(aes(
      ifelse(
        min.thresh == -Inf, max.thresh-1, ifelse(
          max.thresh == Inf, min.thresh+1, (min.thresh+max.thresh)/2)),
      -4, label=paste0("q=",q)),
      vjust=0,
      size=2.5,
      data=p.roc)+
    scale_color_manual(leg, values=err.colors)+
    scale_size_manual(leg, values=err.sizes)+
    scale_y_continuous("Label Errors")+
    scale_x_continuous(
      "Constant added to predicted values",
      limits=c(0, 12),
      breaks=p.roc[["min.thresh"]])
  g.list <- list(auc=g, aum=g.aum)
  for(plot.type in names(g.list)){
    out.png <- sprintf(
      "figure-more-than-one-%s-%s.png",
      m, plot.type)
    png(
      out.png,
      width=if(plot.type=="auc")2.5 else 4.5,
      height=if(plot.type=="auc")2.5 else 2, units="in", res=200)
    print(g.list[[plot.type]])
    dev.off()
    f.tex <- sub("png", "tex", out.png)
    tikz(f.tex, width=3, height=3, standAlone = TRUE)
    print(g.list[[plot.type]])
    dev.off()
    system(paste("pdflatex", f.tex))
  }
  g <- ggplot()+
    theme_bw()+
    scale_fill_manual(values=c(positive="black", negative="red"))+
    geom_polygon(aes(
      FPR, TPR, group=paste(seg, model), fill=area),
      alpha=0.2,
      data=p.poly)+
    geom_path(aes(
      FPR, TPR),
      data=p.roc)+
    geom_text(aes(
      FPR+0.01, TPR, label=sprintf("q=%d",q), hjust=hjust, vjust=vjust),
      size=3,
      data=p.roc)+
    geom_point(aes(
      FPR, TPR),
      size=0.5,
      data=p.roc)+
    coord_equal(
      xlim=c(0,1.1), ylim=c(-0.05, 1))+
    scale_y_continuous("True Positive Rate", breaks=c(0,0.5,1))+
    scale_x_continuous("False Positive Rate", breaks=c(0,0.5,1))+
    geom_text(aes(
      0.5, 0.5, label=sprintf("AUC=%.2f", auc)),
      data=p.auc)+
    theme(legend.position="none")
  limits.vec <- c(0, 12)
  type.breaks <- c(
    FP="$\\text{FPT}_{\\mathbf{\\hat{y}}}(c)$",
    FN="$\\text{FNT}_{\\mathbf{\\hat{y}}}(c)$",
    "min(FP,FN)"="$M_{\\mathbf{\\hat{y}}}(c)$")
  g.aum <- ggplot()+
    geom_vline(aes(
      xintercept=max.thresh),
      data=p.roc,
      color="grey")+
    theme_bw()+
    theme(
      panel.grid.minor=element_blank(),
      panel.spacing=grid::unit(0, "lines"))+
    coord_cartesian(expand=FALSE)+
    geom_rect(aes(
      xmin=min.thresh, xmax=max.thresh,
      ymin=0, ymax=value),
      color="grey",
      alpha=0.75,
      fill="grey",
      data=p.rect)+
    geom_text(aes(
      6, 1.5, label=sprintf("AUM=%.0f", aum)),
      data=p.auc)+
    geom_segment(aes(
      min.thresh, value,
      color=Variable, size=Variable,
      xend=max.thresh, yend=value),
      data=p.fp.fn)+
    geom_text(aes(
      ifelse(
        min.thresh == -Inf, max.thresh-1, ifelse(
          max.thresh == Inf, min.thresh+1, (min.thresh+max.thresh)/2)),
      ifelse(q %% 2, -2, -1),
      label=sprintf("q=%d",q)),
      vjust=-0.5,
      size=3,
      data=p.roc)+
    scale_color_manual(leg, values=err.colors, breaks=names(type.breaks), labels=type.breaks)+
    scale_size_manual(leg, values=err.sizes, breaks=names(type.breaks), labels=type.breaks)+
    scale_x_continuous(
      "Constant $c$ added to predicted values",
      breaks=seq(0, 12, by=2),
      limits=limits.vec)+
    scale_y_continuous(
      "Total label errors over all
$n$ labeled training examples",
breaks=seq(0, 10, by=2),
limits=c(NA, 11))
  g.list <- list(
    auc=g,
    "aum-nomath"=g.aum,
    aum=g.aum+
      geom_text(aes(
        thresh, 5.5,
        vjust=fcase(
          thresh %in% c(7,Inf), -0.5,
          default=1.2),
        label=sprintf(
          "$\\tau(\\mathbf{\\hat{y}})_{%d}=%s$",
          q, ifelse(
            is.finite(thresh),
            paste(thresh),
            paste(
              ifelse(thresh<0, "-", ""),
              "\\infty"))
        )),
        angle=90,
        data=p.roc[, data.table(
          q=c(0, .I),
          thresh=c(-Inf, max.thresh))])
  )
  for(plot.type in names(g.list)){
    f.tex <- sprintf(
      "figure-more-than-one-%s-%s.tex",
      m, plot.type)
    s <- 0.8
    tikz(
      f.tex,
      width=if(plot.type!="auc")5*s else 3*s,
      height=3*s,
      standAlone = TRUE)
    print(g.list[[plot.type]])
    dev.off()
    system(paste("pdflatex", f.tex))
  }
}

