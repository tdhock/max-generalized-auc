source("packages.R")
data(neuroblastoma,package="neuroblastoma")
data(neuroblastomaProcessed, package="penaltyLearning")

e <- function(label, profile.id, chromosome){
  Label <- paste0("\n", label)
  data.table(
    label, Label,
    profile.id=factor(profile.id),
    chromosome=factor(chromosome))
}
select.dt <- rbind(
  e("positive", 4, 2),
  e("negative", 513, 3)
  ##e("negative", 7, 4)
  ##e("negative", 76, 2)
  ##e("negative", 409, 4)
  ##e("negative", 485, 2)
  ##e("negative", 490, 2) #more separated
)[
, Label := factor(Label, Label)
][]
makeSeq <- function(DT){
  DT[, SeqID := paste0(profile.id,".",chromosome)][]
}
nb.dt <- makeSeq(data.table(neuroblastoma$profiles))
label.dt <- makeSeq(data.table(neuroblastoma$ann))
some.profiles <- nb.dt[select.dt, on=c("profile.id","chromosome")]
some.labels <- label.dt[select.dt, on=c("profile.id","chromosome")]
some.err <- makeSeq(neuroblastomaProcessed$errors[select.dt, {
  .(profile.id, chromosome,
    fp, fn, possible.fp, possible.fn,
    min.log.lambda=-max.log.lambda,
    max.log.lambda=-min.log.lambda,
    errors, labels,
    n.segments,
    label, 
    Label)
}, on=list(profile.id, chromosome)
][, `:=`(
  fp.diff = c(NA, diff(fp)),
  fn.diff = c(NA, diff(fn))
), by=label])
some.pred <- some.err[
  is.na(fp.diff) | fp.diff!=0 | fn.diff!=0
][
, pred.log.lambda := round(ifelse(
  min.log.lambda== -Inf,
  max.log.lambda-0.5,
  (min.log.lambda+max.log.lambda)/2),1)
][]
err.sizes <- c(
  "min(FP,FN)"=1,
  FP=3,
  FN=2)
err.colors <- c(
  "min(FP,FN)"="black",
  FP="red",
  FN="deepskyblue")
some.err.tall <- melt(
  some.err,
  measure.vars=c("fp","fn"),
  variable.name="var.lower"
)[
, variable := toupper(var.lower)
][]
leg <- "Error type"
pred.dt <- rbind(
  data.table(SeqID="513.3", pred=2.55),
  data.table(SeqID="4.2", pred=-3))
gg.err <- ggplot()+
  theme_bw()+
  theme(
    legend.position=c(0.2,0.25),
    legend.background=element_rect(color="black",fill="grey90"),
    panel.spacing=grid::unit(0, "lines"))+
  facet_grid(SeqID ~ ., labeller=label_both)+
  geom_vline(aes(
    xintercept=pred.log.lambda),
    data=some.pred,
    size=1,
    color="grey")+
  geom_vline(aes(
    xintercept=pred),
    data=pred.dt,
    size=1,
    color="violet")+
  geom_text(aes(
    pred.log.lambda, 0.5,
    label=n.segments-1),
    data=some.pred,
    color="black")+
  geom_text(aes(
    x=2.5,y=0.5,label="changes predicted"),
    data=some.pred[n.segments==4],
    color="black")+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value,
    color=variable, size=variable),
    data=some.err.tall)+
  scale_y_continuous(
    "Label errors",
    breaks=c(0,1),
    limits=c(-0.2, 1.2))+
  scale_color_manual(leg,values=err.colors)+
  scale_size_manual(leg,values=err.sizes)+
  scale_x_continuous(
    "Predicted value f(x)",
    breaks=seq(-10, 10, by=1))+
  coord_cartesian(xlim=c(-3, 5))
png("figure-aum-convexity-new-profiles.png", 3.5, 2.5, units="in", res=200)
print(gg.err)
dev.off()

some.profiles[, data.i := seq(1,.N),by=SeqID]
(nb.segs <- some.profiles[, {
  cum.vec <- cumsum(c(0, logratio))
  d <- diff(position)/2
  between <- position[-1]-d
  data.start.pos <- c(position[1]-d[1], between)
  data.end.pos <- c(between, position[.N]+d[.N-1])
  max.segments <- max(some.pred$n.segments)
  fit <- jointseg::Fpsn(logratio, max.segments)
  end.t <- t(fit$t.est)
  end.dt <- data.table(
    end=as.integer(end.t),
    n.segments=as.integer(col(end.t))
  )[!is.na(end)]
  end.dt[, start := c(0, end[-.N])+1, by=n.segments]
  end.dt[, mean := (cum.vec[end+1]-cum.vec[start])/(end-start+1)]
  end.dt[, `:=`(
    start.pos=data.start.pos[start],
    end.pos=data.end.pos[end]
  )]
}, by=SeqID][])
getY <- function(rel){
  M <- -1
  m <- -2
  d <- M-m
  rel*d+m
}
some.pred[
, yrel := seq(1,.N)/.N
, by = SeqID
][, `:=`(
  ymax=getY(yrel-0.05),
  ymin=getY(yrel-0.2)
)]
(nb.changes <- nb.segs[, .SD[-1, .(
  change=start-0.5
)], by=.(SeqID,n.segments)
][some.pred, on=c("SeqID","n.segments"),nomatch=0L])
some.label.data <- some.profiles[some.labels, .(
  min.i=min(data.i),
  max.i=max(data.i),
  label=label[1],
  annotation=annotation[1]
), on=.(SeqID,position>min, position<max), by=.EACHI
][, .(SeqID, min.i, max.i, label, annotation)]
err.list <- penaltyLearning::labelError(
  some.pred,
  some.label.data,
  nb.changes,
  change.var="change",
  label.vars=c("min.i","max.i"),
  problem.vars="SeqID")
err.abbrev <- c(
  "false negative"="FN",
  "false positive"="FP",
  "correct"="None")
nb.changes[, x := max(change), by=SeqID]
gg <- ggplot()+
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing=grid::unit(0, "lines"))+
  scale_fill_manual(
    "label",
    values=c(
      positive="violet",
      negative="orange"))+
  geom_rect(aes(
    xmin=min.i, xmax=max.i,
    ymin=-Inf, ymax=Inf,
    fill=label),
    alpha=0.5,
    data=some.label.data)+
  geom_rect(aes(
    xmin=min.i, xmax=max.i,
    ymin=ymin, ymax=ymax),
    fill=NA,
    color="grey",
    size=0.5,
    data=err.list$label.errors[, Error := err.abbrev[status]])+
  geom_rect(aes(
    xmin=min.i, xmax=max.i,
    ymin=ymin, ymax=ymax,
    color=Error,
    linetype=Error),
    fill=NA,
    size=1,
    data=err.list$label.errors[, Error := err.abbrev[status]])+
  geom_text(aes(
    min.i, -0.9, label=paste(
      label, "label\n",
      ifelse(
        label=="positive",
        "1 change expected",
        "0 changes expected"))),
    hjust=0,
    size=2.5,
    vjust=0,
    data=some.label.data)+
  geom_text(aes(
    min.i, ymin, label=Error,
    color=Error),
    hjust=1.1,
    vjust=0,
    data=err.list$label.errors)+
  geom_point(aes(
    data.i, logratio),
    color="grey50",
    shape=1,
    data=some.profiles)+
  facet_grid(SeqID ~ ., labeller=label_both)+
  geom_segment(aes(
    change, ymax,
    xend=change, yend=ymin),
    size=0.5,
    data=nb.changes)+
  geom_text(aes(
    Inf, ymin, label=sprintf(
      "%d change%s predicted ",
      n.segments-1,
      ifelse(n.segments-1==1,"","s"))),
    hjust=1,
    size=2.5,
    vjust=0,
    data=err.list$label.errors)+
  scale_linetype_manual(
    "Error",
    values=c(
      None=0,
      FN=3,
      FP=1))+
  scale_color_manual(
    values=err.colors)+
  scale_x_continuous(
    "Data index along time/space sequence",
    limits=c(-10,NA))+
  scale_y_continuous(
    "Data value",
    breaks=seq(-1,1,by=0.5))
png("figure-aum-convexity-new.png", 4.5, 2.5, units="in", res=400)
print(gg)
dev.off()
