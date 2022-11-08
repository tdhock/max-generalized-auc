source("packages.R")

auc.improved <- readRDS("auc.improved.rds")

norm.list <- list(
  l1=function(x)sum(abs(x)),
  l0=function(x)sum(x != 0))
norm.dt <- auc.improved[, {
  diff.wide <- roc[[1]][, lapply(.SD, diff), .SDcols=c("FPR","TPR")]
  diff.tall <- melt(diff.wide, measure=c("FPR", "TPR"))
  data.table(norm.name=names(norm.list))[, {
    norm.fun <- norm.list[[norm.name]]
    diff.tall[value!=0, .(
      norm.value=as.numeric(norm.fun(value))
    ), by=.(variable, sign=sign(value))]
  }, by=norm.name]
}, by=.(fold, set.name, initialization, pred.name)]
norm.wide <- dcast(
  norm.dt, 
  fold + set.name + initialization + variable + 
  sign + norm.name ~ pred.name , value.var="norm.value")
norm.wide[abs(improved - initial)>1e-5] ##??

moves.dt <- auc.improved[order(-min.thresh), {
  diff.wide <- roc[[1]][, lapply(.SD, diff), .SDcols=c("FPR","TPR")]
  diff.wide[, `:=`(
    fp.move = fcase(
      FPR > 0, "right",
      FPR < 0, "left"),
    tp.move = fcase(
      TPR > 0, "up",
      TPR < 0, "down"))
    ][, move := fcase(
      !is.na(fp.move) & !is.na(tp.move), paste0(tp.move, "+", fp.move),
      is.na(fp.move), tp.move,
      is.na(tp.move), fp.move
    )]
  diff.wide[, .(
    moves=as.numeric(.N),
    FPR=sum(FPR),
    TPR=sum(TPR)
  ), by=move]
}, by=.(fold, set.name, initialization, pred.name)]
moves.dt[is.na(move)]
moves.tall <- melt(moves.dt, measure=c("moves", "FPR", "TPR"))
moves.wide <- dcast(
  moves.tall, 
  fold + set.name + initialization + variable + move ~ pred.name)
moves.wide[order(initial-improved)]

auc.wide <- dcast(
  auc.improved,
  fold + set.name + initialization ~ pred.name , value.var="auc")
best <- auc.wide[initialization=="min.error"][order(initial-improved)][1]

on.vec <- c("fold", "set.name", "initialization")
auc.improved[best, on=on.vec]
moves.wide[best, .(
  move, variable, initial, improved, diff=round(initial-improved, 6)
), on=on.vec]

roc.dt <- auc.improved[, {
  roc[[1]][, .(
    thresh=c(-Inf,max.thresh), FPR=c(1,FPR), TPR=c(1,TPR)
  )]
}, by=.(fold, set.name, initialization, pred.name)]
roc.best <- roc.dt[best, on=on.vec]

regular.roc <- roc.dt[, {
  reg.dt <- data.table(
    FPR=cummin(FPR), TPR=cummin(TPR)
  )
  for(XPR in c("FPR","TPR")){
    reg.dt[, count := .N, by=XPR]
    reg.dt[, keep := TRUE]
    reg.dt[count>1, keep := c(TRUE, rep(FALSE,.N-2), TRUE), by=XPR]
    reg.dt <- reg.dt[keep==TRUE]
  }
  reg.dt
}, by=.(fold, set.name, initialization, pred.name)]
regular.auc <- regular.roc[, {
  AUC.WeightedROC <- WeightedROC::WeightedAUC(.SD)
  AUC.geometry <- geometry::polyarea(c(FPR,1), c(TPR,0))
  if(!isTRUE(all.equal(AUC.WeightedROC, AUC.geometry))){
    print(rbind(AUC.WeightedROC, AUC.geometry))
    print(.SD)
    browser()
  }
  data.table(auc.regular=AUC.WeightedROC)
}, by=.(fold, set.name, initialization, pred.name)]
auc.both <- auc.improved[regular.auc, on=.NATURAL]
auc.best <- auc.both[best, .(pred.name, auc, auc.regular), on=on.vec]
regular.best <- regular.roc[best, on=on.vec]
auc.best[, `:=`(diff=auc.regular-auc, y=c(0.3, 0.6))]

gg <- ggplot()+
  theme_bw()+
  theme(legend.position="none")+
  coord_equal()+
  geom_label(aes(
    1, y, 
    fill=pred.name,
    label=sprintf(
      "%s Full/color AUC=%.4f\nMonotonic/grey AUC=%.4f\n AUC Difference=%.4f",
      pred.name, auc, auc.regular, diff)),
    hjust=1,
    data=auc.best)+
  geom_path(aes(
    FPR, TPR, group=pred.name),
    size=2,
    color="grey50",
    data=regular.best)+
  geom_path(aes(
    FPR, TPR, color=pred.name),
    size=1,
    data=roc.best)+
  geom_point(aes(
    FPR, TPR, fill=pred.name),
    shape=21,
    data=regular.best)+
  scale_x_continuous(
    "False Positive Rate")+
  scale_y_continuous(
    "True Positive Rate")
png("figure-auc-improved.png", width=4, height=4, units="in", res=200)
print(gg)
dev.off()
