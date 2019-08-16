source("packages.R")

nb.comb <- readRDS("neuroblastomaProcessed.combinations.rds")

worst <- nb.comb$auc[which.max(auc)]
worst.combo <- nb.comb$combos[worst, .(panel, interval), on=list(combo.i)]
nb.comb$segs.min.err[, pred.log.lambda := ifelse(
  min.log.lambda == -Inf, max.log.lambda-worst$size, ifelse(
    max.log.lambda == Inf, min.log.lambda+worst$size, mid.log.lambda))]
nb.comb$segs.min.err[, interval := ifelse(
  is.finite(mid.log.lambda), "finite", "infinite")]
pred.dt <- nb.comb$segs.min.err[worst.combo, on=list(panel, interval)]
L <- penaltyLearning::ROChange(
  nb.comb$some.err, pred.dt, c("panel"))
L$auc

L$auc.polygon[, row := 1:.N]
ggplot()+
  geom_polygon(aes(
    FPR, TPR),
    fill="red",
    color="black",
    alpha=0.5,
    data=L$auc.polygon)+
  geom_text(aes(
    FPR, TPR, label=row),
    data=L$auc.polygon)

sel.dt <- L$auc.polygon[row>1, .(row, first=1)]
setkey(sel.dt, first, row)
L$auc.polygon[, row0 := row]
setkey(L$auc.polygon, row, row0)
cum.poly <- foverlaps(sel.dt, L$auc.polygon, nomatch=0L)
cum.poly[, added := ifelse(i.row==row, "new", "old")]
lim <- c(-0.2, 1.2)
gg <- ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=cum.poly)+
  geom_text(aes(
    FPR, TPR, label=row, color=added),
    data=cum.poly)+
  scale_color_manual(values=c("new"="red", old="black"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("i.row", nrow=2)+
  ##facet_grid(. ~ i.row)+
  coord_equal(xlim=lim, ylim=lim)+
  scale_x_continuous(breaks=seq(0, 1, by=0.5), labels=c("0", "0.5", "1"))+
  scale_y_continuous(breaks=seq(0, 1, by=0.5))
png("figure-neuroblastomaProcessed-combinations-worst.png", 12, 3, units="in", res=100)
print(gg)
dev.off()

gg <- ggplot()+
  geom_point(aes(
    aub, auc),
    color="black",
    shape=21,
    size=5,
    fill=NA,
    data=nb.comb$auc)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ size)
print(gg)
nb.comb$auc[order(aub), .(auc, aub, size, combo.i)]

rfac <- 5
nb.comb$auc[, round.aub := round(aub*rfac)/rfac]
nb.comb$auc[, round.auc := round(auc, 4)]
aub.count <- nb.comb$auc[, list(
  combos=.N
), by=list(aub=round.aub, size, round.auc)]
gg <- ggplot()+
  geom_hline(aes(
    yintercept=yint),
    data=data.table(yint=1),
    color="grey50")+
  geom_point(aes(
    aub, round.auc, fill=combos),
    shape=21,
    size=5,
    data=aub.count)+
  scale_fill_gradient(low="white", high="red")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(size ~ .)+
  geom_text(aes(
    aub, round.auc, label=combos),
    size=3,
    data=aub.count)+
  scale_y_continuous(
    "Area under ROC curve",
    breaks=seq(0, 1.2, by=0.2))+
  scale_x_continuous(
    "Area under both TP and FP curves")
print(gg)
png("figure-neuroblastomaProcessed-combinations-scatter.png", 12, 9, units="in", res=100)
print(gg)
dev.off()

auc.count <- nb.comb$auc[, list(
  combos=.N
), by=list(n.finite, size, round.auc)]
gg <- ggplot()+
  geom_tile(aes(
    n.finite, round.auc, fill=combos),
    data=auc.count)+
  geom_point(aes(
    n.finite, auc),
    color="black",
    shape=21,
    size=5,
    fill=NA,
    data=worst)+
  scale_fill_gradient(low="white", high="red")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ size)+
  geom_text(aes(
    n.finite, round.auc, label=combos),
    size=3,
    data=auc.count)+
  scale_x_continuous(
    "Number of predictions in finite min error interval (other predictions in the infinite min error interval)",
    breaks=unique(auc.count$n.finite))
png("figure-neuroblastomaProcessed-combinations.png", 12, 3, units="in", res=100)
print(gg)
dev.off()

