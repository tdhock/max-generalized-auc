source("packages.R")

aum.csv <- Sys.glob(
  "../neuroblastoma-data/data/*/cv/*/testFolds/*/linear-model-aum.csv")
all.it <- data.table(out.csv=aum.csv)[, {
  data.table::fread(out.csv)
}, by=out.csv]

count.dt <- all.it[, .(
  count=.N
), by=.(data.name, cv.type, test.fold)]
stopifnot(nrow(count.dt)==length(aum.csv))

subtrain.it <- all.it[set=="subtrain"]
subtrain.it[, diff := c(NA, diff(aum)), by=.(init.name, data.name, test.fold, seed)]
subtrain.it[, .(init.name, data.name, test.fold, iteration, aum, diff)]
subtrain.it[diff>1e-6]
gg <- ggplot()+
  ggtitle("check if train AUM decreases")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_line(aes(
    iteration, aum,
    group=paste(seed, init.name)),
    data=subtrain.it)+
  facet_grid(init.name + data.name + test.fold ~ ., scales="free", labeller=label_both)
print(gg)

validation.it <- all.it[set=="validation"]
ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  scale_y_log10()+
  geom_line(aes(
    iteration, aum, color=init.name,
    group=paste(seed, init.name)),
    data=validation.it)+
  geom_point(aes(
    iteration, aum, color=init.name,
    group=paste(seed, init.name)),
    data=validation.it[
     ,
       .SD[which.min(aum)],
       by=.(data.name, test.fold, init.name, seed)])+
  facet_grid(data.name + test.fold ~ ., scales="free", labeller=label_both)

valid.best.ids <- all.it[
  set=="validation",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, cv.type, test.fold, init.name, seed)]
test.best.ids <- all.it[
  set=="test",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, cv.type, test.fold, init.name, seed)]

## model selection.
test.it1 <- all.it[set=="test" & iteration==1]
test.selected <- all.it[set=="test"][valid.best.ids, on=names(valid.best.ids)]
test.best <- all.it[set=="test"][test.best.ids, on=names(test.best.ids)]

## compare with best predictions (no linear model).
test.show <- rbind(
  data.table(iterations="initial", test.it1),
  data.table(iterations="best.linear", test.best),
  data.table(iterations="selected", test.selected))
ifac <- function(x)factor(
  x, c("initial", "selected", "best.linear"))
test.show[, Iterations := ifac(iterations)]
gg <- ggplot()+
  ggtitle("Test AUM, selected=min valid aum, best=min test aum, max it=50")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    aum, Iterations, color=factor(test.fold)),
    shape=1,
    data=test.show)+
  scale_y_discrete(drop=FALSE)+
  facet_grid(
    init.name ~ data.name + cv.type,
    scales="free", labeller=label_both)
print(gg)

gg <- ggplot()+
  ggtitle("Test AUM, selected=min valid aum, best=min test aum, max it=50")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    aum, Iterations, color=factor(test.fold)),
    shape=1,
    data=test.show[init.name=="IntervalRegressionCV"])+
  scale_y_discrete(drop=FALSE)+
  facet_grid(
    init.name ~ data.name + cv.type,
    scales="free", labeller=label_both)
print(gg)

test.show[, neg.auc := -auc]
test.show.tall <- melt(
  test.show[init.name=="IntervalRegressionCV"],
  measure.vars=c("neg.auc", "error.percent", "aum"),
  variable.name="metric")
test.iCV <- dcast(
  test.show.tall,
  data.name + cv.type + test.fold + metric + seed ~ iterations)
test.iCV.tall <- melt(
  test.iCV,
  measure.vars=c("best.linear", "selected"),
  variable.name="iteration")

test.iCV.tall[, improvement := value - initial]
imp.stats <- test.iCV.tall[, .(
  median=median(improvement),
  p.value=tryCatch({
    t.test(initial, value)[["p.value"]]
  }, error=function(e){
    NA_real_
  })
), by=.(data.name, cv.type, test.fold, iteration, metric)][!is.na(p.value)][median < 0][order(p.value)]
top10 <- imp.stats[metric=="aum" & iteration=="best.linear"][1:min(.N, 10)]
some.types <- unique(top10[, .(data.name, cv.type, test.fold)])

test.show[, Data.name := paste0("\n", data.name)]
gg <- ggplot()+
  ggtitle("Optimizing train AUM can reduce test AUM if number of iterations is chosen correctly")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    aum, Iterations),
    shape=1,
    data=test.show[init.name=="IntervalRegressionCV"][some.types, on=names(some.types)])+
  scale_y_discrete(drop=FALSE)+
  facet_grid(
    . ~ Data.name + cv.type + test.fold,
    scales="free", labeller=label_both)+
  xlab("Test AUM, each dot is a different random seed/initialization for IntervalRegressionCV")
png("figure-linear-model-test-analyze.png", width=20, height=2.5, units="in", res=100)
print(gg)
dev.off()
