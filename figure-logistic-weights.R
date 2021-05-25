source("packages.R")

zip.X.list <- list()
zip.y.list <- list()
for(set in c("train", "test")){
  f <- sprintf("zip.%s.gz", set)
  if(!file.exists(f)){
    u <- paste0("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/", f)
    download.file(u, f)
  }
  zip.dt <- data.table::fread(f)
  y.vec <- zip.dt[[1]]
  is.01 <- y.vec %in% 0:1
  y01.dt <- data.table(label=y.vec[is.01])
  y01.dt[, cum := 1:.N, by=label]
  max.dt <- y01.dt[, .(max=max(cum)), by=label]
  keep <- y01.dt$cum <= min(max.dt[["max"]])
  zip.y.list[[set]] <- y01.dt[keep, label]
  zip.X.list[[set]] <- as.matrix(zip.dt[is.01, -1, with=FALSE][keep,])
}
(y.tab <- sapply(zip.y.list, table))

train.set.list <- list(
  full=list(X=zip.X.list[["train"]], y=zip.y.list[["train"]]))
some.props <- c(0.01, 0.05)
prop.pos.vec <- sort(unique(c(some.props, 1-some.props, 0.5)))
##want p/(p + n) = 0.05 => 0.05*(p+n) = p => 0.05p + 0.05n = p => 0.05n = 0.95p => p = 0.05 / 0.95n
min.prop.pos <- min(prop.pos.vec)
min.n.pos <- as.integer(min.prop.pos/(1-min.prop.pos) * y.tab["0", "train"])
min.total <- min.n.pos + y.tab["0", "train"]
c(min.n.pos, y.tab["0", "train"])/min.total
N.obs <- 1000
train.y.dt <- data.table(label=zip.y.list[["train"]])
train.y.dt[, i := 1:.N]
test.y <- zip.y.list[["test"]]
result.dt.list <- list()
for(prop.pos in prop.pos.vec){
  prop.dt <- rbind(
    data.table(prop=prop.pos, label=1),
    data.table(prop=1-prop.pos, label=0))
  prop.dt[, class.N := as.integer(N.obs*prop) ]
  prop.dt[, weight := 1/class.N]
  for(seed in 1:3){
    cat(sprintf("prop=%f seed=%d\n", prop.pos, seed))
    set.seed(seed)
    index.dt <- prop.dt[train.y.dt, on="label"][, .(
      i=.SD[sample(1:.N), i[1:class.N] ]
    ), by=.(label, weight, class.N)]
    seed.i <- index.dt[["i"]]
    seed.y <- zip.y.list[["train"]][seed.i]
    seed.X <- zip.X.list[["train"]][seed.i,]
    weight.list <- list(
      identity=rep(1, length(seed.y)),
      balanced=index.dt[["weight"]])
    for(weight.name in names(weight.list)){
      weight.vec <- weight.list[[weight.name]]
      fit <- glmnet::cv.glmnet(seed.X, seed.y, weight.vec, family="binomial")
      seed.pred <- predict(fit, zip.X.list[["test"]])
      roc.df <- WeightedROC::WeightedROC(seed.pred, test.y)
      seed.pred.class <- ifelse(0<seed.pred, 1, 0)
      accuracy <- mean(seed.pred.class == test.y)
      auc <- WeightedROC::WeightedAUC(roc.df)
      result.dt.list[[paste(prop.pos, seed, weight.name)]] <- data.table(
        prop.pos, seed, weight.name, accuracy, auc)
    }
  }
}
(result.dt <- do.call(rbind, result.dt.list))

result.tall <- melt(result.dt, measure.vars=c("accuracy", "auc"))
result.tall[, percent.positive.labels := factor(prop.pos*100)]
ggplot()+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_point(aes(
    percent.positive.labels, value, color=weight.name),
    data=result.tall)

result.stats <- result.tall[, .(
  max=max(value),
  q75=quantile(value, 0.75),
  median=median(value),
  q25=quantile(value, 0.25),
  min=min(value),
  seeds=.N
), by=.(variable, prop.pos, percent.positive.labels, weight.name)]
gg <- ggplot()+
  ggtitle(paste0(
    "cv.glmnet run on data sets with same number of observations, N=",
    nrow(seed.X),
    "\nand with different proportions of positive labels"))+
  facet_grid(variable ~ ., labeller = label_both, scales="free")+
  geom_ribbon(aes(
    prop.pos, ymin=min, ymax=max, fill=weight.name),
    alpha=0.5,
    data=result.stats)+
  geom_line(aes(
    prop.pos, median, color=weight.name),
    data=result.stats)+
  scale_x_continuous(
    "Proportion positive labels in train set",
    breaks=unique(result.stats[["prop.pos"]]))+
  ylab("Accuracy or AUC of predictions
on a test set of 50% positive
and 50% negative labels")
png("figure-logistic-weights.png", width=10, height=3, units="in", res=200)
print(gg)
dev.off()
