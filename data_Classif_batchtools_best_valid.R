library(data.table)
library(ggplot2)
best.dt <- fread("data_Classif_batchtools_best_valid.csv")
loss2show <- rev(c(
  Logistic="Logistic/Cross-entropy\n(classic baseline)",
  SquaredHinge="All Pairs Squared Hinge\n(recent alternative)",
  AUM="AUM=Area Under Min(FP,FN)\n(proposed complex loss)",
  NULL))
best.wide <- dcast(
  best.dt,
  data.name + loss ~ .,
  list(mean, sd, length),
  value.var="auc"
)[, `:=`(
  Loss = factor(loss2show[loss], loss2show),
  Data = data.name,
  lo = auc_mean-auc_sd,
  hi = auc_mean+auc_sd
)][
, mid := (min(lo)+max(hi))/2
, by=Data][]
gg <- ggplot()+
  theme(
    plot.margin=grid::unit(c(0,1,0,0), "lines"),
    panel.spacing=grid::unit(1.5, "lines"))+
  geom_point(aes(
    auc_mean, Loss),
    shape=1,
    data=best.wide)+
  geom_segment(aes(
    lo, Loss,
    xend=hi, yend=Loss),
    data=best.wide)+
  geom_text(aes(
    auc_mean, Loss,
    hjust=ifelse(auc_mean<mid, -0.1, 1.1),
    vjust=ifelse(loss=="AUM", 0.5, 1.3),
    label=sprintf(
      "%.4f±%.4f", auc_mean, auc_sd)),
    size=3,
    data=best.wide)+
  facet_grid(. ~ Data, labeller=label_both, scales="free")+
  scale_x_continuous(
    "Max validation AUC (Mean ± SD over 4 random initializations)")
png("data_Classif_batchtools_best_valid.png", width=10, height=1.6, units="in", res=200)
print(gg)
dev.off()
