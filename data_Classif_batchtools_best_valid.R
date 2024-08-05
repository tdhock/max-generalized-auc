library(data.table)
library(ggplot2)
best.dt <- fread("data_Classif_batchtools_best_valid.csv")
best.wide <- dcast(
  best.dt,
  data.name + loss ~ .,
  list(mean, sd, length),
  value.var="auc"
)[, `:=`(
  Data = data.name,
  lo = auc_mean-auc_sd,
  hi = auc_mean+auc_sd
)][
, mid := (min(lo)+max(hi))/2
, by=Data][]
  
gg <- ggplot()+
  theme(
    panel.spacing=grid::unit(1, "lines"))+
  geom_point(aes(
    auc_mean, loss),
    shape=1,
    data=best.wide)+
  geom_segment(aes(
    lo, loss,
    xend=hi, yend=loss),
    data=best.wide)+
  geom_text(aes(
    auc_mean, loss,
    hjust=ifelse(auc_mean<mid, 0, 1),
    label=sprintf(
      "%.4f±%.4f", auc_mean, auc_sd)),
    size=3,
    vjust=1.2,
    data=best.wide)+
  facet_grid(. ~ Data, labeller=label_both, scales="free")+
  scale_x_continuous(
    "Max validation AUC (Mean ± SD over 4 random initializations)")
png("data_Classif_batchtools_best_valid.png", width=8, height=1.5, units="in", res=200)
print(gg)
dev.off()
