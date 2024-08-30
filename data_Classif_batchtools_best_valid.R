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
  list(mean, sd, length, min, max),
  value.var=c("auc","step_number")
)[, `:=`(
  Loss = factor(loss2show[loss], loss2show),
  Data = data.name,
  lo = auc_mean-auc_sd,
  hi = auc_mean+auc_sd
)][
, mid := (min(lo)+max(hi))/2
, by=Data][]
dput(RColorBrewer::brewer.pal(3,"Dark2"))
loss.colors <- c("black", "#D95F02", "#7570B3")
names(loss.colors) <- loss2show
p <- function(Data,x,y){
  data.table(Data,x,y)
}
blank.dt <- rbind(
  p("CIFAR10",0.8,100),
  p("MNIST",0.99,c(10,100000)),
  p("STL10",0.8,30))
gg <- ggplot()+
  theme_bw()+
  theme(
    plot.margin=grid::unit(c(0,1,0,0), "lines"),
    legend.key.spacing.y=grid::unit(1, "lines"),
    axis.text.x=element_text(angle=30, hjust=1),
    panel.spacing=grid::unit(1.5, "lines"))+
  geom_blank(aes(x, y), data=blank.dt)+
  geom_point(aes(
    auc_mean, step_number_mean,
    color=Loss),
    shape=1,
    data=best.wide)+
  geom_segment(aes(
    auc_min, step_number_mean,
    color=Loss,
    xend=auc_max, yend=step_number_mean),
    data=best.wide)+
  geom_segment(aes(
    auc_mean, step_number_min,
    color=Loss,
    xend=auc_mean, yend=step_number_max),
    data=best.wide)+
  facet_wrap("Data", nrow=1, labeller=label_both, scales="free")+
  scale_color_manual(
    values=loss.colors)+
  scale_y_log10(
    "Gradient descent epochs\n(using best learning rate)")+
  scale_x_continuous(
    "Best validation AUC (dot=mean, segments=range over 4 random initializations)")
print(gg)
## 100000 max steps
## lr=10^seq(-4,5)
png("data_Classif_batchtools_best_valid_scatter.png", width=10, height=2, units="in", res=200)
print(gg)
dev.off()

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
