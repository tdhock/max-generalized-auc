library(data.table)
library(ggplot2)
constant.dt <- fread("data_Classif_batchtools.csv")
constant.dt[1]
best.constant <- constant.dt[
  set_name=="validation",
  .SD[which.max(auc)],
  keyby=.(data.name,N,loss,seed)
]
select.dt <- best.constant[, .(data.name,N,loss,seed,lr,step_number)]
select.best <- constant.dt[select.dt, on=.NATURAL]
fwrite(best.constant, "data_Classif_batchtools_best_valid.csv")

best.constant[
, algorithm := paste0(loss,", constant")
][]
ls.dt <- fread("data_Classif_line_search.csv")
best.ls <- ls.dt[
  set=="validation",
  .SD[which.max(auc)],
  keyby=.(data.name,N,seed)
][
, algorithm := "AUM, line search"
][]

constant.dt[step_number==0 & loss=="AUM" & lr==1, .(data.name,seed,set_name,auc)]
ls.dt[step.number==0 & objective=="max.auc" & set.obj=="subtrain", .(data.name,seed,set,auc)]

best.both <- rbind(
  best.ls[, .(data.name, N, seed, algorithm, auc)],
  best.constant[, .(data.name, N, seed, algorithm, auc)])

best.wide <- dcast(
  best.both,
  data.name + N + algorithm ~ .,
  list(mean,sd,length),
  value.var="auc")
ggplot()+
  geom_segment(aes(
    auc_mean+auc_sd, algorithm,
    xend=auc_mean-auc_sd, yend=algorithm),
    data=best.wide)+
  geom_point(aes(
    auc_mean, algorithm),
    data=best.wide)+
  facet_grid(. ~ data.name + N, scales="free")

ls.one <- ls.dt[data.name=="CIFAR10" & seed==1]
ggplot()+
  facet_wrap( ~ objective + set.obj, labeller=label_both, scales="free")+
  geom_line(aes(
    step.number, auc, color=set),
    data=ls.one)+
  scale_x_log10()+
  scale_y_log10()
    
