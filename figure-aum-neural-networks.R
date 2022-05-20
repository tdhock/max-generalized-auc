library(data.table)
library(ggplot2)
steps.dt <- data.table::fread(
  "figure-aum-neural-networks-data.csv"
)[epoch<10]
steps.dt[,which(is.na(out_value))[1]]
steps.dt[, out_num := as.numeric(out_value)]
steps.dt[, iteration := epoch*(max(step)+1)+step]
## analysis of missing values, seems to happen when loss
## underflows. Suggests to use sum rather than mean.
count.values <- function(DT){
  dcast(
    DT,
    loss + lr+set_name ~ data_set,
    fun.aggregate=length)
}
count.values(steps.dt[is.na(out_num)])
count.values(steps.dt[out_value==0])
steps.dt[is.na(out_num)]
steps.dt[is.na(out_num) & !is.na(out_value)]
str(steps.dt)

one <- steps.dt[
  data_set=="MNIST" & out_name=="AUC" & set_name=="validation"]
ggplot()+
  facet_grid(lr ~ seed, labeller=label_both)+
  geom_line(aes(
    iteration, out_num, color=loss),
    data=one)

## check if all initializations were the same.
init.dt <- dcast(
  steps.dt[iteration==0&out_name=="AUC"],
  seed+data_set+set_name~.,
  fun.aggregate=list(min,max),
  value.var="out_num")
init.dt[1e-5 < out_num_max-out_num_min]#should be empty.

valid.auc <- steps.dt[set_name=="validation" & out_name=="AUC"]
by.vars <- c("loss","seed","data_set")
valid.vars <- c(by.vars,"lr")
test.vars <- c(valid.vars, "iteration")
selected.dt <- valid.auc[, .SD[which.max(out_num)], by=by.vars]
select.test <- selected.dt[,test.vars,with=FALSE]
all.test.auc <- steps.dt[set_name=="test" & out_name=="AUC"]
select.test.auc <- all.test.auc[select.test, on=names(select.test)]
ggplot()+
  facet_grid(. ~ data_set, labeller=label_both)+
  geom_point(aes(
    out_num, loss),
    data=select.test.auc)
wide.test.auc <- dcast(
  select.test.auc,
  loss + data_set ~ .,
  fun.aggregate=list(mean,sd),
  value.var="out_num")
show.names <- c(
  balanced="logistic.weighted",
  logistic="logistic.unweighted",
  AUM="AUM.count",
  AUM_rate="AUM.rate")
wide.test.auc[, loss.name := ifelse(
  loss %in% names(show.names), show.names[loss], loss)]
levs <- wide.test.auc[data_set=="MNIST"][order(out_num_mean), loss.name]
wide.test.auc[, loss.fac := factor(loss.name, levs)]
gg <- ggplot()+
  facet_grid(. ~ data_set, labeller=label_both)+
  geom_point(aes(
    out_num_mean, loss.fac),
    shape=1,
    data=wide.test.auc)+
  geom_segment(aes(
    out_num_mean+out_num_sd, loss.fac,
    xend=out_num_mean-out_num_sd, yend=loss.fac),
    data=wide.test.auc)+
  scale_y_discrete(
    "Loss function")+
  scale_x_continuous(paste(
    "Test AUC",
    "(Mean +/- SD over 4 random initializations of neural network weights)"))
png(
  "figure-aum-neural-networks-test-auc.png",
  width=7, height=1.3, units="in", res=200)
print(gg)
dev.off()

p.wide <- dcast(select.test.auc, data_set+seed ~ loss,value.var="out_num")
p.tall <- melt(p.wide, measure=c("AUM","balanced","logistic"))
p.tall[, {
  t.test(
    AUM_rate, value, alternative="greater", paired=TRUE
  )[c("estimate","p.value")]
}, keyby=.(data_set,variable)]
    
select.valid <- selected.dt[,valid.vars,with=FALSE]
select.valid.auc <- valid.auc[select.valid,on=names(select.valid)]
gg <- ggplot()+
  coord_cartesian(ylim=c(0.5,1))+
  geom_line(aes(
    iteration,out_num,color=loss),
    data=select.valid.auc)+
  facet_grid(data_set ~ seed, labeller=label_both)
print(gg)
png(
  "figure-aum-neural-networks-best-valid-auc-curves.png",
  width=10, height=4, units="in", res=200)
print(gg)
dev.off()

subtrain.loss <- steps.dt[set_name=="subtrain" & out_name=="loss"]
selected.subtrain <- subtrain.loss[select.valid,on=names(select.valid)]
selected.subtrain[iteration>100 & data_set=="FashionMNIST" & loss=="AUM_rate"]
ggplot()+
  geom_line(aes(
    iteration, out_num, color=factor(seed)),
    data=selected.subtrain)+
  facet_grid(loss ~ data_set, labeller=label_both, scales="free")+
  scale_y_log10()
