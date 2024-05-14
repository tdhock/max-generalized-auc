library(ggplot2)
library(data.table)
library(atime)
load("data_Classif.RData")
loss.dt.list <- list()
it.step.list <- list()
parse_expr.name <- function(df)nc::capture_first_df(
  df, 
  expr.name=list(
    "lm ",
    nc::field("maxIt","=",".*?"),
    ",",
    nc::field("seed","=",".*")))
show.units <- c(
  "mean.it.per.step","steps",
  ##"seconds",
  NULL)
for(data.name in names(result.list)){
  atime.result <- result.list[[data.name]]
  atime.ref <- atime::references_best(atime.result)
  it.step.list[[data.name]] <- data.table(
    data.name,
    atime.ref$measurements[unit %in% show.units]
  )
  loss.dt <- parse_expr.name(atime.result$measurements[
  , fit[[1]]$loss, by=.(N,expr.name)
  ][
  , N.i := .GRP, by=N
  ])
  loss.dt[step.number==0]#verify initialization same.
  min.dt <- loss.dt[, .SD[which.min(aum)], by=.(set,N)]
  loss.dt.list[[data.name]] <- data.table(data.name, loss.dt)
  ggplot()+
    geom_hline(aes(
      yintercept=aum, color=maxIt),
      data=min.dt)+
    geom_line(aes(
      step.number, aum, color=maxIt, size=maxIt),
      data=loss.dt)+
    facet_grid(N ~ set)+
    scale_color_manual(values=c(
      linear="black",
      quadratic="red",
      min.aum="deepskyblue"))+
    scale_size_manual(values=c(
      linear=3,
      min.aum=2,
      quadratic=1))
  ggplot()+
    scale_y_log10()+
    geom_hline(aes(
      yintercept=aum, color=maxIt),
      data=min.dt)+
    geom_line(aes(
      step.number, aum, color=maxIt, size=maxIt),
      data=loss.dt)+
    facet_wrap( ~ N+set,ncol=2,scales="free")+
    scale_color_manual(values=c(
      linear="black",
      quadratic="red",
      min.aum="deepskyblue"))+
    scale_size_manual(values=c(
      linear=3,
      min.aum=2,
      quadratic=1))
}
maxIt.colors <- c(
  max="black",
  quadratic="#1B9E77",
  min.aum="#D95F02",
  linear="#7570B3")
(ldt <- rbindlist(loss.dt.list))[
, maxIt.fac := factor(maxIt,names(maxIt.colors))
][]
(it.step <- parse_expr.name(rbindlist(it.step.list))[
, Data := data.name
][
, maxIt.fac := factor(maxIt,names(maxIt.colors))
][])
if(FALSE){
  ggplot()+
    scale_y_log10()+
    geom_line(aes(
      step.number, aum, color=maxIt, size=maxIt),
      data=ldt)+
    facet_grid(N.i ~ data.name+set,scales="free")+
    scale_color_manual(values=c(
      linear="black",
      quadratic="red",
      min.aum="deepskyblue"))+
    scale_size_manual(values=c(
      linear=3,
      min.aum=2,
      quadratic=1))
}

select.dt <- ldt[
, seeds := .N
, by=.(data.name,N,set,step.number,maxIt.fac)
][
  seeds==max(seeds) & maxIt.fac!="quadratic"
][, .(
  algos=length(unique(maxIt.fac))
), by=.(data.name,N)
][
  algos==2
][
, .SD[which.max(N)]
, by=data.name]
dput(RColorBrewer::brewer.pal(3,"Dark2"))
(it.step.stats <- dcast(
  it.step,
  Data + N + maxIt.fac + unit ~ .,
  list(mean,sd,length),
  value.var="empirical"
)[empirical_length==max(empirical_length)])
leg <- "line search\niterations"
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  geom_vline(aes(
    xintercept=N),
    color="grey",
    size=1,
    data=data.table(select.dt)[, Data := data.name])+
  scale_y_log10(
    "")+
  scale_x_log10(
    "Number of rows of training data")+
  geom_line(aes(
    N, empirical_mean, color=maxIt.fac),
    size=1,
    data=it.step.stats)+
  scale_color_manual(leg,values=maxIt.colors)+
  scale_fill_manual(leg,values=maxIt.colors)+
  geom_ribbon(aes(
    N,
    ymin=empirical_mean-empirical_sd,
    group=maxIt.fac,
    fill=maxIt.fac,
    ymax=empirical_mean+empirical_sd),
    alpha=0.4,
    color=NA,
    data=it.step.stats)+
  ##facet_grid(unit~Data,scales="free")
  facet_wrap(~unit+Data,scales="free",nrow=length(show.units))
it.max <- it.step[, .(
  iterations=(N-1)*N/2,
  maxIt.fac=factor("max",names(maxIt.colors)),
  unit="mean.it.per.step"
), by=.(Data,N)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  geom_line(aes(
    N, iterations, color=maxIt.fac),
    data=it.max)+
  geom_vline(aes(
    xintercept=N),
    color="grey",
    size=1,
    data=data.table(select.dt)[, Data := data.name])+
  scale_y_log10(
    "")+
  scale_x_log10(
    "Number of rows of training data")+
  geom_line(aes(
    N, empirical, color=maxIt.fac, group=paste(maxIt.fac,seed)),
    data=it.step)+
  scale_color_manual(leg,values=maxIt.colors,drop=FALSE)+
  facet_wrap(~unit+Data,scales="free",nrow=length(show.units))
png("data_Classif_figure_units_all.png", width=20, height=2.5*length(show.units), units="in", res=300)
print(gg)
dev.off()

show.data <- c("MNIST","FashionMNIST","CIFAR10","STL10")
some <- function(DT)DT[Data %in% show.data]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  geom_line(aes(
    N, iterations, color=maxIt.fac),
    data=some(it.max))+
  geom_vline(aes(
    xintercept=N),
    color="grey",
    size=1,
    data=some(data.table(select.dt)[, Data := data.name]))+
  scale_y_log10(
    "")+
  scale_x_log10(
    "Number of rows of training data")+
  geom_line(aes(
    N, empirical_mean, color=maxIt.fac),
    size=1,
    data=some(it.step.stats))+
  scale_color_manual(leg,values=maxIt.colors,drop=FALSE)+
  scale_fill_manual(leg,values=maxIt.colors,guide="none",drop=FALSE)+
  geom_ribbon(aes(
    N,
    ymin=empirical_mean-empirical_sd,
    group=maxIt.fac,
    fill=maxIt.fac,
    ymax=empirical_mean+empirical_sd),
    alpha=0.4,
    color=NA,
    data=some(it.step.stats))+
  facet_grid(unit~Data,scales="free")
png("data_Classif_figure_units.png", width=8, height=1.5*length(show.units), units="in", res=300)
print(gg)
dev.off()

select.loss <- ldt[
  select.dt, on=.(data.name,N)
][
  seed==1
][
, Data := sprintf(
  "%s\nN=%d",
  data.name, N)
][maxIt!="quadratic"]
select.min <- select.loss[
  set=="validation",
  .SD[which.min(aum)],
  by=.(Data,maxIt)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  scale_y_log10(
    "AUM")+
  scale_x_log10(
    "Number of gradient descent steps with line search")+
  ## geom_line(aes(
  ##   step.number, aum, group=paste(set,maxIt)),
  ##   color="white",
  ##   size=2.5,
  ##   data=select.loss)+
  geom_line(aes(
    step.number, aum,
    size=set,
    color=maxIt,
    linetype=set),
    data=select.loss)+
  ## geom_point(aes(
  ##   step.number, aum,
  ##   color=maxIt,
  ##   fill=step),
  ##   shape=21,
  ##   data=data.table(select.min,step="min"))+
  geom_point(aes(
    step.number, aum,
    color=maxIt),
    shape=21,
    fill="white",
    data=data.table(select.min,step="min"))+
  facet_grid(.~Data,scales="free")+
  scale_fill_manual(values=c(min="white"))+
  scale_linetype_manual(values=c(
    validation="solid",
    subtrain="dotted"))+
  scale_color_manual(values=maxIt.colors)+
  scale_size_manual(values=c(
    subtrain=1,
    validation=0.7))
png("data_Classif_figure_subtrain_validation_all.png", width=16, height=3, units="in", res=300)
print(gg)
dev.off()

show <- function(DT)DT[data.name %in% show.data]
set.linetype.vec <- c(
  validation="solid",
  subtrain="dotted")
select.loss[, Set := factor(set, names(set.linetype.vec))]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  scale_y_log10(
    "AUM")+
  scale_x_log10(
    "Number of gradient descent steps with line search")+
  ## geom_line(aes(
  ##   step.number, aum, group=paste(set,maxIt)),
  ##   color="white",
  ##   size=2.5,
  ##   data=select.loss)+
  geom_line(aes(
    step.number, aum,
    size=Set,
    color=maxIt,
    linetype=Set),
    data=show(select.loss))+
  ## geom_point(aes(
  ##   step.number, aum,
  ##   color=maxIt,
  ##   fill=step),
  ##   shape=21,
  ##   data=data.table(show(select.min),step="min"))+
  geom_point(aes(
    step.number, aum,
    color=maxIt),
    fill="white",
    shape=21,
    data=data.table(show(select.min),step="min"))+
  ##facet_grid(.~Data,scales="free")+
  facet_wrap(.~Data,scales="free",nrow=1)+
  scale_fill_manual(values=c(min="white"))+
  scale_linetype_manual(
    values=set.linetype.vec)+
  scale_color_manual(
    leg,
    values=maxIt.colors)+
  scale_size_manual(values=c(
    subtrain=1,
    validation=0.7))
png("data_Classif_figure_subtrain_validation.png", width=8, height=2, units="in", res=300)
print(gg)
dev.off()

select.max <- select.loss[
  set=="validation",
  .SD[which.max(auc)],
  by=.(Data,maxIt)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(1,"lines"))+
  scale_y_log10(
    "Area Under ROC Curve (AUC)")+
  scale_x_log10(
    "Number of gradient descent steps with line search")+
  geom_line(aes(
    step.number, auc,
    size=Set,
    color=maxIt,
    linetype=Set),
    data=show(select.loss))+
  geom_point(aes(
    step.number, auc,
    color=maxIt),
    fill="white",
    shape=21,
    data=data.table(show(select.max),step="max"))+
  facet_wrap(.~Data,scales="free",nrow=1)+
  scale_fill_manual(values=c(min="white"))+
  scale_linetype_manual(
    values=set.linetype.vec)+
  scale_color_manual(
    leg,
    values=maxIt.colors)+
  scale_size_manual(values=c(
    subtrain=1,
    validation=0.7))
png("data_Classif_figure_subtrain_validation_AUC.png", width=8, height=2, units="in", res=300)
print(gg)
dev.off()

system("cd /projects/genomic-ml/ && unpublish_data max-generalized-auc && publish_data projects/max-generalized-auc")
