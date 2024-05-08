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
for(data.name in names(result.list)){
  atime.result <- result.list[[data.name]]
  atime.ref <- atime::references_best(atime.result)
  it.step.list[[data.name]] <- data.table(
    data.name,
    atime.ref$measurements[unit %in% c("mean.it.per.step","steps")]
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
(ldt <- rbindlist(loss.dt.list))
it.step <- parse_expr.name(rbindlist(it.step.list)[
, Data := data.name
])

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

select.dt <- ldt[
, seeds := .N
, by=.(data.name,N,set,step.number,maxIt)
][
  seeds==max(seeds)
][, .(
  algos=length(unique(maxIt))
), by=.(data.name,N)
][
  algos>1
][
, .SD[which.max(N)]
, by=data.name]
linetype.scale <- scale_linetype_manual(values=c(
  quadratic=2,
  linear=1,
  min.aum=3))
(it.step.stats <- dcast(
  it.step,
  Data + N + maxIt + unit ~ .,
  list(mean,sd,length),
  value.var="empirical"
)[empirical_length==max(empirical_length)])
gg <- ggplot()+
  geom_vline(aes(
    xintercept=N),
    color="grey",
    size=1,
    data=data.table(select.dt)[, Data := data.name])+
  scale_y_log10(
    "")+
  scale_x_log10(
    "Number of rows of training data")+
  theme_bw()+
  linetype.scale+
  geom_line(aes(
    N, empirical_mean, linetype=maxIt),
    size=1,
    data=it.step.stats)+
  geom_ribbon(aes(
    N,
    ymin=empirical_mean-empirical_sd,
    group=maxIt,
    ymax=empirical_mean+empirical_sd),
    alpha=0.4,
    fill="black",
    color=NA,
    data=it.step.stats)+
  facet_grid(unit~Data,scales="free")
png("data_Classif_figure_units.png", width=8, height=3, units="in", res=300)
print(gg)
dev.off()

select.loss <- ldt[
  select.dt, on=.(data.name,N)
][
  seed==1
][
, Data := sprintf(
  "Data: %s, N=%d",
  data.name, N)
][]
select.min <- select.loss[
  set=="validation",
  .SD[which.min(aum)],
  by=Data]
gg <- ggplot()+
  theme_bw()+
  scale_y_log10(
    "AUM")+
  scale_x_log10(
    "Number of gradient descent steps with line search")+
  geom_line(aes(
    step.number, aum,
    size=set,
    linetype=maxIt,
    color=set),
    data=select.loss)+
  geom_point(aes(
    step.number, aum,
    fill=step,
    color=set),
    shape=21,
    data=data.table(select.min,step="min"))+
  facet_grid(.~Data,scales="free")+
  scale_fill_manual(values=c(min="white"))+
  linetype.scale+
  scale_color_manual(values=c(
    subtrain="red",
    validation="black"))+
  scale_size_manual(values=c(
    subtrain=2,
    validation=1))
png("data_Classif_figure_subtrain_validation.png", width=8, height=3, units="in", res=300)
print(gg)
dev.off()
