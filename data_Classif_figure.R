library(ggplot2)
library(data.table)
library(atime)
load("data_Classif.RData")

plot(ref.list$FishSonar)
plot(ref.list$N)

loss.dt.list <- list()
it.step.list <- list()
for(data.name in names(result.list)){
  atime.result <- result.list[[data.name]]
  atime.ref <- ref.list[[data.name]]
  it.step.list[[data.name]] <- data.table(
    data.name,
    atime.ref$measurements[unit %in% c("mean.it.per.step","steps")]
  )
  loss.dt <- atime.result$measurements[
  , fit[[1]]$loss, by=.(N,expr.name)
  ][
  , maxIt := sub(".*=", "", expr.name)
  ][
   , N.i := .GRP, by=N
  ][]
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
ldt <- rbindlist(loss.dt.list)
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

select.dt <- ldt[, .(
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
it.step <- rbindlist(it.step.list)[
, maxIt := sub(".*=", "", expr.name)
][
, Data := data.name
][]
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
    N, empirical, linetype=maxIt),
    size=1,
    data=it.step)+
  facet_grid(unit~Data,scales="free")
png("data_Classif_figure_units.png", width=6, height=3, units="in", res=300)
print(gg)
dev.off()

select.loss <- ldt[
  select.dt, on=.(data.name,N)
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
png("data_Classif_figure_subtrain_validation.png", width=6, height=3, units="in", res=300)
print(gg)
dev.off()
