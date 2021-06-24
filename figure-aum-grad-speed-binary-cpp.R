source("packages.R")

timing.dt <- data.table::fread("figure-aum-grad-speed-binary-cpp-data.csv")
timing.stats <- timing.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(N, prediction.order, Algorithm=sub(" All", "\nAll", algorithm))]

algo.colors <- c(
  "Squared Hinge\nAll Pairs"="#A6CEE3",
  "squared hinge each example"="#1F78B4",
  "Logistic"="#B2DF8A", #"#33A02C","#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  "AUM"="black"
)
unsorted <- timing.stats[prediction.order == "unsorted"]
max.N <- unsorted[, max(N)]
gg <- ggplot()+
  ggtitle("Binary classification")+
  theme(legend.position='none')+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=Algorithm),
    alpha=0.5,
    data=unsorted)+
  geom_line(aes(
    N, median, color=Algorithm),
    data=unsorted)+
  scale_x_log10(
    "Number of predictions",
    limits=c(10, max.N*10))+
  scale_y_log10(
    "Computation time in seconds,
median line, min/max band, 10 timings") +
  geom_dl(aes(N, median, label = Algorithm, color = Algorithm), 
          method = "right.polygons",
          data = unsorted[Algorithm == "Squared Hinge\nAll Pairs",]) +
  geom_dl(aes(N, median, label = Algorithm, color = Algorithm), 
          method = "right.polygons",
          data = unsorted[Algorithm != "Squared Hinge\nAll Pairs",])
png("figure-aum-grad-speed-binary-cpp-algos.png", width=5, height=4, res=200, units="in")
print(gg)
dev.off()

gg <- ggplot()+
  facet_grid(. ~ algorithm, labeller=label_both)+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=prediction.order),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N, median, color=prediction.order),
    data=timing.stats)+
  scale_x_log10(
    "Number of predictions",
    limits=c(10, max.N*10))+
  scale_y_log10(
    "Computation time in seconds,
median line, min/max band, 10 timings")
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-aum-grad-speed-binary-cpp.png", width=10, height=3, res=200, units="in")
print(dl)
dev.off()
