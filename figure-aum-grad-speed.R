source("packages.R")

csv.file.vec <- c(
  "Changepoint detection"="figure-aum-grad-speed-data.csv",
  "Binary classification"="figure-aum-grad-speed-binary-cpp-data.csv")
for(Problem in names(csv.file.vec)){
  csv.file <- csv.file.vec[[Problem]]
  timing.dt <- data.table::fread(csv.file)
  print(timing.dt)
}


timing.dt <- data.table::fread("figure-aum-grad-speed-data.csv")
timing.stats <- timing.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(N, pred.type, algorithm)]
some.stats <- timing.stats[pred.type=="pred.rnorm" & algorithm != "sort"]
some.stats[, Algorithm := c(
  "squared.hinge.each.example"="Squared\nHinge\nEach\nExample",
  aum="AUM")[algorithm]
  ]
algo.colors <- c(
  "Squared Hinge\nAll Pairs"="#A6CEE3",
  "Squared\nHinge\nEach\nExample"="#1F78B4",
  "Logistic"="#B2DF8A", #"#33A02C","#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  "AUM"="black"
)
gg <- ggplot()+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=Algorithm),
    alpha=0.5,
    data=some.stats)+
  geom_line(aes(
    N, median, color=Algorithm),
    data=some.stats)+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  scale_x_log10(
    "Number of predicted values",
    limits=c(10, 12000),
    breaks=c(10, 100, 1000, timing.stats[, max(N)]))+
  scale_y_log10(paste0("Computation time in seconds,
median line, min/max band over ",timing.stats[1, times], " timings"))+
  ggtitle("Changepoint detection")
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-aum-grad-speed-random.png", width=5, height=4, res=200, units="in")
print(dl)
dev.off()
gg <- ggplot()+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=algorithm),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N, median, color=algorithm),
    data=timing.stats)+
  facet_grid(. ~ pred.type)+
  scale_x_log10(
    "Number of predicted values",
    breaks=c(10, 100, 1000, timing.stats[, max(N)]))+
  scale_y_log10(paste0("Computation time in seconds,
median line, min/max band
over ",timing.stats[1, times], " timings"))
png("figure-aum-grad-speed.png", width=7, height=3, res=200, units="in")
print(gg)
dev.off()
