source("packages.R")
timing.dt <- data.table::fread("figure-aum-grad-speed-data.csv")
timing.stats <- timing.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(N, pred.type, algorithm)]
some.stats <- timing.stats[pred.type=="pred.rnorm"]
some.stats[, Algorithm := gsub("[.]", "\n", algorithm)]
gg <- ggplot()+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=Algorithm),
    alpha=0.5,
    data=some.stats)+
  geom_line(aes(
    N, median, color=Algorithm),
    data=some.stats)+
  scale_x_log10(
    "Number of predicted values",
    limits=c(10, 12000),
    breaks=c(10, 100, 1000, timing.stats[, max(N)]))+
  scale_y_log10(paste0("Computation time in seconds,
median line, min/max band
over ",timing.stats[1, times], " timings"))
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-aum-grad-speed-random.png", width=4, height=3, res=200, units="in")
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
