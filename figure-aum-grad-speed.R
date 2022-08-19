source("packages.R")

replacements <- c(
  "squared.hinge.each.example"="Squared\nHinge\nEach\nExample",
  "Squared Hinge All Pairs"="Squared Hinge\nAll Pairs",
  aum="AUM")
row_fun <- function(...){
  form.list <- as.list(match.call()[-1])
  sym <- sapply(form.list, is.symbol)
  names(form.list)[sym] <- paste(form.list[sym])
  form.list[sym] <- NA
  make_row <- function(){}
  formals(make_row) <- form.list
  body(make_row) <- as.call(lapply(c("data.frame", names(form.list)), as.symbol))
  make_row
}
f <- row_fun(Problem, file.csv, col.name, col.value)
csv.file.dt <- rbind(
  f("Changepoint detection","figure-aum-grad-speed-data.csv","pred.type","pred.rnorm"),
  f("Binary classification","figure-aum-grad-speed-binary-cpp-data.csv","prediction.order","unsorted"))

problem.dt.list <- list()
for(file.i in 1:nrow(csv.file.dt)){
  csv.file.row <- csv.file.dt[file.i,]
  timing.dt <- data.table::fread(csv.file.row[["file.csv"]])
  problem.dt.list[[file.i]] <- timing.dt[
    algorithm != "sort" &
    get(csv.file.row[["col.name"]]) == csv.file.row[["col.value"]],
    data.table(
      csv.file.row, N, seconds,
      Algorithm=ifelse(
        algorithm%in%names(replacements),
        replacements[algorithm],
        algorithm))
  ]
}
(problem.dt <- do.call(rbind, problem.dt.list))

algo.colors <- c(
  "Squared Hinge\nAll Pairs"="#A6CEE3",
  "Squared\nHinge\nEach\nExample"="#1F78B4",
  "Logistic"="#B2DF8A", #"#33A02C","#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  "AUM"="black"
)
problem.stats <- problem.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(Problem, N, Algorithm)]
mydl <- function(data){
  geom_dl(aes(N, median, label = Algorithm, color = Algorithm), 
          method = list(cex=0.8, "right.polygons"),
          data = data)
}
breaks <- 10^seq(0, 6)
dl <- ggplot()+
  theme(legend.position="none")+
  ##facet_wrap(. ~ Problem, labeller=label_both, scales="free")+
  facet_grid(. ~ Problem, labeller=label_both, scales="free", space="free")+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=Algorithm),
    alpha=0.5,
    data=problem.stats)+
  geom_line(aes(
    N, median, color=Algorithm),
    data=problem.stats)+
  geom_blank(aes(
    N*5, median, color=Algorithm),
    data=problem.stats)+
  geom_blank(aes(
    x,y),
    data=data.table(x=100, y=10^c(-5, -2)))+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  mydl(problem.stats[Algorithm != "Squared Hinge\nAll Pairs",])+
  mydl(problem.stats[Algorithm == "Squared Hinge\nAll Pairs",])+
  scale_x_log10(
    "n = number of predicted values = size of gradient vector",
    breaks=breaks,
    labels=sprintf("%.e", breaks))+
  scale_y_log10(paste0("Computation time in seconds,
median line, min/max band over ",problem.stats[1, times], " timings"),
breaks=10^seq(-6, 0))
png("figure-aum-grad-speed-both.png", width=7, height=3.2, res=200, units="in")
print(dl)
dev.off()

binary.stats <- problem.stats[Problem=="Binary classification"]
dl <- ggplot()+
  theme(legend.position="none")+
  geom_ribbon(aes(
    N, ymin=min, ymax=max, fill=Algorithm),
    alpha=0.5,
    data=binary.stats)+
  geom_line(aes(
    N, median, color=Algorithm),
    data=binary.stats)+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  mydl(binary.stats[Algorithm != "Squared Hinge\nAll Pairs",])+
  mydl(binary.stats[Algorithm == "Squared Hinge\nAll Pairs",])+
  scale_x_log10(
    "n = number of predicted values = size of gradient vector",
    breaks=breaks,
    limits=c(1e1, 2e6),
    labels=sprintf("%.e", breaks))+
  scale_y_log10(paste0("Computation time in seconds,
median line, min/max band over ",problem.stats[1, times], " timings"),
breaks=10^seq(-6, 0))
dl
png("figure-aum-grad-speed-binary.png", width=7, height=3.4, res=200, units="in")
print(dl)
dev.off()

timing.dt <- data.table::fread("figure-aum-grad-speed-data.csv")
timing.stats <- timing.dt[, .(
  max=max(seconds),
  median=median(seconds),
  min=min(seconds),
  times=.N
), by=.(N, pred.type, algorithm)]
some.stats <- timing.stats[pred.type=="pred.rnorm" & algorithm != "sort"]
some.stats[, Algorithm := replacements[algorithm] ]
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
    "Number of predicted values in gradient computation",
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
