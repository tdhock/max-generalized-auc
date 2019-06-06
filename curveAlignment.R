source("packages.R")

future::plan("multiprocess")

some.probs <- data.table(prob.dir=paste0(
  "ATAC_JV_adipose/samples/AC1/MSC",
  c(83, 91),
  "/problems/chrX-37148256-49242997"))
some.labels <- some.probs[, {
  labels.bed <- paste0(
    "../feature-learning-benchmark/data/",
    prob.dir, "/labels.bed")
  fread(
    labels.bed,
    col.names=c(
      "chrom",
      "labelStart",
      "labelEnd",
      "annotation"))
}, by=list(prob.dir)]
min.labelStart <- min(some.labels$labelStart)
max.labelEnd <- max(some.labels$labelEnd)
label.range <- max.labelEnd-min.labelStart
expand <- label.range/20

profiles.dt <- some.probs[, {
  from <- min.labelStart-expand
  to <- max.labelEnd+expand
  s <- sub(":", "-", prob.dir)
  bg <- paste0(
    "../feature-learning-benchmark/data/",
    s, "/coverage.bedGraph")
  gz <- paste0(bg, ".gz")
  dt <- fread(
    gz,
    col.names=c("chrom", "chromStart", "chromEnd", "coverage"))
  fwrite(dt, bg, col.names=FALSE, sep="\t")
  dt[from < chromEnd & chromStart < to]
}, by=list(prob.dir)]
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
ggplot()+
  ggtitle(
    "Noisy coverage data and labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(prob.dir ~ ., scales="free")+
  geom_tallrect(aes(
    xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
    data=some.labels,
    color="grey")+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    data=profiles.dt,
    color="grey50")+
  scale_x_continuous(breaks=seq(4e4, 5e4, by=5))

win <- function(windowStart, windowEnd){
  data.table(windowStart, windowEnd)
}
win.dt <- rbind(
  win(43447, 43457),
  win(43502, 43512))*1000
win.dt[, window := 1:.N]
setkey(profiles.dt, chromStart, chromEnd)
setkey(some.labels, labelStart, labelEnd)
setkey(win.dt, windowStart, windowEnd)
win.profiles <- foverlaps(profiles.dt, win.dt, nomatch=0L)
win.labels <- foverlaps(some.labels, win.dt, nomatch=0L)
gg <- ggplot()+
  ggtitle(
    "Noisy coverage data and labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(prob.dir ~ window, scales="free", labeller=label_both)+
  geom_tallrect(aes(
    xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
    data=win.labels,
    color="grey")+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    data=win.profiles,
    color="grey50")+
  scale_x_continuous(breaks=seq(4e4, 5e4, by=5))
print(gg)

win.err.list <- list()
win.segs.list <- list()
win.selection.list <- list()
for(seq.i in 1:nrow(some.probs)){
  s <- some.probs[seq.i]
  pdir <- paste0(
    "../feature-learning-benchmark/data/",
    sub(":", "-", s$prob.dir))
  L <- PeakSegPipeline::problem.target(pdir, 1)
  plabels <- win.labels[prob.dir==s$prob.dir]
  plabels[, chromStart := labelStart]
  plabels[, chromEnd := labelEnd]
  selection.dt <- data.table(penaltyLearning::modelSelection(
    L$models, "total.loss", "peaks"))
  win.selection.list[[paste(seq.i)]] <- data.table(
    prob.dir=s$prob.dir,
    selection.dt)
  for(model.i in 1:nrow(selection.dt)){
    model <- selection.dt[model.i]
    pen.str <- paste(model$penalty)
    pen.info <- PeakSegDisk::problem.PeakSegFPOP(pdir, pen.str)
    seg.dt <- data.table(prob.dir=s$prob.dir, model, pen.str, pen.info$segments)
    setkey(seg.dt, chromStart, chromEnd)
    over.dt <- foverlaps(seg.dt, win.dt, nomatch=0L)
    peak.dt <- over.dt[status=="peak"]
    e <- PeakError::PeakErrorChrom(peak.dt, plabels)
    win.err.list[[paste(seq.i, model.i)]] <-
      data.table(
        prob.dir=s$prob.dir,
        window=plabels$window,
        model, pen.str, e)
    win.segs.list[[paste(seq.i, model.i)]] <- over.dt
  }
}
win.selection <- do.call(rbind, win.selection.list)
win.segs <- do.call(rbind, win.segs.list)
win.err <- do.call(rbind, win.err.list)
win.segs[, segStart := ifelse(chromStart<windowStart, windowStart, chromStart)]
win.segs[, segEnd := ifelse(windowEnd<chromEnd, windowEnd, chromEnd)]

possible <- fread(
  "../feature-learning-benchmark/labeled_problems_possible_errors.csv")
possible[, prob.dir := sub(":", "-", prob.dir)]
eval.dt <- possible[win.selection, nomatch=0L, on=list(prob.dir)]

out.list <- list(
  evaluation=eval.dt,
  segments=win.segs,
  errors=win.err,
  problems=some.probs,
  profiles=win.profiles,
  labels=win.labels)
saveRDS(out.list, "curveAlignment.rds")

