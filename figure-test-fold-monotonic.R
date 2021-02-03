library(data.table)
set.name <- "systematic"
set.path <- file.path("../neuroblastoma-data/data", set.name)
cv.type <- "R-3.6.0-profileSize"
folds.csv <- file.path(set.path, "cv", cv.type, "folds.csv")
fold.dt <- data.table::fread(folds.csv)
test.seqs <- fold.dt[fold==1]
errors.csv <- file.path(set.path, "errors.csv")
if(!file.exists(errors.csv)){
  errors.csv.xz <- paste0(errors.csv, ".xz")
  system(paste("unxz", errors.csv.xz))
}
err.dt <- data.table::fread(errors.csv)
err.test <- err.dt[test.seqs, on="sequenceID"]
err.tall <- data.table::melt(
  err.test,
  measure=c("fp", "fn"),
  id=c("sequenceID", "n.segments"))
err.tall[, diff := c(NA, diff(value)), by=.(sequenceID, variable)]
err.tall[!is.na(diff) & diff != 0, .(
  count=.N
), by=sequenceID][order(count)]
err.tall[sequenceID == "508_chr2"]

d <- function(data.name, cv.type, test.fold){
  data.table(data.name, cv.type, test.fold)
}
data.dt <- rbind(
  d("ATAC_JV_adipose", "equal_labels", 4),
  d("H3K27ac-H3K4me3_TDHAM_BP", "equal_labels", 2),
  d("H3K4me3_XJ_immune", "equal_labels", 2),
  d("H3K4me3_XJ_immune", "equal_labels", 4),
  d("systematic", "R-3.6.0-profileSize", 1))
meta.dt.list <- list()
for(data.i in 1:nrow(data.dt)){
  data.row <- data.dt[data.i]
  set.path <- file.path("../neuroblastoma-data/data", data.row$data.name)
  folds.csv <- file.path(set.path, "cv", data.row$cv.type, "folds.csv")
  fold.dt <- data.table::fread(folds.csv)
  inputs.csv <- file.path(set.path, "inputs.csv")
  inputs.dt <- data.table::fread(inputs.csv)
  inputs.mat <- as.matrix(inputs.dt[, -1, with=FALSE])
  keep <- apply(is.finite(inputs.mat), 2, all)
  errors.csv <- file.path(set.path, "evaluation.csv")
  if(!file.exists(errors.csv)){
    errors.csv.xz <- paste0(errors.csv, ".xz")
    system(paste("unxz", errors.csv.xz))
  }
  err.dt <- data.table::fread(errors.csv)
  train.seqs <- fold.dt[fold != data.row$test.fold]
  err.train <- err.dt[train.seqs, on="sequenceID"]
  err.tall <- data.table::melt(
    err.train,
    measure=c("fp", "fn"),
    id="sequenceID")
  err.tall[, diff := c(NA, diff(value)), by=.(sequenceID, variable)]
  break.dt <- err.tall[!is.na(diff) & diff != 0, .(
    count=.N
  ), by=sequenceID][order(count)]
  meta.dt.list[[data.i]] <- data.table(
    data.row,
    features=sum(keep),
    n.train=nrow(train.seqs),
    mean.breaks=mean(break.dt$count))
}
(meta.dt <- do.call(rbind, meta.dt.list))
meta.dt[, .(data.name, test.fold, features, n.train, mean.breaks)]
meta.tall <- data.table::melt(
  meta.dt,
  measure=c("features", "n.train", "mean.breaks"))
(meta.stats <- meta.tall[, .(
  min=min(value),
  max=max(value)
), by=variable])

