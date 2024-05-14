library(data.table)
scaled.csv.vec <- grep(
  "seed=",
  Sys.glob("data_Classif_scaled/*.csv"),
  invert=TRUE,
  value=TRUE)
for(csv.i in seq_along(scaled.csv.vec)){
  scaled.csv <- scaled.csv.vec[[csv.i]]
  scaled.dt <- fread(scaled.csv)
  names(scaled.dt)
  scaled.ncol <- ncol(scaled.dt)-2
  for(seed in 1:4){
    set.seed(seed)
    initial.weights <- rnorm(scaled.ncol)
    out.dt <- data.table(weight=initial.weights)
    out.csv <- sub(".csv", paste0("_seed=",seed,".csv"), scaled.csv)
    fwrite(out.dt,out.csv)
  }
}
