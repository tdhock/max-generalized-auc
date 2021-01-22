source("packages.R")
abbrev.vec <- c(
  "prediction vector"="optimized",
  "linear model weights"="train")
both.list <- list()
for(optimization.variable in names(abbrev.vec)){
  fname.rds <- paste0("figure-aum-", abbrev.vec[[optimization.variable]], "-data.rds")
  result.list <- readRDS(fname.rds)
  for(data.type in c("roc", "auc")){
    dt <- result.list[[data.type]]
    dt[, model := sub("improved", "optimized", pred.name)]
    both.list[[data.type]][[optimization.variable]] <- data.table(
      optimization.variable, dt)
  }
}
for(data.type in names(both.list)){
  both.list[[data.type]] <- do.call(rbind, both.list[[data.type]])
}
both.list$auc[, `:=`(x=c(0.25), y=ifelse(model=="initial", 0.75, 0.5))]
gg <- ggplot()+
  facet_grid(. ~ optimization.variable, labeller=label_both)+
  geom_path(aes(
    FPR, TPR, color=model),
    data=both.list$roc)+
  geom_point(aes(
    FPR, TPR, color=model),
    fill="white",
    shape=21,
    data=both.list$auc)+
  geom_segment(aes(
    x, y,
    xend=FPR, yend=TPR,
    color=model),
    size=0.25,
    data=both.list$auc)+
  geom_label(aes(
    x, y, color=model,
    label=sprintf(
      "%s aum=%.2f\nerrors=%d auc=%.2f",
      model, aum, errors, auc)),
    size=3,
    hjust=0,
    vjust=1,
    data=both.list$auc)+
  coord_equal()+
  guides(color="none")+
  theme(panel.spacing=grid::unit(0.5, "cm"))
png("figure-aum-train-both.png", width=5.8, height=3.3, units="in", res=200)
print(gg)
dev.off()
