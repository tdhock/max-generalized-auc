source("packages.R")

result.list <- readRDS("figure-binary-test-auc-data.rds")

test.loss <- do.call(rbind, result.list)
gg <- ggplot()+
  geom_point(aes(
    auc, loss.name),
    data=test.loss)+
  facet_grid(select ~ ., labeller=label_both)+
  xlab("Test AUC")
png("figure-binary-test-auc.png", width=4, height=4, units="in", res=200)
print(gg)
dev.off()
