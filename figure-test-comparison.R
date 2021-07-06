library(dplyr)
library(ggplot2)
library(data.table)


if(!file.exists("figure-test-comparison.csv"))
{
  source("figure-test-comparison-data.R")
}

test.aum.dt.combined <- data.table::fread("figure-test-comparison.csv")


png("figure-test-auc-comparison.png", width = 26, height = 3, res = 200, units = "in")


ggplot(data = test.aum.dt.combined) +
  geom_point(aes(x = `Test AUC`, y = algorithm), size = 5) +
  #ggtitle(c(data.name, cv.type, test.fold)) +
  facet_grid(.~data.name + new.test.fold, scales = "free") +
  theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25))

dev.off()

png("figure-test-aum-comparison.png", width = 26, height = 3, res = 200, units = "in")


ggplot(data = test.aum.dt.combined) +
  geom_point(aes(x = `Test AUM`, y = algorithm), size = 5) +
  #ggtitle(c(data.name, cv.type, test.fold)) +
  facet_grid(.~data.name + new.test.fold, scales = "free") +
  theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25))

dev.off()

