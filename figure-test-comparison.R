library(dplyr)
library(ggplot2)
library(data.table)

# testFold.vec <- Sys.glob("../AUM meeting figures/28Jan2021/*/linear-model-post-cv-aum.csv")
# precv.vec <- Sys.glob("../AUM meeting figures/28Jan2021/*/linear-model-cv-aum.csv")
# postcv.vec <- Sys.glob("../AUM meeting figures/28Jan2021/*/linear-model-post-cv-aum.csv")

testFold.vec <- c(file.path("neuroblastoma-data/data/ATAC_JV_adipose/cv/equal_labels/testFolds/4"),
                  file.path("neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/2"),
                  file.path("neuroblastoma-data/data/systematic/cv/R-3.6.0-profileSize/testFolds/1"),
                  file.path("neuroblastoma-data/data/H3K4me3_XJ_immune/cv/equal_labels/testFolds/2"),
                  file.path("neuroblastoma-data/data/H3K4me3_XJ_immune/cv/equal_labels/testFolds/4"))

#png.testFold.vec <- sapply(testFold.vec, function(x){file.path(paste0(substring(x, 1, nchar(x) - 4), ".png"))})



# lapply(1:3, function(x){pdf_convert(testFold.vec[x], png.testFold.vec[x] )})
# for(index in 1:length(testFold.vec))
# {
#   try(pdftools::pdf_convert(testFold.vec[index], filenames = png.testFold.vec[index]), silent = TRUE)
# }

aum.list <- list()
auc.list <- list()
test.aum.dt.list <- list()

for(curr.index in 1:length(testFold.vec))
{
      
      testFold.path <- testFold.vec[curr.index]
      #Filename and directory initializations
      cv.path <- dirname(dirname(testFold.path))
      folds.csv <- file.path(cv.path, "folds.csv")
      cv.type <- basename(cv.path)
      test.fold <- basename(testFold.path)
      data.dir <- dirname(dirname(cv.path))
      data.name <- basename(data.dir)
      data.list <- list()
      
      #Initialize a generic plot title for later use
      plot.title <- paste0("Data Name = ", data.name, ", cv.type = ", cv.type, ", test.fold = ", test.fold)
      cv.algos <- c("random")
      n.seeds <- 4
      
      cv.dt <- data.table::fread(file.path(testFold.path, "linear-model-cv-aum.csv")) %>% 
        mutate(seed = as.factor(seed))
      post.cv.dt <- data.table::fread(file.path(testFold.path, "linear-model-post-cv-aum.csv")) %>%
        mutate(seed = as.factor(seed))
      
      #Make the data.table of the mean aum/auc with respect to each seed, iteration,
      #data set, and cross validation algorithm
      mean.cv.dt <- cv.dt[, .(mean.aum=mean(aum), mean.auc=mean(auc)), by=.(seed, iteration, set, cv.algo)]
      
      #Make the data.table of the min/max aum/auc with respect to each seed, iteration,
      #data set, and cross validation algorithm, and include their corresponding iterations
      min.cv.dt <- cv.dt[, .(min.aum.iteration = which.min(aum), min.aum = min(aum),
                             max.auc.iteration = which.max(auc), max.auc = max(auc)), 
                         by=.(seed, val.fold, set, cv.algo)]
      
      min.mean.cv.dt <- mean.cv.dt[, .(min.aum.iteration = which.min(mean.aum), min.aum = min(mean.aum),
                                       max.auc.iteration = which.max(mean.auc), max.auc = max(mean.auc)), 
                                   by=.(seed, set, cv.algo) ]
      
      
      
      best.linear.dt <- post.cv.dt[set == "test", .(`Test AUM` = min(aum),
                                                    `Test AUC` = max(auc),
                                                    algorithm = "best.linear",
                                                    aum.iteration = which.min(aum),
                                                    auc.iteration = which.max(auc)),
                                   by = .(seed)]
      
      #Create the selected iteration data.table for later use in 
      #training on the whole training dataset
      selected.iter.dt <- min.mean.cv.dt[set == "validation"]
      
      #Create selected iteration data table lists with respect to picking
      #the minimum aum or the minimum auc.
      selected.aum.dt.list <- list()
      selected.auc.dt.list <- list()
      for(curr.cv.algo in cv.algos)
      {
        for(curr.seed in 1:n.seeds)
        {
          curr.dt <- post.cv.dt[seed == curr.seed & set == "test"]
          min.aum.iteration <- selected.iter.dt[seed == curr.seed & cv.algo == curr.cv.algo]$min.aum.iteration
          max.auc.iteration <- selected.iter.dt[seed == curr.seed & cv.algo == curr.cv.algo]$max.auc.iteration
          
          selected.aum.dt.list[[paste0(curr.seed, curr.cv.algo)]] <- data.table(`Test AUM` = curr.dt$aum[min.aum.iteration],
                                                                                `Test AUC` = curr.dt$auc[min.aum.iteration],
                                                                                algorithm = paste0("Min.Valid.AUM"),
                                                                                seed = curr.seed,
                                                                                aum.iteration = min.aum.iteration,
                                                                                auc.iteration = min.aum.iteration)
          
          selected.auc.dt.list[[paste0(curr.seed, curr.cv.algo)]] <- data.table(`Test AUM` = curr.dt$aum[max.auc.iteration],
                                                                                `Test AUC` = curr.dt$auc[max.auc.iteration],
                                                                                algorithm = paste0("Max.Valid.AUC"),
                                                                                seed = curr.seed,
                                                                                aum.iteration = max.auc.iteration,
                                                                                auc.iteration = max.auc.iteration)
        }
      }
      
      selected.aum.dt <- do.call(rbind, selected.aum.dt.list)
      selected.auc.dt <- do.call(rbind, selected.auc.dt.list)
      
      # selected.dt <-
      #   test.post.cv.dt[, .(aum = last(aum),
      #                       type = "selected"), 
      #                   by = seed]
      
      #Create the data.table containing the initial aum/auc for
      #the dataset
      initial.dt <-
        post.cv.dt[set == "test", .(`Test AUM` = first(aum),
                                    `Test AUC` = first(auc),
                                    algorithm = "Initial",
                                    aum.iteration = 1,
                                    auc.iteration = 1),  
                   by = seed]
      
      #Combine intitial, selected, and best.linear results for plotting purposes
      test.aum.dt <- rbind(selected.aum.dt, selected.auc.dt, initial.dt)
      
      test.aum.dt.list[[curr.index]] <- test.aum.dt %>% mutate(data.name = data.name, 
                                                               cv.type = cv.type,
                                                               test.fold = test.fold)
      # aum.list[[curr.index]] <- (ggplot(data = test.aum.dt) +
      #         geom_point(aes(x = aum, y = type, color = seed)) +
      #         ggtitle(plot.title))
      # 
      # auc.list[[curr.index]] <- (ggplot(data = test.aum.dt) +
      #         geom_point(aes(x = auc, y = type, color = seed)) +
      #         ggtitle(c(data.name, cv.type, test.fold)))
}

test.aum.dt.combined <- do.call(rbind, test.aum.dt.list) %>%
  mutate(new.test.fold = paste0("Test Fold ", test.fold))

png("figure-test-auc-comparison.png", width = 26, height = 3, res = 200, units = "in")


ggplot(data = test.aum.dt.combined) +
  geom_point(aes(x = `Test AUC`, y = algorithm, color = seed), size = 5) +
  #ggtitle(c(data.name, cv.type, test.fold)) +
  facet_grid(.~data.name + new.test.fold, scales = "free") +
  theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25))

# ggsave(
#   "figure-test-auc-comparison.png",
# ggplot(data = test.aum.dt.combined) +
#   geom_point(aes(x = `Test AUC`, y = algorithm, color = seed)) +
#   #ggtitle(c(data.name, cv.type, test.fold)) +
#   facet_grid(.~data.name + test.fold, scales = "free") +
#   theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25)),
# width = 26,
# height = 3,
# dpi = 1800
# )

dev.off()

png("figure-test-aum-comparison.png", width = 26, height = 3, res = 200, units = "in")

# ggsave(
#   "figure-test-aum-comparison.png",
# ggplot(data = test.aum.dt.combined) +
#   geom_point(aes(x = `Test AUM`, y = algorithm, color = seed)) +
#   #ggtitle(c(data.name, cv.type, test.fold)) +
#   facet_grid(.~data.name + test.fold, scales = "free") +
#   theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25)),
# width = 26,
# height = 3,
# dpi = 1800
# )

ggplot(data = test.aum.dt.combined) +
  geom_point(aes(x = `Test AUM`, y = algorithm, color = seed), size = 5) +
  #ggtitle(c(data.name, cv.type, test.fold)) +
  facet_grid(.~data.name + new.test.fold, scales = "free") +
  theme(panel.spacing=grid::unit(1, "cm"), text = element_text(size=25))

dev.off()

#grid.arrange(aum.list[[1]], aum.list[[2]], ncol = 2)
