library(dplyr)
library(ggplot2)
library(data.table)
source("figure-test-comparison-data.R")

#Set the neuroblastoma dataset here
nb.data.dir <- file.path("neuroblastoma-data/data")

data.dir.vec <- c(file.path("ATAC_JV_adipose/cv/equal_labels/testFolds/4"),
                  file.path("H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/2"),
                  file.path("systematic/cv/R-3.6.0-profileSize/testFolds/1"),
                  file.path("H3K4me3_XJ_immune/cv/equal_labels/testFolds/2"),
                  file.path("H3K4me3_XJ_immune/cv/equal_labels/testFolds/4"))

testFold.vec <- sapply(data.dir.vec, function(x){file.path(nb.data.dir, x)})

test.aum.dt.list <- list()

if(!file.exists("figure-test-comparison.csv"))
{
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
    
    if(!file.exists(file.path(testFold.path, "linear-model-cv-aum.csv")) ||
       !file.exists(file.path(testFold.path, "linear-model-post-cv-aum.csv")))
    {
      OneFold(testFold.path)
    }
    
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
    for(curr.cv.algo in cv.algos)
    {
      for(curr.seed in 1:n.seeds)
      {
          curr.dt <- post.cv.dt[seed == curr.seed & set == "test"]
          min.aum.iteration <- selected.iter.dt[seed == curr.seed & cv.algo == curr.cv.algo]$min.aum.iteration
          max.auc.iteration <- selected.iter.dt[seed == curr.seed & cv.algo == curr.cv.algo]$max.auc.iteration
          
        for(aum.auc in c("aum", "auc"))
        {
          if(aum.auc == "aum")
          {
            selected.aum.dt.list[[paste0(curr.seed, curr.cv.algo, aum.auc)]] <- data.table(`Test AUM` = curr.dt$aum[min.aum.iteration],
                                                                                           `Test AUC` = curr.dt$auc[min.aum.iteration],
                                                                                           algorithm = paste0("Min.Valid.AUM"),
                                                                                           seed = curr.seed,
                                                                                           aum.iteration = min.aum.iteration,
                                                                                           auc.iteration = min.aum.iteration)
          } else {
            selected.aum.dt.list[[paste0(curr.seed, curr.cv.algo, aum.auc)]] <- data.table(`Test AUM` = curr.dt$aum[max.auc.iteration],
                                                                                           `Test AUC` = curr.dt$auc[max.auc.iteration],
                                                                                           algorithm = paste0("Max.Valid.AUC"),
                                                                                           seed = curr.seed,
                                                                                           aum.iteration = max.auc.iteration,
                                                                                           auc.iteration = max.auc.iteration)
          }
          
          
          
        }
      }
    }
    
    selected.aum.dt <- do.call(rbind, selected.aum.dt.list)
    
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
    test.aum.dt <- rbind(selected.aum.dt, initial.dt)
    
    test.aum.dt.list[[curr.index]] <- test.aum.dt %>% mutate(data.name = data.name, 
                                                             cv.type = cv.type,
                                                             test.fold = test.fold)
  }
  
  test.aum.dt.combined <- do.call(rbind, test.aum.dt.list) %>%
    mutate(new.test.fold = paste0("Test Fold ", test.fold))
  
  data.table::fwrite(test.aum.dt.combined, "figure-test-comparison.csv")
} else {
  test.aum.dt.combined <- data.table::fread("figure-test-comparison.csv")
}

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

