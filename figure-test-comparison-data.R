OneFold <- function(testFold.path)
{
  #############################################################################
  ############################# Initializations ###############################
  #############################################################################
  
  #needed libraries
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(future.apply)
  library(directlabels)
  
  #Filename and directory initializations
  cv.path <- dirname(dirname(testFold.path))
  folds.csv <- file.path(cv.path, "folds.csv")
  cv.type <- basename(cv.path)
  test.fold <- basename(testFold.path)
  data.dir <- dirname(dirname(cv.path))
  data.name <- basename(data.dir)
  data.list <- list()
  
  #Initialize a generic plot title for later use
  plot.title <- paste0("Data Name = ", data.name, ", cv.type = ", cv.type, ",\n test.fold = ", test.fold)
  
  #Read the "inputs.csv", "outputs.csv", and "evaluations.csv", for each data 
  for(f in c("inputs", "outputs", "evaluation")){
    f.csv.xz <- file.path(data.dir, paste0(f, ".csv.xz"))
    if(file.exists(f.csv.xz)){
      system(paste("xz -dk", f.csv.xz))
    }
    f.csv <- file.path(data.dir, paste0(f, ".csv"))
    f.dt <- data.table::fread(f.csv)
    data.list[[f]] <- f.dt
  }
  
  
  ## replace positive fp/fn at end with 0 to avoid AUM=Inf.
  data.list[["evaluation"]][min.log.lambda==-Inf & 0<fn, fn := 0]
  data.list[["evaluation"]][max.log.lambda==Inf & 0<fp, fp := 0]
  
  #random cross validation algorithm will randomly assign
  #subtrain/validation folds, where preset will use the test
  #folds already in place.
  cv.algos <- c("random")
  cv.dt.list <- list()
  
  #############################################################################
  ############################# Model Fitting #################################
  #############################################################################
  
  for( curr.cv.algo in cv.algos )
  {
    ## read folds.csv for the specific data type
    folds.dt <- data.table::fread(folds.csv)
    
    #set test folds 
    folds.dt[fold == test.fold, set := "test"]
    folds.dt[fold != test.fold, set := "subtrain"]
    
    #If the current cross validation algorithm is "preset",
    #use the test fold assignments already in place to
    #create subtrain/validation folds
    if( curr.cv.algo == "preset")
    {
      cv.folds <- unique(folds.dt[fold != test.fold]$fold)
    }
    
    #If the current cross validation algorithm is "random",
    #create my own randomized subtrain/validation folds.
    if( curr.cv.algo == "random")
    {
      #initialize my own validation folds
      set.seed(1)
      n.val.folds <- 2
      
      cv.folds <- 1:n.val.folds
      val.fold.assignments <- sample(rep(cv.folds, l = nrow(folds.dt[set == "subtrain"])))
    }
    
    #Loop through every validation folds in the total folds in
    #the cross validation algorithm
    for(val.fold in cv.folds)
    {
      
      #Adjust the folds with respect to the current validation fold
      #and the current cross validation algorithm
      if( curr.cv.algo == "preset")
      {
        folds.dt[set != "test",] <- folds.dt[set != "test",] %>%
          mutate(set = ifelse(fold == val.fold, "validation", "subtrain"))
      }
      if( curr.cv.algo == "random")
      {
        folds.dt[fold == test.fold, set := "test"]
        folds.dt[fold != test.fold, set := "subtrain"]
        folds.dt[set == "subtrain"]$set <- val.fold.assignments
        folds.dt[set != "test",] <- folds.dt[set != "test",] %>% 
          mutate(set = ifelse(set == val.fold, "validation", "subtrain"))
      }
      
      seqs.train <- folds.dt[["sequenceID"]]
      X.all <- scale(data.list$inputs[, -1])
      rownames(X.all) <- data.list$inputs$sequenceID
      X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
      set.list <- list()
      for(s in unique(folds.dt$set)){
        set.list[[s]] <- rownames(X.finite) %in% folds.dt[s==set, sequenceID]
      }
      X.list <- lapply(set.list, function(i)X.finite[i, ])
      neg.t.X.subtrain <- -t(X.list[["subtrain"]])
      y.train <- data.list[["outputs"]][
        seqs.train,
        cbind(min.log.lambda, max.log.lambda),
        on="sequenceID"]
      keep <- apply(is.finite(y.train), 1, any)
      X.train <- X.finite[seqs.train, ]
      init.fun.list <- list(
        IntervalRegressionCV=function(){
          fit <- penaltyLearning::IntervalRegressionCV(
            X.train[keep, ],
            y.train[keep, ])  
          fit[["param.mat"]]
        }
      )
      
      #Set hyperparameters for the training algorithm
      n.seeds <- 2
      num.iterations <- 2
      
      #Fit model. Output is a data.table that includes information
      #relating to the AUC, AUM, and more with respect to every seed,
      #iteration, cross validation algorithm, current validation fold, and
      #fold type (subtrain, validation, test) and more.
      fit.model <- function(X.train, y.train, set.list, num.iterations, seed)
      {
        iteration.dt.list <- list()
        for (init.name in names(init.fun.list))
        {
          init.fun <- init.fun.list[[init.name]]
          set.seed(seed)
          int.weights <- init.fun()
          
          weight.vec <- int.weights[-1]
          intercept <- int.weights[1]
          computeAUM <- function(w, i, is.set) {
            pred.pen.vec <- (X.finite %*% w) + i
            pred.dt <- data.table(
              sequenceID = rownames(pred.pen.vec),
              pred.log.lambda = as.numeric(pred.pen.vec)
            )
            set.dt <- pred.dt[is.set]
            penaltyLearning::ROChange(data.list$evaluation, set.dt, "sequenceID")
          }
          for (iteration in 1:num.iterations) {
            summary.dt.list <- list()
            set.roc.list <- list()
            for (set in names(set.list)) {
              set.roc.list[[set]] <-
                computeAUM(weight.vec, intercept, set.list[[set]])
              summary.dt.list[[set]] <-
                with(set.roc.list[[set]],
                     data.table(set,
                                thresholds[threshold == "predicted"],
                                auc,
                                aum))
            }
            summary.dt <- do.call(rbind, summary.dt.list)
            iteration.dt.list[[paste(seed, init.name, iteration)]] <-
              data.table(seed = paste(seed), init.name, iteration, summary.dt)
            cat(
              sprintf(
                "it=%d seed=%d init=%s cv.algo=%s\n",
                iteration,
                seed,
                init.name,
                curr.cv.algo
              )
            )
            g.dt <- set.roc.list[["subtrain"]][["aum.grad"]]
            ## If aum.grad has some problems with no changes in error then
            ## they may be missing.
            g.vec <- rep(0, ncol(neg.t.X.subtrain))
            names(g.vec) <- colnames(neg.t.X.subtrain)
            g.vec[g.dt[["sequenceID"]]] <- g.dt[["lo"]]
            is.differentiable <- all(g.dt[["lo"]] == g.dt[["hi"]])
            direction.vec <- neg.t.X.subtrain %*% g.vec
            take.step <- function(s) {
              weight.vec + s * direction.vec
            }
            
            #line search
            set.aum.list <- list()
            for (step.size in 10 ^ seq(-10, 0, by = 0.5)) {
              new.weight.vec <- take.step(step.size)
              for (set in "subtrain") {
                set.roc <- computeAUM(new.weight.vec, 0, set.list[[set]])
                set.aum.list[[paste(step.size, set)]] <-
                  data.table(step.size,
                             set,
                             aum = set.roc$aum,
                             intercept = set.roc$thresholds[threshold == "min.error", (max.thresh +
                                                                                         min.thresh) / 2])
              }#line search
            }
            set.aum <- do.call(rbind, set.aum.list)
            best.dt <- set.aum[, .SD[min(aum) == aum], by = set]
            weight.vec <- take.step(best.dt[["step.size"]])
            intercept <- best.dt[["intercept"]]
          }#iteration
        }
        output.dt <- data.table(do.call(rbind, iteration.dt.list),
                                data.name,
                                cv.type,
                                test.fold,
                                is.differentiable)
        
        output.dt
        
      }
      
      for(curr.seed in 1:n.seeds)
      {
        #Create the data.tables used during the validation phase
        cv.dt.list[[paste(curr.seed, val.fold, curr.cv.algo)]] <- data.table( fit.model(X.train, 
                                                                                        y.train, 
                                                                                        set.list,
                                                                                        num.iterations = num.iterations, 
                                                                                        seed = curr.seed),
                                                                              val.fold = paste(val.fold),
                                                                              cv.algo = curr.cv.algo)
        
      }#seed
      cat(sprintf("Validation fold %d complete\n", val.fold))
    }#validation fold
    
  }
  
  #Create the cross validation data.table and save it to a file
  cv.dt <- do.call(rbind, cv.dt.list)
  cv.csv <- file.path(testFold.path, "linear-model-cv-aum.csv")
  data.table::fwrite(cv.dt, cv.csv)
  
  #############################################################################
  ################################## Plotting #################################
  #############################################################################
  
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
  
  
  for(curr.cv.algo in cv.algos)
  {
    for(curr.set in unique(cv.dt$set))
    {
      #Create an auc line graph for every dataset, including a mean line graph
      out.png.path <- file.path(testFold.path, paste0("linear-model-precv-", curr.set, 
                                                      "-", curr.cv.algo, "-aum-line-graph.png"))
      png(out.png.path)
      print(ggplot(data = cv.dt[set == curr.set & cv.algo == curr.cv.algo]) +
              geom_line(aes(x = iteration, y = aum, color = val.fold)) +
              geom_point(data = min.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = min.aum.iteration, y = min.aum, color = val.fold)) +
              geom_line(data = mean.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = iteration, y = mean.aum, color = "mean"), size = 1) +
              geom_point(data = min.mean.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = min.aum.iteration, y = min.aum), size = 3) +
              facet_wrap(.~seed, labeller = label_both) +
              ggtitle(paste0(curr.set, " AUM with for each validation fold\n", plot.title)))
      dev.off()
      
      
      #Create an auc line graph for every dataset, including a mean line graph
      out.png.path <- file.path(testFold.path, paste0("linear-model-precv-", curr.set,
                                                      "-", curr.cv.algo, "-auc-line-graph.png"))
      png(out.png.path)
      print(ggplot(data = cv.dt[set == curr.set & cv.algo == curr.cv.algo]) +
              geom_line(aes(x = iteration, y = auc, color = val.fold)) +
              geom_point(data = min.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = max.auc.iteration, y = max.auc, color = val.fold)) +
              geom_line(data = mean.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = iteration, y = mean.auc, color = "mean"), size = 1) +
              geom_point(data = min.mean.cv.dt[set == curr.set & cv.algo == curr.cv.algo], mapping = aes(x = max.auc.iteration, y = max.auc), size = 3) +
              facet_wrap(.~seed, labeller = label_both) +
              ggtitle(paste0(curr.set, " AUC with for each validation fold\n", plot.title)))
      dev.off()
    }
  }
  
  
  #############################################################################
  ############################# Post CV Initialization ########################
  #############################################################################
  
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set := "subtrain"]
  iteration.dt.list <- list()
  
  seqs.train <- folds.dt[["sequenceID"]]
  X.all <- scale(data.list$inputs[, -1])
  rownames(X.all) <- data.list$inputs$sequenceID
  X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
  set.list <- list()
  for(s in unique(folds.dt$set)){
    set.list[[s]] <- rownames(X.finite) %in% folds.dt[s==set, sequenceID]
  }
  X.list <- lapply(set.list, function(i)X.finite[i, ])
  neg.t.X.subtrain <- -t(X.list[["subtrain"]])
  y.train <- data.list[["outputs"]][
    seqs.train,
    cbind(min.log.lambda, max.log.lambda),
    on="sequenceID"]
  keep <- apply(is.finite(y.train), 1, any)
  X.train <- X.finite[seqs.train, ]
  init.fun.list <- list(
    IntervalRegressionCV=function(){
      fit <- penaltyLearning::IntervalRegressionCV(
        X.train[keep, ],
        y.train[keep, ])  
      fit[["param.mat"]]
    }
  )
  
  #############################################################################
  ############################## Post CV Fitting ##############################
  #############################################################################
  
  post.cv.dt.list <- list()
  for(curr.seed in 1:n.seeds)
  {
    post.cv.dt.list[[paste(curr.seed)]] <- fit.model(X.train, 
                                                     y.train,
                                                     set.list,
                                                     num.iterations = num.iterations,
                                                     seed = curr.seed)
  }#seed
  
  #Create post cross validation data table
  post.cv.dt <- do.call(rbind, post.cv.dt.list)
  post.cv.csv <- file.path(testFold.path, "linear-model-post-cv-aum.csv")
  data.table::fwrite(post.cv.dt, post.cv.csv)
  
  #Create data.table with best linear model results with respect to
  #the test set
  best.linear.dt <- post.cv.dt[set == "test", .(aum = min(aum),
                                                auc = max(auc),
                                                type = "best.linear",
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
      
      selected.aum.dt.list[[paste0(curr.seed, curr.cv.algo)]] <- data.table(aum = curr.dt$aum[min.aum.iteration],
                                                                            auc = curr.dt$auc[min.aum.iteration],
                                                                            type = paste0("selected.", curr.cv.algo,".aum"),
                                                                            seed = curr.seed,
                                                                            aum.iteration = min.aum.iteration,
                                                                            auc.iteration = min.aum.iteration)
      
      selected.auc.dt.list[[paste0(curr.seed, curr.cv.algo)]] <- data.table(aum = curr.dt$aum[max.auc.iteration],
                                                                            auc = curr.dt$auc[max.auc.iteration],
                                                                            type = paste0("selected.", curr.cv.algo, ".auc"),
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
    post.cv.dt[set == "test", .(aum = first(aum),
                                auc = first(auc),
                                type = "initial",
                                aum.iteration = 1,
                                auc.iteration = 1),  
               by = seed]
  
  #Combine intitial, selected, and best.linear results for plotting purposes
  test.aum.dt <- rbind(best.linear.dt, selected.aum.dt, selected.auc.dt, initial.dt)
  
  
  #Create a dotplot with every cv algorithm,
  out.png.path <- file.path(testFold.path, "linear-model-test-aum-comparison.png")
  png(out.png.path)
  print(ggplot(data = test.aum.dt) +
          geom_point(aes(x = aum, y = type, color = seed)) +
          ggtitle(plot.title))
  
  dev.off()
  
  out.png.path <- file.path(testFold.path, "linear-model-test-auc-comparison.png")
  png(out.png.path)
  
  print(ggplot(data = test.aum.dt) +
          geom_point(aes(x = auc, y = type, color = seed)) +
          ggtitle(plot.title))
  
  dev.off()
  
  
  out.png.path <- file.path(testFold.path, "linear-model-postcv-subtrain-aum-line-graph.png")
  png(out.png.path)
  
  print(ggplot(post.cv.dt[set == "subtrain",]) +
          geom_line(aes(x = iteration, y = aum)) +
          facet_wrap(.~seed, labeller = label_both) +
          ggtitle(paste0("subtrain AUM for every seed after CV"), plot.title))
  
  dev.off()
  
  
  out.png.path <- file.path(testFold.path, "linear-model-postcv-subtrain-auc-line-graph.png")
  png(out.png.path)
  
  print(ggplot(post.cv.dt[set == "subtrain",]) +
          geom_line(aes(x = iteration, y = auc)) +
          facet_wrap(.~seed, labeller = label_both) +
          ggtitle(paste0("subtrain AUC for every seed after CV"), plot.title))
  
  dev.off()
  
  
  out.png.path <- file.path(testFold.path, "linear-model-postcv-test-aum-line-graph.png")
  png(out.png.path)
  
  print(ggplot(post.cv.dt[set == "test",]) +
          geom_line(aes(x = iteration, y = aum)) +
          geom_point(data = best.linear.dt, mapping = aes(x = aum.iteration, y = aum, color = "best.linear")) +
          # geom_point(data = selected.aum.dt[type == "selected.preset.aum",], mapping = aes(x = aum.iteration, y = aum, color = "selected.preset.aum")) +
          # geom_point(data = selected.auc.dt[type == "selected.preset.auc",], mapping = aes(x = aum.iteration, y = aum, color = "selected.preset.auc")) +
          geom_point(data = selected.aum.dt[type == "selected.random.aum",], mapping = aes(x = aum.iteration, y = aum, color = "selected.random.aum")) +
          geom_point(data = selected.auc.dt[type == "selected.random.auc",], mapping = aes(x = aum.iteration, y = aum, color = "selected.random.auc")) +
          facet_wrap(.~seed, labeller = label_both) +
          ggtitle(paste0("test AUM for every seed after CV"), plot.title))
  
  dev.off()
  
  out.png.path <- file.path(testFold.path, "linear-model-postcv-test-auc-line-graph.png")
  png(out.png.path)
  
  print(ggplot(post.cv.dt[set == "test",]) +
          geom_line(aes(x = iteration, y = auc)) +
          geom_point(data = best.linear.dt, mapping = aes(x = aum.iteration, y = auc, color = "best.linear")) +
          # geom_point(data = selected.aum.dt[type == "selected.preset.aum",], mapping = aes(x = aum.iteration, y = auc, color = "selected.preset.aum")) +
          # geom_point(data = selected.auc.dt[type == "selected.preset.auc",], mapping = aes(x = aum.iteration, y = auc, color = "selected.preset.auc")) +
          geom_point(data = selected.aum.dt[type == "selected.random.aum",], mapping = aes(x = aum.iteration, y = auc, color = "selected.random.aum")) +
          geom_point(data = selected.auc.dt[type == "selected.random.auc",], mapping = aes(x = aum.iteration, y = auc, color = "selected.random.auc")) +
          facet_wrap(.~seed, labeller = label_both) +
          ggtitle(paste0("test AUM for every seed after CV"), plot.title))
  
  dev.off()
  
}