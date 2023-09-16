KFoldCustomOSCC <- function(MeanMarkerExpression, k = 5, nRepeats = 1,
                            metaData, nFeatures = 2){
  set.seed(1994)
  
  res <- c()
  res_test <- c()
  truth_list <- c()
  predictions_list <- c()
  cur_best <- 0
  cur_best_features <- c()
  
  count = 1
  for(i in 1:nRepeats){
    res_cur <- c()
    res_cur_test <- c()
    cv <- SKM::cv_kfold_strata(as.factor(metaData$Gensini_bin), k = k)
    
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- MeanMarkerExpression[idx, ]
      
      testRes <- colTest(cur_x_train,
                         condition = metaData[match(rownames(cur_x_train),metaData$sample_id), 
                                              "Gensini_bin"])
      x <- rownames(cur_x_train)
      sigFeatures <- testRes$cluster[1:nFeatures]
      cur_x_train <- cur_x_train[, sigFeatures]
      colnames(cur_x_train) <- janitor::make_clean_names(colnames(cur_x_train))
      cur_x_train$response <- factor(metaData[match(x,metaData$sample_id), 
                                                 "Gensini_bin"])
      
      cur_x_train <- na.omit(cur_x_train)
      cur_y_train <- cur_x_train$response
      cur_x_train <- cur_x_train[, 1:ncol(cur_x_train)-1]
      
      cur_x_val <- MeanMarkerExpression[idx_val, sigFeatures]
      x_val <- rownames(cur_x_val)
      cur_x_val$response <- factor(metaData[match(x_val,metaData$sample_id), 
                                            "Gensini_bin"])
      cur_x_val <- na.omit(cur_x_val)
      cur_y_val <- cur_x_val$response
      cur_x_val <- cur_x_val[, 1:ncol(cur_x_val)-1]
      
      model <-  support_vector_machine(
        as.matrix(cur_x_train),
        factor(cur_y_train),
        kernel = "linear",
        verbose = F,
        seed = 1994,
        cost = 1
      )
      
      predictions <- predict(model, as.matrix(cur_x_val))
      predictions_list[[count]] <- predictions$probabilities[, 2]
      
      roc_pred_val <- ROCR::prediction(predictions$probabilities[, 2],
                                       cur_y_val)
      auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
      auc_j <- auc_perf@y.values[[1]]
      
      res_cur[[j]] <- auc_j
      
      if(auc_j >= cur_best){
        cur_best <- auc_j
        cur_best_features <- c(cur_best_features,colnames(cur_x_train))
        message(paste('The best AUC score is: ', round(cur_best, 2)))
      }
    }
    res[[i]] <- mean(unlist(res_cur))
   
  }
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  return(cur_best_features)
}


runKFeatures <- function(data, featureRange = c(2,5,10,15,20,25), clinicalData){
  
  # define a list to hold the plots
  myplots <- vector('list', length = length(featureRange)*2)
  
  for (i in length(featureRange)) {
    # generate top n features
    nFeatures <-featureRange[[i]]
    t <- KFoldCustomOSCC(MeanMarkerExpression = data, k = 5, nRepeats = 20,
                         metaData = clinicalData,
                         nFeatures = nFeatures)
    t <- sort(table(t), decreasing=TRUE)[1:nFeatures] %>% names()
    
    # Generate classification data
    x <- rownames(data)
    colnames(data) <- janitor::make_clean_names(colnames(data))
    
    data_final <- MeanMarkerExpression[, t]
    colnames(data_final) <- janitor::make_clean_names(colnames(data_final))
    
    data_final$response <- factor(clinicalData[match(x,clinicalData$sample_id), 
                                               "Gensini_bin"])
    data_final <- na.omit(data_final)
    
    # build classification model
    n <- ncol(data_final)
    train_x <- data_final[, 1:n-1] %>%
      as.matrix()
    train_y <- data_final$response
    # train_y <- factor(if_else(train_y == 'N', 0, 1))
    
    model_svm <- support_vector_machine(
      as.matrix(train_x),
      factor(train_y),
      kernel = "linear",
      verbose = F,
      seed = 1994)
    
    # do prediction
    predictions <- predict(model_svm, as.matrix(train_x))
    
    roc.pred.test <- ROCR::prediction(predictions$probabilities[, 2], 
                                      train_y)
    perf <- ROCR::performance(roc.pred.test, "tpr", "fpr")
    
    
    plt_dat = data.frame(
      FPR = perf@x.values[[1]],
      TPR = perf@y.values[[1]]
    )
    
    auc_perf <- ROCR::performance(roc.pred.test, measure = "auc")
    auc <- auc_perf@y.values[[1]]
    
    # generate AUC plot
    p.auc <- ggplot(plt_dat, aes(x = FPR, y = TPR)) +
      geom_line(colour = "blue") +
      labs(x = perf@x.name, y = perf@y.name) +
      geom_abline(slope = 1, intercept = 0) + theme_bw() + 
      theme(
        plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5),
        plot.subtitle = element_text(color = "red", size = 16, hjust = 0.5),
        axis.title.x = element_text(color="Black", size=16),
        axis.title.y = element_text(color="Black", size=16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16) #change legend text font size)
      )  + ggtitle(paste("AUC for top ", nFeatures, 'features: ' ,
                         round(auc, digits = 2)))
    # generate model weights
    sv <- model_svm[["fitted_model"]][["SV"]]
    coef <- model_svm[["fitted_model"]][["coefs"]]
    
    w = t(sv) %*% coef %>%
      as.data.frame()
    b = -model_svm[["fitted_model"]][["rho"]]
    
    w$Feature <- rownames(w)
    colnames(w)[[1]] <- 'Weight'
    
    # plot model weights
    w$Feature <- gsub('fox_p3', 'foxp3', w$Feature)
    w$Feature <- gsub('pan_ck', 'panck', w$Feature)
    w$Feature <- stri_replace_last(w$Feature, fixed = "_", ": ") 
    w$Feature <- gsub('_', ' ', w$Feature)
    
    p.weights <- ggplot(w, aes(x=reorder(Feature, Weight), y=Weight, fill=Weight))+
      geom_bar(stat="identity") +
      scale_fill_gradient2(low="darkblue", high="darkred") +
      ggtitle("Plot of feature weights for SVM model")+
      labs(x = "Feature",
           y = "Weight") +
      theme(axis.title = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 15),
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
      coord_flip()
    
    # myplots[[2*i]] <- p.weights
    p <- cowplot::plot_grid(p.auc, p.weights)
  }
  
  return(list(Features = t, Plots = p))
}

getProp <- function(cells, feature = "cellType", 
         imageID = "imageID", logit = TRUE) {
  if (is.data.frame(cells)) {
    df <- cells[, c(imageID, feature)]
  }
  
  if (is(cells, "SingleCellExperiment") | is(cells, "SpatialExperiment")) {
    df <- as.data.frame(SummarizedExperiment::colData(cells))[, c(imageID, feature)]
  }
  
  if (is(cells, "SegmentedCells")) {
    cellSummary <- cellSummary(cells, bind = TRUE)
    df <- as.data.frame(cellSummary[, c(imageID, feature)])
  }
  
  
  tab <- table(df[, imageID], df[, feature])
  tab <- sweep(tab, 1, rowSums(tab), "/")
  
  
  if (logit) {
    try(tab[tab == 0] <- 0.001, silent = TRUE)
    try(tab[tab == 1] <- 0.999, silent = TRUE)
    tab <- gtools::logit(tab)
  }
  as.data.frame.matrix(tab)
}

KFoldCustom <- function(train_x, train_y, k = 5, nRepeats = 1,
                        method = 'glm', alpha = 1){
  set.seed(1994)
  
  res <- c()
  truth_list <- c()
  predictions_list <- c()
  count = 1
  for(i in 1:nRepeats){
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      if(method == 'glm'){
        model <-  glmnet::cv.glmnet(as.matrix(cur_x_train), 
                                    cur_y_train, family = 'binomial',
                                    alpha = 1, nfolds = 5, type.measure = 'auc')
      } else if(method == 'rf') {
        model <-  SKM::random_forest(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994
        )
      } else {
        model <-  SKM::support_vector_machine(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994)
      }
      
      # make predictions
      if(method == 'glm'){
        predictions <- predict(model, as.matrix(cur_x_val), type = 'response', 
                               s = 'lambda.min')
        predictions_list[[count]] <- predictions[, 1]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions[, 1],
                                         cur_y_val)
      } else{
        predictions <- predict(model, as.matrix(cur_x_val))
        predictions_list[[count]] <- predictions$probabilities[, 2]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions$probabilities[, 2],
                                         cur_y_val)
      }
      auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
      auc_j <- auc_perf@y.values[[1]]
      res_cur[[j]] <- auc_j
      
      count = count + 1
    }
    res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  cv_dat <- list(predictions = predictions_list,
                 labels = truth_list)
  
  return(list(cv_res = unlist(res), cv_dat = cv_dat))
}

