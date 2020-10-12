give_test_train_auc = function(
  mat_train, 
  recu_train, 
  mat_test, 
  recu_test, 
  genes,
  resample_n,
  test_train_perma_tune_length
  ) {
  
  sub_mat_train = mat_train
  sub_recu_train = recu_train
  sub_mat_test = mat_test
  sub_recu_test = recu_test
  gene_sig = genes
  print(gene_sig)
  
  x_train = t(varianceStabilizingTransformation(t(sub_mat_train[, gene_sig]), fitType="mean"))
  y_train = factor(sub_recu_train$recuYes, levels = c("Yes", "No"))
  
  fold_n = 3
  repeats_n = 5
  fold_repeats = fold_n*repeats_n
  sampling_type = "down"
  print(paste("fold_n:", fold_n, "repeats_n:", repeats_n, "sampling_type:", sampling_type))
  
  # tune_length = "0.001_2.001_0.1"
  # tune_length = "1_1_1"
  tune_length=test_train_perma_tune_length
  cat("test train give auc have test_train_perma_tune_length ", tune_length, "\n")
  
  if (length(strsplit(tune_length, "_")[[1]]) == 3) {
    model_n = length(
      seq(
        as.numeric(strsplit(tune_length, '_')[[1]][1]), 
        as.numeric(strsplit(tune_length, '_')[[1]][2]), 
        as.numeric(strsplit(tune_length, '_')[[1]][3]))
    )
    
    seeds = vector(mode='list', length = fold_repeats + 1)
    set.seed(123)
    for (i in 1:fold_repeats) seeds[[i]] = sample.int(model_n*model_n*100, model_n*model_n)  # 21*21, 21 is length(seq(0.001, 2.001, 0.1))
    set.seed(123)
    seeds[[fold_repeats+1]] = sample.int(1000, 1)
    
  } else if (length(strsplit(tune_length, "_")[[1]]) == 2) {  # only one model with set alpha and tau.
    model_n=1
    seeds = vector(mode='list', length = fold_repeats + 1)
    set.seed(123)
    for (i in 1:fold_repeats) seeds[[i]] = sample.int(model_n*model_n*100, model_n*model_n)  # 21*21, 21 is length(seq(0.001, 2.001, 0.1))
    set.seed(123)
    seeds[[fold_repeats+1]] = sample.int(1000, 1)
    
  }
  
  
  train_control = trainControl(
    method = "repeatedcv",              
    number = fold_n,       # 3 fold CV.
    repeats = repeats_n,      # 5 repeats.
    sampling = sampling_type,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    seeds = seeds,
    allowParallel = TRUE
  )
  
  cat("Selecting alpha tau on train. \n")
  

  
  set.seed(123)
  perma_train = train(       # inside this train, it uses 10 repeats of 5 cv to do AUC
    x = x_train,
    y = y_train,
    method = perma_caret,
    tuneLength = tune_length,    # seq(0.01, 5.01, 0.5)   # 0.5
    # tuneLength = 0.5,    # seq(0.01, 5.01, 0.5)   # 0.5
    metric = "ROC",
    trControl = train_control
  )
  
  # Save RDS
  # saveRDS(perma_train, "perma_train.RDS")
  # print(perma_train)
  # print(perma_train$bestTune)
  # perma_train$bestTune
  # best_tune_row = as.numeric(rownames(perma_train$bestTune))
  # Best tune resulst.
  # print(perma_train$results[best_tune_row, ])
  print(perma_train$bestTune)
  
  # test
  x_test = t(varianceStabilizingTransformation(t(sub_mat_test[, gene_sig]), fitType="mean"))
  y_test = factor(sub_recu_test$recuYes, levels = c("Yes", "No"))
  
  cat("Use selected alpha tau on train and test.", "\n")
  perma_test = predict(perma_train, x_test, type = 'prob')
  perma_train = predict(perma_train, x_train, type = 'prob')
  
  
  # Convert yes to 1. no to 0.
  y_test_roc = NULL
  y_test_roc[y_test == "No"] = 0
  y_test_roc[y_test == "Yes"] = 1
  
  y_train_roc = NULL
  y_train_roc[y_train == "No"] = 0
  y_train_roc[y_train == "Yes"] = 1
  
  test_auc = roc(y_test_roc, perma_test$Yes)$auc
  library(pROC)
  
  # ROC plot smooth.
  rocTest = smooth(roc(y_test_roc, perma_test$Yes), method="fitdistr")
  test_roc = paste0("Test: ", round(rocTest$auc, 2))
  
  rocTest_df = data.frame(spe=rocTest$specificities, sens=rocTest$sensitivities)
  rocTest_df$type = "Test"
  
  rocTrain = smooth(roc(y_train_roc, perma_train$Yes), method='fitdistr')
  train_roc = paste0("Train: ", round(rocTrain$auc, 2))
  rocTrain_df = data.frame(spe=rocTrain$specificities, sens=rocTrain$sensitivities)
  rocTrain_df$type = "Train"
  
  cat("smoothed AUC test: ", test_roc, " train: ", train_roc, " resample_n ", resample_n, "\n")
  
  roc_tt = rbind(rocTest_df, rocTrain_df)
  roc_plot = ggplot(roc_tt, aes(spe, sens, colour=type)) +
    geom_line(size=1) +
    scale_x_reverse() +
    geom_abline(slope = 1, intercept = 1, linetype=2) +
    labs(x="Specificity", y="Sensitivity") +
    scale_colour_manual(name='AUC', labels=c(test_roc, train_roc), values = c("red", 'blue')) +
    theme_bw(base_size = 22) +
    theme(legend.position = c(0.8, 0.2)) 
    

  # saveRDS(roc_plot, paste0("auc_plot_", resample_n, ".RDS"))
  
  # ggsave(paste0("auc_plot_geneN_", length(gene_sig), "_resample_", resample_n, ".pdf"), roc_plot)
  
  
  # ROC plot not smooth.
  rocTest = roc(y_test_roc, perma_test$Yes)
  test_roc = paste0("Test: ", round(rocTest$auc, 2))
  rocTest_df = data.frame(spe=rocTest$specificities, sens=rocTest$sensitivities)
  rocTest_df$type = "Test"
  
  rocTrain = roc(y_train_roc, perma_train$Yes)
  train_roc = paste0("Train: ", round(rocTrain$auc, 2))
  rocTrain_df = data.frame(spe=rocTrain$specificities, sens=rocTrain$sensitivities)
  rocTrain_df$type = "Train"
  
  cat("not-smoothed AUC test: ", test_roc, " train: ", train_roc, " resample_n ", resample_n, "\n")
  
  
  roc_tt = rbind(rocTest_df, rocTrain_df)
  roc_plot = ggplot(roc_tt, aes(spe, sens, colour=type)) +
    geom_line() +
    scale_x_reverse() +
    geom_abline(slope = 1, intercept = 1, linetype=2) +
    labs(x="Specificity", y="Sensitivity") +
    scale_color_manual(name='AUC', labels=c(test_roc, train_roc), values = c("red", "blue")) +
    theme_bw(base_size = 22) +
    theme(legend.position = c(0.8, 0.2))
    
  

  # print(class(roc_plot))
  # saveRDS(roc_plot, paste0("auc_plot_not_smooth_", resample_n, ".RDS"))
  # ggsave(paste0("auc_plot_not_smooth_", resample_n, ".pdf"), roc_plot)
  
  #
  
  train_auc = roc(y_train_roc, perma_train$Yes)$auc
  # pdf('train_auc_plot.pdf')
  # plot(roc(y_train_roc, perma_train$Yes))
  # dev.off()
  
  return(c(test_auc, train_auc))
  
}
