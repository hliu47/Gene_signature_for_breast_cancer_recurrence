perma_give_genes = function(
  mat_input,
  recu_input,
  perma_tune_length
) {
  cat("Give perma gene have perma_tune_length ", perma_tune_length, "\n")
  
  sub_mat_perma = mat_input
  sub_recu_perma = recu_input

  ## VST
  mat_vst_perma =t(varianceStabilizingTransformation(t(sub_mat_perma), blind = FALSE, fitType="mean"))       # varianceStabilizingTransformation transformation.

  gene_all = colnames(mat_vst_perma)
  
  best_gene_len = 1
  auc_gene_all = list()
  
  # Control max length of gene to be 98.
  best_gene_len_largest = ifelse(length(gene_all) < 101, length(gene_all)-2, 98)
  
  cat("Largest best_gene_len will be ", best_gene_len_largest, "\n")
  
  while(best_gene_len <= best_gene_len_largest) {
#    cat("\n")
#    cat('best_gene_len is', best_gene_len, '\n', sep = " ")

    if (best_gene_len == 1) gene_list = as.list(gene_all)
#    cat('gene_list length is', length(gene_list), '\n', sep = " ")
    
    auc_table = foreach(i = 1:length(gene_list), .combine = rbind) %dopar% {
      gene_evaluate = gene_list[[i]]
      
      # print(gene_evaluate)
      
      x = as.matrix(mat_vst_perma[, gene_evaluate])
      colnames(x) = gene_evaluate
      dim(x)
      y = sub_recu_perma$recu1
      table(y)
      y[y==0] = 'No'
      y[y==1] = 'Yes'
      table(y)
      y = factor(y)
      
      train_ctrl_rep = 1
      train_ctrl_number = 3
      train_ctrl_rep_num = train_ctrl_rep*train_ctrl_number
      
      seeds = vector(mode='list', length = train_ctrl_rep_num + 1)
      set.seed(123)
      for (i in 1:train_ctrl_rep_num) seeds[[i]] = sample.int(1000, 1)  # only 1 model tunelength=0.5
      set.seed(123)
      seeds[[train_ctrl_rep_num+1]] = sample.int(1000, 1)
      
      train_control = trainControl(
        method = "repeatedcv",
        # method = "cv",
        # method = 'boot',
        # method = 'none',
        sampling = "down",  # rose smote up down
        number = train_ctrl_number,   # 3-fold.  
        repeats = train_ctrl_rep,  # 1 repeat.      # 1 repeat of 3-fold CV.
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        seeds = seeds,
        allowParallel = FALSE
      )
      
      set.seed(123)
      train_gene = train(       
        x = x,
        y = y,
        method = perma_caret,
        # method = "knn",
        tuneLength = perma_tune_length,   # refer to write_perma_into_Caret.R.
        metric = "ROC",
        # metric = "Kappa",
        trControl = train_control   
      ) 
      # cat("\n")
      # print("triain_gene:")
      # print(train_gene)
      
      train_gene_auc = train_gene$results$ROC
      
      if (train_gene_auc < 0.5) train_gene_auc = 1 - train_gene_auc
      
      gene_eval_auc = append(train_gene_auc, gene_evaluate)   # last position is AUC.
      
      return(gene_eval_auc)
    }
    
    auc_gene = auc_table[auc_table[, 1] == max(auc_table[, 1]), ]  # take max AUC.
    
    if (class(auc_gene) != 'matrix') {   
      auc_gene_all[[best_gene_len]] = auc_gene
      gene_not_eval = gene_all[!(gene_all %in% auc_gene[-1])]
      gene_list = lapply(gene_not_eval, append, auc_gene[-1])
      best_gene_len = best_gene_len + 1
    } else if (class(auc_gene) == "matrix") {   
      # It means adding more genes can not increase AUC anymore.
      auc_gene = auc_gene[1, ]   # take first.
      auc_gene_all[[best_gene_len]] = auc_gene
      best_gene_len = best_gene_len_largest + 1  # should stop.
      
    }
    
  }
  

  res_all = data.frame(auc = rep(NA, length(auc_gene_all)), n_gene = rep(NA, length(auc_gene_all)))
  for (i in 1:length(auc_gene_all)) {
    res_all$auc[i] = as.numeric(as.character(auc_gene_all[[i]][1]))
    res_all$n_gene[i] = i
  }
  
  # Desirability
  library(desirability)
  rocVal = res_all$auc
  d_roc = dMax(0.5, 1, 1)  # scale of 1
  d_size = dMin(1, 100, 2) # scale of 2
  # combine the two with a geometric mean
  d_all = dOverall(d_roc, d_size)
  D = predict(d_all, data.frame(rocVal, res_all$n_gene))
  res_all$D = D
  
  # print(res_all)
  
  # Take best auc
  best_gene = res_all[res_all$auc == max(res_all$auc), ]
  if (length(best_gene$auc) > 1) best_gene = best_gene[1, ]  # if same auc take first
  perma_gene_choose_auc = auc_gene_all[[best_gene$n_gene]][-1] 
  
  # Take best D
  best_gene = res_all[res_all$D == max(res_all$D), ]
  if (length(best_gene$D) > 1) best_gene = best_gene[1, ]  # if same D take first.
  perma_gene_choose_d = auc_gene_all[[best_gene$n_gene]][-1]
  
  # print(res_all)
  # print(auc_gene_all)
  
  return(list(auc = perma_gene_choose_auc, D = perma_gene_choose_d))
}


