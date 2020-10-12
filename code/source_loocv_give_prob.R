# LOOCV to give each sample prob.

source('../0.0.1_perma_sources/Helper_functions.R')   # calculate kappa.
# source('../0.0.1_perma_sources/rf_caret_my.R')
source("../0.0.1_perma_sources/write_perma_into_Caret.R")
source("../0.0.1_perma_sources/perma_generic.R")
source("../0.0.1_perma_sources/perma_give_prob_myOwn_2.R")


# in parallel.
library(doMC)
library(doParallel)
library(foreach)
n_cores = parallel::detectCores()
registerDoMC(n_cores - 2)
#registerDoMC(1)


loocv_give_prob = function(sub_mat_vst, sub_recu, gene_sig, tune_length) {
  
  res_all = foreach(sample_i = 1:dim(sub_mat_vst)[1], .combine = rbind, .verbose = FALSE) %dopar% {
    # without sample_i is for train
    # print(sample_i)
    x = as.matrix(sub_mat_vst[-sample_i, gene_sig])
    colnames(x) = gene_sig
    
    y = sub_recu$recu1[-sample_i]
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
      # method = "LOOCV",       # LOOCV has no resampling, it just run once.
      method = "repeatedcv",
      number = train_ctrl_number,
      repeats = train_ctrl_rep,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      savePredictions = "all",
      sampling = "down",
      seeds = seeds,
      allowParallel = FALSE
    )
    
    set.seed(123456)
    within_train = train(
      x = x,
      y = y,
      method = perma_caret,
      tuneLength = tune_length,
      # tuneLength = "0.001_2.2_0.2",     # seq(0.001, 2.5, 0.5)
      metric = "ROC",
      trControl = train_control
    )

    within_train_roc = roc(within_train$pred$obs, within_train$pred$Yes, quiet = TRUE)

    best_auc_threshold = coords(
      within_train_roc,
      "best",
      best.method = "closest.topleft",
      transpose = FALSE,
      ret = 'threshold'
    )

    best_auc_direction = within_train_roc$direction

    # predict sample_i
    x_i = matrix(sub_mat_vst[sample_i, gene_sig], nrow = 1)
    colnames(x_i) = gene_sig

    pred_sample_i = predict(within_train, x_i, type = 'prob')

    if (best_auc_direction == '<') {
      pred_class_i = ifelse(pred_sample_i$Yes >= best_auc_threshold, "Yes", "No")
    } else {
      pred_class_i = ifelse(pred_sample_i$Yes <= best_auc_threshold, "Yes", "No")
    }
    
    res_i = data.frame(
      sample = rownames(sub_mat_vst)[sample_i],
      pred_class_within_train = as.character(pred_class_i),
      true_class = sub_recu$recu1[sample_i],
      years = sub_recu$days[sample_i]/365.25,
      pred_prob_Yes = pred_sample_i$Yes,
      withinTrain_repeatedcv_AUC = within_train_roc$auc,
      withinTrain_ROC_threshold = as.character(best_auc_threshold),
      withinTrain_ROC_direc = best_auc_direction
    ) 
    
    return(res_i)
    
  }
  
  
  return(res_all)
  
}

