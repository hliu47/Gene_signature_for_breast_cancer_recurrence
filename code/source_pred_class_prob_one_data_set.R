library(caret)
library(ggplot2)
library(mlbench)
library(DESeq2)
library(pROC)
library(e1071)
library(ROCR)
library(survival)
library(glmnet)
library(txtplot)


source("./source_loocv_give_prob.R")


pred_class_prob_one_data_set = function(gene_type, data_type, data_set_mat, data_set_recu, pam50_subtype, node_status) {
  
  my_meth = "together"
  my_meth_by = "AUC"
  
  # subtype and node_status.
  sam_sele = data_set_recu[, pam50_subtype] & data_set_recu[, "nodePos"] == node_status & !is.na(data_set_recu[, 'nodePos'])
  
  print(table(sam_sele, useNA="always"))
  
  sub_mat = data_set_mat[sam_sele, ]
  sub_recu = data_set_recu[sam_sele, ]
  
  print(table(sub_recu$recu1))
  
  # use genes from together_best_genes
  my_genes = readRDS(paste0("../1_together_deseq_cox_perma/together_best_genes_max_by", 
                            my_meth_by, "_test_", pam50_subtype, "_node_", node_status, "_.RDS"))
  print(my_genes)
  
  # Use other gene set.
  # Bao et al. J Transl Med (2019) 17:380 https://doi.org/10.1186/s12967-019-2126-6
  if (gene_type == "baoGenes") {
    bao_genes = c("ITPRIPL1", "SIAH2", "KCNH8", "KRT19", "NDRG2", "STAC2", "TPD52", "EZR", "PCDHGA12", "HIF3A", "PCDHGA3", "ECRG4", "CCND2")   # C2orf40 -> ECRG4
    my_genes = bao_genes
  } else if (gene_type == "randomGenes") {
    my_genes = colnames(sub_mat)[100:106]
  } 
  
  cat('length of my_genes is', length(my_genes), '\n')
  cat(my_genes, '\n')
  
  sub_mat = sub_mat[, my_genes]
  
  cat('sub_mat dim is ', "\n")
  print(dim(sub_mat))
  cat('sub_recu has recu table ', "\n")
  print(table(sub_recu$recu1))
  
  ## VST
  # mat_vst = log(sub_mat + 1)      # log transformation.
  if (data_type %in% c("tcga", 'gh')) {
    sub_mat_vst =t(varianceStabilizingTransformation(t(sub_mat), fitType = 'mean', blind = FALSE))       # varianceStabilizingTransformation transformation
  } else {
    sub_mat_vst = sub_mat
  }
  
  #######    LOOCV for each sample to give prob then use all sample prob to get all sample class.
  #######    Set one alpha tau in within train, then test a sequence of alpha tau after LOOCV for all sample.
  
  gene_sig = my_genes
  
  ###############################################
  ### LOOCV to give each sample prob.
  ###############################################
  
  # ####### Test which alpha tau is the bset.
  # ####### Uncomment below.
  # res_all_many_tune_length = NULL
  # for (tune_length in seq(0.001, 2.1, 0.1)) {    # 0.6 is the best.
  #   tryCatch({
  #     # print(tune_length)
  #     res_all = loocv_give_prob(sub_mat_vst, sub_recu, gene_sig, tune_length)
  #     res_all_pred_prob_roc = roc(res_all$true_class, res_all$pred_prob_Yes)
  #     res_all_pred_prob_auc = as.character(res_all_pred_prob_roc$auc)
  #     # print(res_all_pred_prob_auc)
  #     res_all_many_tune_length = rbind(res_all_many_tune_length, c(tune_length, res_all_pred_prob_auc))
  #   }, error=function(err){
  #     print(paste("errorrrrr:", err, "for", tune_length))
  #   })
  # }
  # 
  # colnames(res_all_many_tune_length) = c("Tune_length_alpha_tau", "AUC")
  # print(res_all_many_tune_length)   # 0.6 is the best.
  # saveRDS(res_all_many_tune_length, paste0("res_all_many_tune_length", data_type, '.RDS'))
  # print("Test tune_length DONE.")
  # ######
  
  
  # After above test, select the best alpha tau input as tune_length.
  best_tune_length = 0.6  # g5 is 0.6.
  res_all = loocv_give_prob(sub_mat_vst, sub_recu, gene_sig, tune_length=best_tune_length)   
  print(roc(res_all$true_class, res_all$pred_prob_Yes))
  
  saveRDS(res_all, paste("pred_class", pam50_subtype, my_meth, "by", 
                         my_meth_by, "dataType", data_type, 
                         "geneType", gene_type, 'node_status', node_status, ".RDS", sep = "_"))
  
  
}








