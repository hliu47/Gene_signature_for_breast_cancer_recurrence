source("./sources/give_perma_genes_together.R")
source("./sources/combine_perm_gene_and_test_train_auc.R")
source("./sources/deseq_give_genes.R")
source("./sources/cox_give_genes.R")
source("./sources/perma_give_genes.R")
source("./sources/give_test_auc.R")

source('../0.0.1_perma_sources/Helper_functions.R')   # calculate kappa.
source("../0.0.1_perma_sources/write_perma_into_Caret.R")
source("../0.0.1_perma_sources/perma_generic.R")
source("../0.0.1_perma_sources/perma_give_prob_myOwn_2.R")


together_deseq_cox_perma_subtype = function(subtype, rep_train_test_n, run_one_or_many, rep_start) {
  
  
  pam50_subtype = subtype
  
  sam_sele = tcga_recu[, pam50_subtype]
  sub_mat = tcga_mat[sam_sele, ]
  sub_recu = tcga_recu[sam_sele, ]
  
  # print(dim(sub_mat))
  
  # for test on mac, on aws it will use all.
  # if (grepl('darwin', version$os)) sub_mat = sub_mat[, 1:1000]

  # repeated_train_test_n = ifelse(grepl("darwin", version$os), rep_train_test_n, 100)
  repeated_train_test_n = rep_train_test_n
  # train_ratio = 3/4  # 3/4 for train, 1/4 for test
  train_ratio = 2/3  # 3/4 for train, 1/4 for test
  cat("train_ratio is ", train_ratio, "\n")
  
  set.seed(1234)
  train_index = createDataPartition(
    sub_recu$recuYes, 
    times = repeated_train_test_n,
    p = train_ratio      # 3/4 for train, 1/4 for test.
  )
  # for test
  
  
  res_gene_auc_all = list()
  if (run_one_or_many == "one") rep_seq = repeated_train_test_n 
  if (run_one_or_many == "many") rep_seq = seq(rep_start, repeated_train_test_n, 1) 
  #if (run_one_or_many == "some") rep_seq = repeated_train_test_n   # c(3,5,7)
  # print(rep_seq)  
  for (resample_n in rep_seq) {
    # resample_n = 5
    print(resample_n)
    tryCatch({
      # resample_n = 1
      cat("\n")
      cat("\n")
      cat("resample is on ", resample_n, "\n")
      
      # Train
      sub_mat_train = sub_mat[train_index[[resample_n]], ]
      sub_recu_train = sub_recu[train_index[[resample_n]], ]
      print(dim(sub_mat_train))
      print(table(sub_recu_train$recuYes, useNA='always'))
      
      # Test.
      sub_mat_test = sub_mat[-train_index[[resample_n]], ]
      sub_recu_test = sub_recu[-train_index[[resample_n]], ]
      print(dim(sub_mat_test))
      print(table(sub_recu_test$recuYes, useNA="always"))
      
      cat("Producing gene_sig ...", "\n")
    
      deseq_p = 0.05  # in DESeq2, padj
      deseq_f = 2     # in DESeq2, log2 fold change.  # for lum_a 1, for basal 2.
      #perma_tune_length = -1.8  # can only be one number. In selecting perma gene.
      # test_train_perma_tune_length="1_1_1"  # in format of 1_1_1   After perma gene selected, calculate AUC on train and test.
      # test_train_perma_tune_length_seq = c("-2_2_1")  # both_alpha_tau: start_end_by
      # test_train_perma_tune_length_seq = c("-2_-10")  # alpha_tau

      
      #################### For selecting perma gene. #########################################
      perma_tune_length_start=-10   # for both alpha tau.
      perma_tune_length_end=10
      perma_tune_length_by=1
      
      perma_tune_length_grid=expand.grid(
        seq(perma_tune_length_start,
            perma_tune_length_end,
            perma_tune_length_by), 
        seq(perma_tune_length_start,
            perma_tune_length_end,
            perma_tune_length_by))
      perma_tune_length_seq = apply(perma_tune_length_grid, 1, function(x) paste(x[1], x[2], sep = "_"))
      
      # Specify alpha tau
      #perma_tune_length_seq = c("1_1")  #  alpha_tau
      #######################################################################################

      ################## For computing AUC on test and train. ###############################
      test_train_perma_tune_length_start=-2   # for both alpha tau.
      test_train_perma_tune_length_end=2
      test_train_perma_tune_length_by=1
      
      test_train_perma_tune_length_grid=expand.grid(
        seq(test_train_perma_tune_length_start,
            test_train_perma_tune_length_end,
            test_train_perma_tune_length_by), 
        seq(test_train_perma_tune_length_start,
            test_train_perma_tune_length_end,
            test_train_perma_tune_length_by))
      test_train_perma_tune_length_seq = apply(test_train_perma_tune_length_grid, 1, function(x) paste(x[1], x[2], sep = "_"))
      
      # Specify alpha tau
      test_train_perma_tune_length_seq = "1_1"  #  alpha_tau
      #######################################################################################
      
      perma_tunLen_all = list()
      for (perma_tune_length in perma_tune_length_seq) {   # In selecting perma genes.
        for (test_train_perma_tune_length in test_train_perma_tune_length_seq) {
          tryCatch({
            
            cat('perma_tune_length ', perma_tune_length, '\n')
            combine_res = combine_perma_gene_and_test_train_auc(
              sub_mat_train, 
              sub_recu_train,
              sub_mat_test,
              sub_recu_test,
              deseq_p,
              deseq_f,
              perma_tune_length,
              test_train_perma_tune_length,
              resample_n)
            
            test_train_auc = combine_res$test_train_auc
            gene_sig = combine_res$gene_sig
            
            res_gene_auc_byAUC = c(test_train_auc, gene_sig$auc)     # test_train_auc first is test second is train.
            
            print(res_gene_auc_byAUC)
            perma_tunLen_all[[as.character(paste(perma_tune_length, test_train_perma_tune_length, sep = "-"))]] = list(res=res_gene_auc_byAUC)
            cat("\n")
            
            #print(perma_tunLen_all)
            #tunLen = "-2_2_0.01"  # node Pos lum_a -1.8 
            #saveRDS(perma_tunLen_all, paste0('perma_tunLen_all_tunLen_', tunLen, ".RDS"))
            
          }, error=function(e){
            cat('error in perma_tune_length', conditionMessage(e), '\n')
          })
        }
      }
      
      
      # After get best result use thie below.
      res_gene_auc_all[[as.character(resample_n)]] = list(byAUC = res_gene_auc_byAUC)
      
      # print(res_gene_auc_all)
      
      # save it every resample.
#       saveRDS(res_gene_auc_all, 
#               paste('res_gene_auc_all', pam50_subtype,
# 		    'rep_start', rep_start, 
#                     'repeatTrainTest', repeated_train_test_n,
#                     'trainRatio', format(train_ratio, digits = 2), 
#                     ".RDSrds", sep = "_"))
      # garbage collection free up memory.
      gc()

      # save to s3 every resample. prevent ec2 shutdown by aws.
      # if (!grepl('darwin', version$os))  system("aws s3 cp . s3://hliulab/together_perma_gene/ --recursive")
      
      
    }, error = function(e){
      cat("ErrorForOneResampleDoesntMatter:", conditionMessage(e), "\n")
    })
    
    
  }
  
  
  
}
