combine_perma_gene_and_test_train_auc = function(
  sub_mat_train, 
  sub_recu_train,
  sub_mat_test,
  sub_recu_test,
  deseq_p,
  deseq_f,
  perma_tune_length,
  test_train_perma_tune_length,
  resample_n) {
  
  
  gene_sig = give_perma_genes_together(
    sub_mat = sub_mat_train,
    sub_recu = sub_recu_train,
    resample_n = resample_n,
    deseq_p=deseq_p,
    deseq_f=deseq_f,
    perma_tune_length=perma_tune_length
  )
  cat("gene_sig done", "\n")
  
  ############  use perma_genes on train then predict on test.
  cat("Using selected gene_sig on train and test to get AUC", "\n")
  
  # gene_sig by auc
  test_train_auc = give_test_train_auc(
    mat_train = sub_mat_train,
    recu_train = sub_recu_train,
    mat_test = sub_mat_test,
    recu_test = sub_recu_test,
    genes = gene_sig$auc,
    resample_n = resample_n,
    test_train_perma_tune_length = test_train_perma_tune_length
  )
  # print("Test done got test AUC.")
  return(list(gene_sig=gene_sig,
              test_train_auc=test_train_auc))
}



