source('./sources/give_heatmap.R')

give_perma_genes_together = function(
  sub_mat,
  sub_recu,
  resample_n,
  deseq_p,
  deseq_f,
  perma_tune_length  # can only be one number.
) {
  
  library(EnhancedVolcano)

  ############## Prepare sub_mat.
  # Filter by at least p_s sample have > n_c reads.
  ## For each gene, > p_s cases have raw counts > n_c.
  n_c = 100   # n_count
  p_s = 0.2
  
  n_s = floor(dim(sub_mat)[1] * p_s) # in at least 10% sample..
  filter = apply(sub_mat, 2, function(x) {length(x[x>n_c]) > n_s})     # 2 is column, column is gene.
  table(filter, useNA = "always")
  sub_mat = sub_mat[, filter]
  dim(sub_mat)  # 30329 603 # 31366*603 (8.5.2018)    # 20375 * 603 (2.12.2019)  # 73 17719 (Her2 3.7.2019)
  # sub_mat[1:3, 1:3]
  
  
# ############## DESeq2 give genes.
# cat('\n')
# cat("Runing DESeq2 ...", '\n')
# genes_deseq = deseq_give_genes(
#   mat_input=sub_mat,
#   recu_input=sub_recu,
#   resample_n=resample_n,
#   deseq_p=deseq_p,
#   deseq_f=deseq_f
# )
# cat("DESeq2 done, it gives ", length(genes_deseq), " genes.", '\n')
#
# ############## Cox
# cat("\n")
# print("Running Cox ...")
# cox_genes = cox_give_genes(
#   mat_input = sub_mat[, genes_deseq],      # input raw count
#   recu_input = sub_recu,
#   resample_n = resample_n
# )
# cat("Cox run done, it gives ", length(cox_genes), " genes.", "\n")
# saveRDS(cox_genes, paste0("genes_cox_", resample_n, ".RDS"))
  
   cat("Read cox genes from rds", '\n')
   cox_genes = readRDS(paste0("genes_cox_", resample_n, ".RDS"))
   cat("Cox genes length ", length(cox_genes), "\n")

  ############# perma

  # Only when cox_gene is less than 20, continue.
  if (length(cox_genes)<=30) {
   cat("Running perma using forward selection and taking AUC and D.", '\n')
   perma_genes = perma_give_genes(
    mat_input = sub_mat[, cox_genes],        # input raw count.
    recu_input = sub_recu,
    perma_tune_length=perma_tune_length
  )
   cat("Perma run done, it gives by auc ", length(perma_genes$auc), " genes.", "\n")
   # cat("Perma run done, it gives by D ", length(perma_genes$D), " genes.", "\n")
#    saveRDS(perma_genes, paste0("genes_perma_resample_", resample_n, "_tuneLen_", perma_tune_length, ".RDS"))
#    
#    cat("Read perma genes from RDS", '\n')
#    perma_genes = readRDS(paste0("genes_perma_resample_", resample_n, "_tuneLen_", perma_tune_length, ".RDS"))
#    cat("Perma genes length ", length(perma_genes$auc), "\n")

   # print(perma_genes)


  } else {stop("cox_gene > 30, stop.")}
  
 
  return(perma_genes)
  
}
