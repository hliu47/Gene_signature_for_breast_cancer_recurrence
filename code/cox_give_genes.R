cox_give_genes = function(
  mat_input,
  recu_input,
  resample_n
) {
  
  sub_mat = mat_input
  sub_recu = recu_input
  
  ## VST
  mat_vst =t(varianceStabilizingTransformation(t(sub_mat), fitType="mean", blind = FALSE))     # vst transformation.

  dim(mat_vst)

  ##############################################
  family = "cox" #######
  x = mat_vst
  y = survival::Surv(time = sub_recu$days, event =  sub_recu$recu1)

  # for foldid
  n_obs = dim(x)[1]
  # 3 fold cv. if 10 fold only 111 samples may be too small. if many samples 10 fold is ok.
  fold_n = 5            
  # cat("CV fold is ", fold_n, "\n")
  n_q = n_obs %/% fold_n
  n_r = n_obs %% fold_n
  if(n_r != 0) {
    foldid = as.vector(c(rep(1:fold_n, n_q), seq(1, n_r))) 
  } else if(n_r == 0) {
    foldid = as.vector(c(rep(1:fold_n, n_q))) 
  }

  
  # alpha is for the elastic-net mixing parameter α, with range α∈[0,1]. 
  # α=1 is the lasso (default) and α=0 is the ridge.
  #alpha_1 = 0.001   # alpha=0 is ridge.     
  #alpha_2 = 1  # alpha=1 is lasso.     
  #alpha_by = 0.01      ###########################  Change ###########
  
  set.seed(123)
  alpha_seq = rgamma(100, 0.5)   # 200
  alpha_seq = alpha_seq[alpha_seq>=0]  # 100 > 0. 
  alpha_seq = alpha_seq/max(alpha_seq)
  alpha_seq = round(alpha_seq, 6) 
  alpha_seq = alpha_seq[order(alpha_seq)] 
#  print(alpha_seq)
  # print("alpha=1, Lasso")
  #cat('alpha from', alpha_1, " to ", alpha_2, ' by ', alpha_by, "\n")
  cat("alpha seq len is ", length(alpha_seq), " and min is ", min(alpha_seq), " max is ", max(alpha_seq),  "\n")
  # nlambda = 200
  # cat("nlambda is ", nlambda, "\n")
  # 
  # lambda_1 = -3
  # lambda_2 = -2
  # lambda_by = 0.01
  # cat("lambda from ", lambda_1, ' to ', lambda_2, " by ", lambda_by, "\n")
  
  # Select alpha.
  cvm_alpha_all = foreach (
    #i = seq(alpha_1, alpha_2, by = alpha_by),
    i = alpha_seq,
    # .verbose = T,
    .combine = rbind
  ) %dopar% {
    
    alpha = i
    
    # print(alpha)
    
    # cat(paste("alpha = ", alpha))
    # cat('\n')
    
    set.seed(123)
    # 10 fold cv, measure metric is deviance. 
    cvfit = glmnet::cv.glmnet(
      x = x,
      y = y,   
      family = family, 
      parallel = F, 
      alpha = alpha,
      foldid = foldid
      # nlambda = nlambda  # NOT important
      #  lambda = seq(lambda_1, lambda_2, by = lambda_by)    #### lambda sequence.
    )
  
    # print(cvfit)
    
    lambda_min_cv = cvfit$lambda.min
    lambda_max = max(cvfit$lambda)
    
    coef_min = coef(cvfit, s = "lambda.min")
    gene_keep = which(coef_min != 0)

    # get min deviance
    min_cvm = cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
    
    # get alpha, min_cvm, number of genes selected.
    cvm_temp = c(alpha, min_cvm, length(gene_keep)) 
    
    return(cvm_temp)       # alpha cvm n_genes
    
  }
  
  colnames(cvm_alpha_all) = c('alpha', 'cvm', 'n_genes')
  
#  print(cvm_alpha_all)
  
  cvm_choose_min = cvm_alpha_all[cvm_alpha_all[, 'cvm'] == min(cvm_alpha_all[, 'cvm']), ]
  
  alpha_choose = cvm_choose_min['alpha']
  
  cvm_alpha_all = data.frame(alpha = cvm_alpha_all[, "alpha"], 
                             cvm=cvm_alpha_all[, "cvm"], 
                             n_genes=cvm_alpha_all[, "n_genes"])
  
  min_cvm_alpha= cvm_alpha_all$alpha[with(cvm_alpha_all, cvm==min(cvm))]
  min_cvm_gene=cvm_alpha_all$n_genes[with(cvm_alpha_all, cvm==min(cvm))]
 
  cat("min_cvm_alpha = ", min_cvm_alpha, "\n") 
  
  # Make sure min_cvm_alpha is not 0 nor 1.
  if (min_cvm_alpha == 0 | min_cvm_alpha == 1) stop("min_cvm_alpha is 0 or 1, stop.")
  

  # plot alpha vs cvm
  cox_alpha_plot = 
    ggplot(data.frame(cvm_alpha_all), aes(log(alpha), cvm)) +
    geom_point(colour='black') +
    geom_vline(xintercept = log(min_cvm_alpha), linetype=2) +
    scale_x_continuous(
      sec.axis = sec_axis(
        ~., 
        name='Number of genes', 
        breaks=log(cvm_alpha_all$alpha[seq(1,dim(cvm_alpha_all)[1], 10)]), 
        labels=cvm_alpha_all$n_genes[seq(1, dim(cvm_alpha_all)[1], 10)])) +
    labs(x=expression(log~alpha), y="Partial Likelihood Deviance") +
    theme_bw(base_size = 22)
  
  
  # ggsave(paste0("cox_cvfit_alpha_geneN_", min_cvm_gene, "_resample_", resample_n, ".pdf"), cox_alpha_plot)
  
  # saveRDS(cox_alpha_plot, paste0("cox_cvfit_alpha_", resample_n, ".RDS"))
  
  # run again
  set.seed(123)
  # 10 fold cv, measure metric is deviance. 
  cvfit = glmnet::cv.glmnet(
    x = x,
    y = y,   
    family = family, 
    parallel = F, 
    alpha = alpha_choose,
    foldid = foldid,
    # nlambda = nlambda  # NOT important
    # lambda = exp(seq(lambda_1, lambda_2, by = lambda_by))    #### lambda sequence. exp(raw lambda)
  )
  
  cvfit_df = data.frame(lambda=cvfit$lambda, cvm=cvfit$cvm, n_gene=cvfit$nzero)
  
  cox_plot = 
    ggplot(cvfit_df, aes(log(lambda), cvm)) +
    geom_point(colour='black') +
    geom_vline(xintercept = log(cvfit$lambda.min), linetype=2) +
    scale_x_continuous(
      sec.axis = sec_axis(
        ~., 
        name='Number of genes', 
        breaks=log(cvfit_df$lambda)[seq(1,dim(cvfit_df)[1], 5)], 
        labels=cvfit_df$n_gene[seq(1, dim(cvfit_df)[1], 5)])) +
    labs(x=expression(log~lambda), y="Partial Likelihood Deviance") +
    theme_bw(base_size = 22)
  
    
  # plot(cvfit, se.bands=FALSE)
  # abline(v=log(cvfit$lambda.1se), col='white', lty=1, lwd=2)  # get rid of 1se line.
  # abline(v=log(cvfit$lambda.min), col='black', lty=2, lwd=1)  # emphasize min line.
  # mtext("Number of genes", padj=-5)
  # dev.off()
  
  coef_min = coef(cvfit, s = "lambda.min")
  gene_keep = which(coef_min != 0)
  
  gene_keep_name = rownames(coef_min)[gene_keep]
  
  # ggsave(paste0("cox_cvfit_lambda_geneN_", length(gene_keep_name), "_resample_", resample_n, ".pdf"), cox_plot)
  # saveRDS(cox_plot, paste0("cox_cvfit_lambda_selected_", length(gene_keep_name), "_resample_", resample_n, ".RDS"))
  
  return(gene_keep_name)   # cox give this selected gene names.
  
}



