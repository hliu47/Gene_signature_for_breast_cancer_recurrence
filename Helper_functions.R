# Compute Cohen's kappa statitic. It only takes ROC as input and ouput Kappa and Balanced Acc which is (sens + spec)/2.
compute_kappa = function(res_roc) {   # res_roc should only have one value.
  res_all = coords(
    res_roc,
    "best",
    best.method = "closest.topleft",
    ret = c("sensitivity", "specificity", 'tp', 'tn', 'fp', 'fn', 'threshold')
  )
  # res_all can have two values.
  if (class(res_all) == 'matrix') {
    message('ROC has multiple values, please check ROC.')
    # res_all = res_all[, 1]
  }
  
  # theta below is prevalance    a+c/(a+b+c+d)   # events / all
  theta = (res_all['tp'] + res_all['fn']) / (res_all['tp'] + res_all['fn'] + res_all['tn'] + res_all['fp'])    
  
  sens = res_all['sensitivity']  # sensitivity
  spec = res_all['specificity']  # specificity
  
  res_kappa_up = (2*theta)*(1-theta)*(sens+spec-1)
  res_kappa_down =  (theta)^2 + (1-theta)^2 + (1-2*theta)*(theta*sens - (1-theta)*spec)
  res_kappa = res_kappa_up/res_kappa_down
  res_kappa = as.numeric(res_kappa)
  
  bal_acc = (res_all['sensitivity'] + res_all['specificity'])/2
  bal_acc = as.numeric(bal_acc)
  
  return(c(kappa = res_kappa, bal_acc = bal_acc))
}