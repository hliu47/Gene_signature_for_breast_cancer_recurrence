get_surv_analysis = function(predicted_class) {
  conf_mat = confusionMatrix(
    predicted_class$pred_class,
    predicted_class$true_class_Yes,
    positive = "Yes"
  )
  
  print(conf_mat)
  
  res_spec = conf_mat$byClass['Specificity']
  res_sens = conf_mat$byClass['Sensitivity']
  res_kappa = conf_mat$overall['Kappa']
  res_acc = conf_mat$overall['Accuracy']
  res_posPrev = conf_mat$byClass['Prevalence']
  
  
  # cox fit to get hr
  cox_fit = coxph(survObject ~ pred_class, data = predicted_class)
  # print(2)
  # cox_fit
  sum_cox = summary(cox_fit)
  cox_hr = sum_cox$conf.int[1]
  cox_hr_95l = sum_cox$conf.int[3]
  cox_hr_95u = sum_cox$conf.int[4]
  
  # surv fit to prepare plot.
  surv_fit = survfit(survObject ~ pred_class, data = predicted_class, conf.type = 'log-log')
  # print(3)
  # get some surv data
  km_data = ggsurvplot(surv_fit, data = predicted_class)
  # print(km_data$data.survplot)
  # 10 and 5 year survial rate for Yes and No.
  surv_rate = NULL
  for (class in c("Yes", "No")) {
    # print(class)
    for (year in c(5.1, 10)) {
      # print(year)
      temp = km_data$data.survplot[km_data$data.survplot$pred_class == class & km_data$data.survplot$time <= year, ]
      temp = temp[dim(temp)[1], ]
      temp = temp[c("surv", 'lower', 'upper')]
      rownames(temp) = paste(year, 'year', 'for', class, sep = "_")
      surv_rate = rbind(surv_rate, temp)
      # print(surv_rate)
    }
  }
  # print(4)
  return(list(acc=res_acc, kappa=res_kappa, sens=res_sens, spec=res_spec, 
              cox_hr=cox_hr, cox_hr_95l=cox_hr_95l, cox_hr_95u=cox_hr_95u,
              surv_rate = surv_rate, 
              surv_fit = surv_fit))
}

