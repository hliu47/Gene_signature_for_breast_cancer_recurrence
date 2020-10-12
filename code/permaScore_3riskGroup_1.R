# set working directory with respect to OS.
if (grepl('darwin', version$os)) {
  setwd("/Users/hl/Documents/breast_cancer/write_paper_reDo_3/3_perma_score_by_AUC/")
} else {
  # setwd("/home/ubuntu/haipeng/write_paper_reDo_2/4_perma_score/")
}

library(caret)
library(ggplot2)
library(mlbench)
# library(DESeq2)
library(pROC)
library(e1071)
library(ROCR)
library(survival)
library(glmnet)
library(txtplot)
library(survminer)
library(mclust)

source("../2_KM_on_gene_sig/source_read_in_get_predicted_class.R")


# Parameters
data_type = 'tcga'
gene_type = "myGenes"
sub_type = "subtype_basal"
node_status = 'Negative'
my_meth = 'together'
my_meth_by = 'AUC'

# Get predicted class.
read_in_name = paste(
  "../2_KM_on_gene_sig/pred_class", sub_type, my_meth, "by", my_meth_by, 
  "dataType", data_type, "geneType", gene_type, 'node_status', node_status, ".RDS", sep = "_"
)
pred_class = get_predicted_class(read_in_name)

# pred_class$true_class_Yes = "No"
# pred_class$true_class_Yes[pred_class$true_class == 1] = "Yes"

#
table(pred_class$pred_class, useNA = 'always')
table(pred_class$true_class_Yes, useNA = 'always')

confusionMatrix(
  factor(pred_class$pred_class, levels = c("No", "Yes")),
  factor(pred_class$true_class_Yes, levels = c("No", "Yes")),
  positive = "Yes"
)

# # # remove some sample for Normal, this sample has true recu but pred_prob is very very low.
# if (sub_type == 'subtype_basal') {
#   pred_class = pred_class[pred_class$sample != "TCGA-BH-A1FC", ]  # remove this sample     old:TCGA-A1-A0SK
# }
# # if (sub_type == 'Normal') {
# #   pred_class = pred_class[pred_class$sample != "TCGA-EW-A1P7", ]  # remove this sample
# # }



pred_class_order = pred_class[order(pred_class$pred_prob_Yes), ]

# Use pred_prob_Yes
# Perma score is the number between 0 and 1 and 
# transfor pred_prob_Yes by 
# 1. log transformatin then 
# 2. min-max normalization.
# 
pred_class$log_prob = log(pred_class$pred_prob_Yes)

pred_class$log_minmax = (log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes))) / (
  max(log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes)))
)



#  transform pred_prob by x*sqrt(exp(1/x))
# pred_class$log_minmax = pred_class$pred_prob_Yes * sqrt(exp(1/(pred_class$pred_prob_Yes)))   # x*sqrt(exp(1/x))
# pred_class$log_minmax = (pred_class$pred_prob_Yes)^(log(pred_class$pred_prob_Yes))   # x^log(x)

par(mfrow = c(2,2))
hist(pred_class$pred_prob_Yes, breaks = 50)
hist(pred_class$log_prob, breaks = 50)

hist(pred_class$log_minmax, breaks = 50)
# hist(pred_class$log_minmax, breaks =50)


dens_mclust = densityMclust(
  G = 3,
  pred_class$log_minmax
  # pred_class$log_prob
)
summary(dens_mclust)
# plot(dens_mclust)

plot(dens_mclust, data = pred_class$log_minmax, what = 'density', breaks = 50, xlab = 'Perma score')
# plot(dens_mclust, data = pred_class$log_prob, what = 'density', breaks = 50, xlab = 'Perma score')



# find the cutoff
cut_off = data.frame(data = dens_mclust$data, classification = dens_mclust$classification)
cut_off = cut_off[order(cut_off$data), ]

# if (sub_type == 'Her2') {
#   cut_off[1,2] = 1
#   cut_off[2,2] = 1
# } else if (sub_type == 'Normal') {
#   cut_off[1,2] = 1
#   # cut_off[2,2] = 1
# }
cutoff1_set = max(cut_off$data[cut_off$classification == 1])
cutoff1_set
cutoff2_set = max(cut_off$data[cut_off$classification == 2])
cutoff2_set
plus_minus = (max(cut_off$data) - min(cut_off$data))/50
# 1        2          3
#    0.06      0.24
#    0.07      0.29 
#    0.09      0.33

#    0.07      0.32 is best.  for minmax of log.

# cutoff1 = 2.98
# cutoff2 = 4.45



save_path = paste0("./",sub_type, "_node_", node_status, '_by_', my_meth_by, "_three_parts_by_log_minmax/")
if (!dir.exists(save_path)) dir.create(save_path)

print(pred_class)

hr_high_inter = data.frame(cutoff1=1, cutoff2=1, hr_high_inter=0)

# for (cutoff1 in seq(cutoff1_set - 3*plus_minus, cutoff1_set + 3*plus_minus, plus_minus)) {
for (cutoff1 in seq(0.5057, 0.5057, 0.01)) {
    
  # for (cutoff2 in seq(cutoff2_set - 3*plus_minus, cutoff2_set + 3*plus_minus, plus_minus)) {
  for (cutoff2 in seq(0.706, 0.706, 0.001)) {
      
    # Best tcga myGenes basal node-.
    # cutoff1 = 0.5057
    # cutoff2 = 0.706
    
    dens_plot = 
      ggplot(pred_class, aes(x=log_minmax)) +
      geom_histogram(aes(y=..density..), bins = 20, fill='grey') +
      geom_density(alpha=0.1, fill='#FF6666') +
      scale_x_continuous(name='Perma Risk Score', breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(name='Density') +
      geom_vline(xintercept = c(cutoff1, cutoff2), linetype=2, colour=c('blue', 'red'), size=1) +
      theme_bw(base_size = 22)

    dens_save_name = paste(save_path, 'dens_mclus_3_parts', 
                           data_type, gene_type, sub_type, node_status,
                           format(cutoff1, digits = 4),
                           format(cutoff2, digits = 4), 
                           sep='_', '.pdf')
    ggsave(dens_save_name, dens_plot, scale = 1, width = 15, height = 15, dpi = 300)
    saveRDS(dens_plot, paste0(dens_save_name, '.RDS'))
    
    # assign 3 risk group, 
    pred_class$risk_groups = "Low"
    pred_class$risk_groups[pred_class$log_minmax > cutoff1 & pred_class$log_minmax <= cutoff2] = "Intermediate"
    pred_class$risk_groups[pred_class$log_minmax > cutoff2] = 'High'
    
    pred_class$risk_groups = factor(pred_class$risk_groups, levels = c("Low", "Intermediate", "High"))
    
    
    
    # KM plot for 3 parts.
    
    # prepare surbOject
    pred_class$survObject = with(pred_class, Surv(time = years, event = (true_class_Yes == "Yes")))
    
    
    # cox fit to get hr
    cox_fit = coxph(survObject ~ risk_groups, data = pred_class )
    sum_cox = summary(cox_fit)
    
    # taking Low as reference group.
    cox_hr_inter = sum_cox$conf.int[1, 1]
    cox_hr_95l_inter = sum_cox$conf.int[1, 3]
    cox_hr_95u_inter = sum_cox$conf.int[1, 4]
    
    cox_hr_high = sum_cox$conf.int[2, 1]
    cox_hr_95l_high = sum_cox$conf.int[2, 3]
    cox_hr_95u_high = sum_cox$conf.int[2, 4]
    
    # taking intermediate as reference
    pred_class$risk_groups_inter = factor(pred_class$risk_groups, levels = c("Intermediate", "High", "Low"))
    cox_fit = coxph(survObject ~ risk_groups_inter, data = pred_class )
    sum_cox = summary(cox_fit)
    
    cox_hr_high_interRef = sum_cox$conf.int[1, 1]
    cox_hr_95l_high_interRef = sum_cox$conf.int[1, 3]
    cox_hr_95u_high_interRef = sum_cox$conf.int[1, 4]
    
    if (cox_hr_95l_inter > 1 & cox_hr_95l_high_interRef > 1 ) {
      cat(cox_hr_95l_inter, cox_hr_95l_high_interRef, "\n")
      cat(cutoff1, cutoff2, "\n")
      cat("\n")
    }
    
    if (cox_hr_95l_high_interRef > 0.8 ) {
      cat(cox_hr_95l_inter, cox_hr_95l_high_interRef, "\n")
      cat(cutoff1, cutoff2, "\n")
      cat("\n")
    }
    hr_high_inter = rbind(hr_high_inter, c(cutoff1, cutoff2, cox_hr_95l_high_interRef))
    
    # surv fit to prepare plot.
    surv_fit = survfit(survObject ~ risk_groups, data = pred_class, conf.type = 'log-log')
    
    
    # get some surv data
    km_data = ggsurvplot(surv_fit, data = pred_class)
    
    # 10 and 5 year survial rate for Yes and No.
    surv_rate = NULL
    for (risk_group in c("Low", "Intermediate", "High")) {
      for (year in c(5, 10)) {
        
        # risk_group = "Low"
        # year = 10
        
        temp = km_data$data.survplot[km_data$data.survplot$risk_groups == risk_group & km_data$data.survplot$time <= year, ]
        temp = temp[dim(temp)[1], ]
        temp = temp[c("surv", 'lower', 'upper')]
        rownames(temp) = paste(year, 'year', 'for', risk_group, sep = "_")
        surv_rate = rbind(surv_rate, temp)
      }
    }
    
    
    
    km_plot = ggsurvplot(
      # fit,                     # survfit object with calculated statistics.
      surv_fit,
      # data = BRCAOV.survInfo,  # data used to fit survival curves.
      data = pred_class,
      
      # fun = "cumhaz", # log or cumhaz or event
      
      # linetype = "strata", # Change line type by groups
      # surv.geom = geom_line,
      surv.geom = geom_step,
      # ncensor.plot = TRUE, #to plot the number of censored subjects at time t.
      
      risk.table = TRUE,          # show risk table.
      risk.table.col = "strata",  # Risk table color by groups
      
      pval = TRUE,             # show p-value of log-rank test.
      # conf.int = TRUE,         # show confidence intervals for point estimaes of survival curves.
      xlim = c(0,10),        # present narrower X axis, but not affect survival estimates.
      break.time.by = 1,       # break X axis in time intervals by 500.
      
      # palette = c("#E7B800", "#2E9FDF"), # color
      palette = c('blue', '#A04000', 'red'),
      # palette = c('#239B56', '#A04000'),
      
      xlab = 
        paste('Years  ',
              # "\n",
              sub_type,
              # "\n",
              # 'AUC:', format(res_auc, digits = 4),   # no AUC, because it give direct class, no probability.
              # "Spec:", format(res_spec, digits = 4),
              # "Sens:", format(res_sens, digits = 4),
              # 'Acc:', format(res_acc, digits = 4),
              # 'PosPrev', format(res_posPrev, digit = 4),
              # "kappa:", format(res_kappa, digits = 4),
              
              "\n",
              "5 years surv Low", format(surv_rate['5_year_for_Low', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['5_year_for_Low', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['5_year_for_Low', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              "Inter", format(surv_rate['5_year_for_Intermediate', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['5_year_for_Intermediate', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['5_year_for_Intermediate', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              "High", format(surv_rate['5_year_for_High', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['5_year_for_High', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['5_year_for_High', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              
              "\n",
              "10 years surv Low", format(surv_rate['10_year_for_Low', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['10_year_for_Low', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['10_year_for_Low', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              "Inter", format(surv_rate['10_year_for_Intermediate', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['10_year_for_Intermediate', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['10_year_for_Intermediate', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              "High", format(surv_rate['10_year_for_High', 'surv'], digits = 2),
              '95CI  [', format(round(surv_rate['10_year_for_High', 'lower'], 2), nsmall = 2), ',', format(round(surv_rate['10_year_for_High', 'upper'], 2), nsmall = 2), ']', # 95 CI lower
              
              # "\n",
              # # 'Dev:', format(deviance, digits = 4),
              "\n",
              'Hazard Ratio Low ref: Inter', format(round(cox_hr_inter, 2), nsmall = 2), # Hazard ratio.
              '95CI  [', format(round(cox_hr_95l_inter, 2), nsmall = 2), ',', format(round(cox_hr_95u_inter, 2), nsmall = 2), ']', # 95 CI lower
              'High', format(round(cox_hr_high, 2), nsmall = 2), # Hazard ratio.
              '95CI  [', format(round(cox_hr_95l_high, 2), nsmall = 2), ',', format(round(cox_hr_95u_high, 2), nsmall = 2), ']', # 95 CI lower
              "\n",
              "HR Inter Ref: High", format(round(cox_hr_high_interRef, 2), nsamll = 2),
              '95CI  [', format(round(cox_hr_95l_high_interRef, 2), nsmall = 2), ',', format(round(cox_hr_95u_high_interRef, 2), nsmall = 2), ']' # 95 CI lower
              
              # 
        ),
      ylab = "Recurrence free survival",
      # legend.labs = c('No-Recurrence, predicted', 'Recurrence, predicted'),
      legend = 'right',
      # ggtheme = theme_minimal(), # customize plot and risk table with a theme.
      ggtheme = theme_light(),
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE # show bars instead of names in text annotations in legend of risk table
    )
    
    km_plot$plot = km_plot$plot +
      geom_segment(aes(x = 5, y = 0, xend = 5, yend = 1), linetype = "dashed", size = 0.4) +
      geom_segment(aes(x = 10, y = 0, xend = 10, yend = 1), linetype = "dashed", size = 0.4) 
      # geom_line(aes(color = group)) +
      # geom_smooth(aes(color = group), se=F)
    
    # add cox_hr_inter + cox_hr_high
    cox_hr_inter_high = cox_hr_inter + cox_hr_high
    
    ggsave(
      paste(
        save_path, 
        "km_plot_3_risk_groups", 
        sub_type, 
        format(cutoff1, digits = 4), 
        format(cutoff2, digits = 4), 
        "hr_Inter_High", 
        format(cox_hr_inter_high, digits=5), 
        ".pdf", 
        sep = "_"
        ), 
      print(km_plot, newpage = F), 
      width = 10, 
      height = 15
      )
    saveRDS(
      km_plot, 
      paste(
        save_path, 
        "km_plot_3_risk_groups", 
        sub_type, 
        format(cutoff1, digits = 4), 
        format(cutoff2, digits = 4), 
        "hr_Inter_High", 
        format(cox_hr_inter_high, digits=5), 
        ".RDS", 
        sep = "_"
      ), 
    )
    
  }
}


print(hr_high_inter[order(hr_high_inter$hr_high_inter), ])
