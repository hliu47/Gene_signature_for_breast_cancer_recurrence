# set working directory with respect to OS.
if (grepl('darwin', version$os)) {
  setwd("/Users/hl/Documents/breast_cancer/write_paper_reDo_3/3_perma_score_by_AUC/")
} else {
  # setwd("/home/ubuntu/haipeng/write_paper_reDo_2/4_perma_score/")
}

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
library(survminer)
library(mclust)


source("../2_KM_on_gene_sig/source_read_in_get_predicted_class.R")


# Parameters
data_type = 'tcga'
gene_type = "myGenes"
sub_type = "subtype_basal"
node_status = "Negative"
my_meth = 'together'
my_meth_by = 'AUC'

# Get predicted class.
read_in_name = paste(
  "../2_KM_on_gene_sig/pred_class", sub_type, my_meth, "by", my_meth_by, 
  "dataType", data_type, "geneType", gene_type, "node_status", node_status, ".RDS", sep = "_"
)
pred_class = get_predicted_class(read_in_name)


# prepare data
# sub_type = "subtype_basal"     # need to change
# best 
cutoff1 = 0.5057  # need to change 
cutoff2 = 0.706  # need to change

#         LumA  LumB  Her2  Basal Normal
# cutoff1 0.279 0.18  0.25  0.176 0.4298
# cutoff2 0.781 0.65  0.722 0.612 0.831

# sub_by = 'AUC'
# pred_class = readRDS(paste("../3_KM_on_gene_sig/pred_class", sub_type, "together_by", sub_by, ".RDS", sep = "_"))
# pred_class$true_class_Yes = "No"
# pred_class$true_class_Yes[pred_class$true_class == 1] = "Yes"


# Use pred_prob_Yes
# Perma score is the number between 0 and 1 and 
# transfor pred_prob_Yes by 
# 1. log transformatin then 
# 2. min-max normalization.
# 
pred_class$log_min_max = (log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes))) / (
  max(log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes)))
)

# hist(pred_class$log_min_max, breaks = 50)
# 
# dens_mclust = densityMclust(pred_class$log_min_max, G =3)
# summary(dens_mclust)
# # plot(dens_mclust)
# 
# # find the cutoff
# cut_off = data.frame(data = dens_mclust$data, classification = dens_mclust$classification)
# cut_off = cut_off[order(cut_off$data), ]

# 1        2          3
#    0.06      0.24
#    0.07      0.29 
#    0.09      0.33

#    0.07      0.32 is best.



# pdf(paste('dens_mclus_3_parts', sub_type, cutoff1, cutoff2, sep='_', '.pdf'))
# plot(dens_mclust, data = pred_class$log_min_max, what = 'density', breaks = 50, xlab = 'Perma score')
# abline(v = cutoff1, col = 'blue', lty = 2)
# abline(v = cutoff2, col = 'red', lty = 2)
# title(paste("Density and 3 parts cutoff", cutoff1, cutoff2))
# dev.off()

# assign 3 risk group, 
pred_class$risk_groups = "Low"
pred_class$risk_groups[pred_class$log_min_max > cutoff1 & pred_class$log_min_max <= cutoff2] = "Intermediate"
pred_class$risk_groups[pred_class$log_min_max > cutoff2] = 'High'

pred_class$risk_groups = factor(pred_class$risk_groups, levels = c("Low", "Intermediate", "High"))



### predict distribution of surv for perma_score from 0.01 to 0.99

# prepare surbOject
pred_class$survObject = with(pred_class, Surv(time = years, event = (true_class_Yes == "Yes")))
cox_fit = coxph(survObject ~ log_min_max, data = pred_class )   # log_min_max is perma_score.

# surv_fit = survfit(cox_fit, newdata = data.frame(log_min_max = 0.01))
# a = ggsurvplot(surv_fit, data = pred_class)
# a = plot(surv_fit)

# for 5 and 10 year Recurrence free probability

perma_score_surv = data.frame(
  perma_score = seq(0.001, 0.999, 0.001),
  surv_5_year = NA
)

surv_all = sapply(perma_score_surv$perma_score, function(x) {
  
  # x = 0.5
  surv_fit = survfit(cox_fit, newdata = data.frame(log_min_max = x))
  
  # surv_plot = ggsurvplot(surv_fit, data = pred_class)
  # surv_data = surv_plot$data.survplot
  
  surv_data = data.frame(time=surv_fit$time, surv=surv_fit$surv, lower=surv_fit$lower, upper=surv_fit$upper)
  
  temp1 = surv_data[surv_data$time <=5, ]     # at 5 year 
  temp1 = temp1[dim(temp1)[1], ]
  
  temp2 = surv_data[surv_data$time <=10, ]     # at 10 year 
  temp2 = temp2[dim(temp2)[1], ]
  
  data_ret = c(temp1$surv, temp1$lower, temp1$upper, temp2$surv, temp2$lower, temp2$upper)
  
  return(data_ret)
})

rownames(surv_all) = c('surv_5', 'lower_5', 'upper_5', 'surv_10', 'lower_10', 'upper_10')
colnames(surv_all) = perma_score_surv$perma_score

surv_all = t(surv_all)
surv_all = as.data.frame(surv_all)
surv_all$perma_score = rownames(surv_all)
surv_all[1:5, ]
dim(surv_all)    # 999 7



###############   ggplot

library(reshape2)

temp = surv_all

#########
dim(temp)
temp[1:3, ]
colnames(temp)[1:6] = rep(c('surv', 'lower', 'upper'), 2)
a = temp[, c(1,2,3,7)]
b = temp[, c(4,5,6,7)]
temp = rbind(a, b)
temp$year = 10
temp$year[1:999] = 5

temp = melt(temp, id = c('perma_score', 'year'))



########################################## with conf int
gg_plot = ggplot(temp, aes(x = perma_score, y = value, group = variable, color = year)) +
  geom_line(data = temp[temp$year == 5 & temp$variable == 'surv', ], linetype = 1, size = 1)  +
  geom_line(data = temp[temp$year == 5 & temp$variable == 'lower', ], linetype = 3, size = 0.7)  +
  geom_line(data = temp[temp$year == 5 & temp$variable == 'upper', ], linetype = 3, size = 0.7)  +

  geom_line(data = temp[temp$year == 10 & temp$variable == 'surv', ], linetype = 1, size = 1)  +
  geom_line(data = temp[temp$year == 10 & temp$variable == 'lower', ], linetype = 3, size = 0.7)  +
  geom_line(data = temp[temp$year == 10 & temp$variable == 'upper', ], linetype = 3, size = 0.7)  +

  scale_color_continuous(name = 'Year', breaks = c(5, 10)) +
  guides(color = guide_legend()) +

  # scale_x_discrete(breaks = c(0.07, seq(0.1, 1, 0.1), 0.32)) +
  scale_x_discrete(breaks = seq(0.1, 1, 0.1)) + 
  
  # theme(axis.text.x = element_text(angle=45)) +

  scale_y_continuous(breaks = pretty(seq(0, 1, 0.1), 10)) +
  geom_vline(xintercept = cutoff1 * 1000, linetype = 2, color = 'blue', size = 0.7) +   # 0.07
  geom_vline(xintercept = cutoff2 * 1000, linetype = 2, color = 'blue', size = 0.7) +  # 0.32

  labs(
    title = paste0("5 and 10 year recurrence free survival VS Perma score"),
    x = 'Perma Risk Score',
    y = paste0('Recurrence free survival')) +

  theme_bw()  +
  theme(legend.position = 'right')

ggsave(
  paste("Perma_score_5_10_year_surv_with_confInt", sub_type, ".pdf", sep = "_"),
  # print(gg_plot, newpage = FALSE),
  gg_plot,
  scale = 0.8,
  width = 15,
  height = 10)


####################################### without conf int
gg_plot = ggplot(temp, aes(x = perma_score, y = value, group = variable, color = year)) +
  geom_line(data = temp[temp$year == 5 & temp$variable == 'surv', ], linetype = 1, size = 1)  +
  geom_line(data = temp[temp$year == 10 & temp$variable == 'surv', ], linetype = 1, size = 1)  +
  geom_area(data = temp[temp$year == 5 & temp$variable == 'surv', ], aes(fill='year'), alpha=0.1) +
  scale_color_continuous(name = 'Year', breaks = c(5, 10)) +
  guides(color = guide_legend()) +
  
  # scale_x_discrete(breaks = c(cutoff1, seq(0.1, 1, 0.1), cutoff2)) + 
  scale_x_discrete(breaks = seq(0.1, 1, 0.1)) + 
  # theme(axis.text.x = element_text(angle=45)) +
  
  scale_y_continuous(breaks = pretty(seq(0, 1, 0.1), 10)) + 
  geom_vline(xintercept = cutoff1 * 1000, linetype = 2, color = 'blue', size = 0.7) +   # 0.07
  geom_vline(xintercept = cutoff2 * 1000, linetype = 2, color = 'blue', size = 0.7) +  # 0.32
  
  labs(
    title = paste0("5 and 10 year recurrence free survival VS Perma score"),
    x = 'Perma risk score', 
    y = paste0('Recurrence free survival')) +
  
  theme_bw()  +
  theme(legend.position = 'right') 

ggsave(
  paste("Perma_score_5_10_year_surv", sub_type, ".pdf", sep = "_"), 
  # print(gg_plot, newpage = FALSE), 
  gg_plot,
  scale = 0.8,
  width = 15, 
  height = 10)






