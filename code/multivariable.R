# set working directory with respect to OS.
if (grepl('darwin', version$os)) {
  setwd("/Users/hl/Documents/breast_cancer/write_paper_reDo_3/4_multivariable/")
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
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

source("../2_KM_on_gene_sig/source_read_in_get_predicted_class.R")


# Parameters
# data_type = 'tcga'   # tcga g5
data_type = "g5"
gene_type = "myGenes"
sub_type = "subtype_basal"  # subtype_basal subtype_lum_a
node_status = 'Negative'   # Negative Positive
my_meth = 'together'
my_meth_by = 'AUC'


# Get predicted class.
read_in_name = paste(
  "../2_KM_on_gene_sig/pred_class", sub_type, my_meth, "by", my_meth_by, 
  "dataType", data_type, "geneType", gene_type, "node_status", node_status, ".RDS", sep = "_"
)
pred_class = get_predicted_class(read_in_name)
rownames(pred_class) = pred_class$sample

# Use pred_prob_Yes
# Perma score is the number between 0 and 1 and 
# transfor pred_prob_Yes by 
# 1. log transformatin then 
# 2. min-max normalization.
# 
pred_class$log_min_max = (log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes))) / (
  max(log(pred_class$pred_prob_Yes) - min(log(pred_class$pred_prob_Yes)))
)

###################################################
# 
# Get clinical info
#
############      TCGA        #####################
# data_set = readRDS("../0.1.1_data__tcga/xml_clinical_patient.RDS")
# clin_pred = data_set[, c("bcr_patient_barcode", "days_to_birth", "race_list", 
#                          "history_of_neoadjuvant_treatment", "ethnicity", "menopause_status", 
#                          "stage_event_pathologic_stage", "stage_event_tnm_categories",
#                          "age_at_initial_pathologic_diagnosis")]
# clin_pred = clin_pred[!duplicated(clin_pred$bcr_patient_barcode), ]
# rownames(clin_pred) = clin_pred$bcr_patient_barcode
# clin_pred = clin_pred[rownames(pred_class), ]
# 
# # Make age numeric.
# clin_pred$age_numeric = as.numeric(clin_pred$age_at_initial_pathologic_diagnosis)
# # Get stage.   make it numeric.
# clin_pred$stage = as.numeric(clin_pred$stage_event_pathologic_stage)
# clin_pred$stage[clin_pred$stage == 1] = NA   #

############     G5           ####################
data_name = paste0("../0.2.1_commen_genes_tcga_gh_others/mat_recu_g5.RDS")
data_set = readRDS(data_name)
clin_pred = data_set[, c("id", 'age.at.diagnosis', 'tumor.size', 'nhg', 'endocrine.treated', 'chemo.treated', 'nodePos')]
clin_pred = clin_pred[rownames(pred_class), ]

# Make age numeric
clin_pred$age_numeric = sapply(as.character(clin_pred$age.at.diagnosis), function(x){strsplit(x, ":")[[1]][2]})
clin_pred$age_numeric = as.numeric(clin_pred$age_numeric)
# Make tumor size numeric
clin_pred$tumor_size_numeric = sapply(as.character(clin_pred$tumor.size), function(x){strsplit(x, ":")[[1]][2]})
clin_pred$tumor_size_numeric = as.numeric(clin_pred$tumor_size_numeric)
# Get tumor grade.
clin_pred$grade = sapply(as.character(clin_pred$nhg), function(x){strsplit(x, ": G")[[1]][2]})
clin_pred$grade = as.numeric(clin_pred$grade)
# clin_pred$grade[clin_pred$grade == "NA"] = NA
#
clin_pred$endo = sapply(as.character(clin_pred$endocrine.treated), function(x){strsplit(x, ": ")[[1]][2]})
clin_pred$endo = as.numeric(clin_pred$endo)
#
clin_pred$chemo = sapply(as.character(clin_pred$chemo.treated), function(x){strsplit(x, ": ")[[1]][2]})
clin_pred$chemo = as.numeric(clin_pred$chemo)

####################################################

# combine
pred_class = cbind(pred_class, clin_pred)

# prepare surbOject
pred_class$survObject = with(pred_class, Surv(time = years, event = (true_class_Yes == "Yes")))

###################################################
# 
# Multivariate analysis.
#
###################################################

# For tcga
# covariate_list = vars(age_numeric, stage, log_min_max)

# For g5
covariate_list = vars(age_numeric, tumor_size_numeric, grade, endo, log_min_max)

##################################################
# Factor labeller
factor_labels = c(age_numeric="Age at diagnosis", log_min_max="Perma Score", endo='Endocrine treatment', 
                  chemo='Chemotherapy', tumor_size_numeric='Tumor size', grade='Histologic grade',
                  stage="Pathologic stage")

pred_class %>%
  analyse_multivariate(
    vars(years, true_class),
    covariates = covariate_list) -> 
  res_var

forest_plot(
  res_var,
  factor_labeller = factor_labels,
  endpoint_labeller = c(years='RFS'),
  orderer = ~order(HR, decreasing =TRUE),
  labels_displayed = c("factor", 'n'),
  label_headers = c(endpoint='Endpoint', factor="Factor", n='n'),
  values_displayed = c("HR", "CI", "p"),
  relative_widths = c(1,1.7,1),
  ggtheme = theme_bw(base_size = 17),
  HR_x_breaks = c(0.1, 1, 10, 100, 1000, 10000),
  HR_x_limits = c(0.1, 15000),
  psprintfFormat = "%.5f",
  p_lessthan_cutoff = 0.00001) -> fp1

# fp1 = fp1 + scale_x_continuous(labels = c(0.1, 1))
# ggsave("forest_plot_multivar.pdf", fp1, scale = 0.7, width = 20, height = 10, dpi = 300)

########################################################
#
# Multiple univariate analysis.
#
########################################################
map(covariate_list, function(by) {
  analyse_multivariate(
    pred_class, 
    vars(years, true_class),
    covariates = list(by))
}) %>%
  forest_plot(
    factor_labeller = factor_labels,
    endpoint_labeller = c(years='RFS'),
    orderer = ~order(HR, decreasing =TRUE),
    labels_displayed = c("factor", 'n'),
    label_headers = c(endpoint='Endpoint', factor="Factor", n='n'),
    values_displayed = c("HR", "CI", "p"),
    relative_widths = c(1,1.7,1),
    ggtheme = theme_bw(base_size = 17),
    HR_x_breaks = c(0.1, 1, 10, 100, 1000, 10000),
    HR_x_limits = c(0.1, 15000),
    psprintfFormat = "%.5f",
    p_lessthan_cutoff = 0.00001) -> fp2

# ggsave('forest_plot_multi_uni.pdf', fp2, scale=0.7, width = 20, height = 10, dpi=300)

# Combine multi var and uni var.
fp3 = ggplot() + theme_void()
fp4 = fp3
fp_all = ggarrange(fp4, fp2, fp3, fp1, ncol = 1, nrow = 4, labels = c("", "Univariate Analysis", "", "Multivariate Analysis"), 
                   hjust=-0.1, vjust=-0.3, heights = c(0.3,2,0.3,2))
save_name = paste("forest_plot_together_", data_type, sub_type, "node", node_status, '.pdf', sep='_')
ggsave(save_name, scale = 0.7, width = 20, height = 10, dpi=300)
saveRDS(fp_all, paste0(save_name, ".RDS"))

