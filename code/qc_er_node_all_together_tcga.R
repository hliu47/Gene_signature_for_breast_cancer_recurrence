setwd('/Users/hl/Documents/breast_cancer/write_paper_reDo_4/0.3_add_ER_node_QC/')


library(ggplot2)
library(DESeq2)
library(ggfortify)
library(ggpubr)
library(grid)
library(gridExtra)


# TCGA
# tcga = readRDS('../0.1_data_prepare/tcga_mat_recu_603_20376_10_0.1_.RDS')
# tcga = readRDS("../0.1.1_data__tcga/tcga_hgnc_raw_counts_36862_1081.RDS")
tcga = readRDS("../0.2.1_commen_genes_tcga_gh_others/mat_tcga_hgnc_18877_1081_18Jan2020.RDS")
dim(tcga)  # 18877*1081.
tcga[1:3, 1:2]
#        TCGA.C8.A133.01A.32R.A12D.07 TCGA.C8.A1HJ.01A.11R.A13Q.07
# TSPAN6                         1590                         6321
# TNMD                              7                           53
# DPM1                           1342                         5153

# genes_all = matrix(rownames(tcga), nrow = 439)
# for (i in 1:10000){
#   if (18877 %% i == 0) {
#     print(i)
#   }
# }

colnames(tcga) = sapply(colnames(tcga), function(x) paste(strsplit(x, '\\.')[[1]][1:3], collapse = "-"))

# Recu.
recu = readRDS("../0.1.1_data__tcga/all_status_my_short_1084_7_useThis.RDS")
dim(recu)
recu[1:3, ]
# > recu[1:3, ]
#        patient age recu days Node history_of_neoadjuvant_treatment  nodePos
# 1 TCGA-3C-AAAU  55  Yes 1808    4                               No Positive
# 2 TCGA-3C-AALI  50 <NA> 4005    1                               No Positive
# 3 TCGA-3C-AALJ  62 <NA> 1474    1                               No Positive
sum(duplicated(recu$patient))
sum(duplicated(colnames(tcga)))
common_sample = intersect(recu$patient, colnames(tcga))

recu = recu[recu$patient %in% common_sample, ]
tcga = tcga[, colnames(tcga) %in% common_sample]

recu = recu[order(recu$patient), ]
tcga = tcga[, order(colnames(tcga))]

tcga_recu = cbind(recu, data.frame(t(tcga)))  # sample in row.
dim(tcga_recu) # 1067 18884.
# er ESR1  pr PGR  her2 ERBB2.
gene_of_int = tcga_recu[, c("ESR1", "PGR", "ERBB2")]
gene_of_int$esr1_log = log(gene_of_int$ESR1)
gene_of_int$pgr_log = log(gene_of_int$PGR)
gene_of_int$erbb2_log = log(gene_of_int$ERBB2)

# er
er_gg = 
  ggplot(gene_of_int, aes(esr1_log)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = c(7.40, 8.45), linetype=2, size=1, color='red') +
  labs(title='ER') 
ggsave("tcga_er.pdf", er_gg, scale = 1, width = 10, height = 10, dpi = 300)
saveRDS(er_gg, 'tcga_er.RDS')

# pr
pr_gg =
  ggplot(gene_of_int, aes(pgr_log)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = c(6.1, 7.3), linetype=2, size=1, color='red') +
  labs(title = "PR")
ggsave("tcga_pr.pdf", pr_gg, scale = 1, width = 10, height = 10, dpi = 300)
saveRDS(pr_gg, 'tcga_pr.RDS')

# her2
her2_gg = 
  ggplot(gene_of_int, aes(erbb2_log)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = c(11.3, 11.9), linetype=2, size=1, color='red') +
  labs(title = "HER2")
ggsave("tcga_her2.pdf", her2_gg, scale = 1, width = 10, height = 10, dpi = 300)
saveRDS(her2_gg, 'tcga_her2.RDS')



# er
gene_of_int$er_status = 0
gene_of_int$er_status = sapply(gene_of_int$esr1_log, function(x) ifelse(x>=8.45, 1, ifelse(x>=7.40 & x<8.45, "NotSure", 0)))

# pr
gene_of_int$pr_status = 0
gene_of_int$pr_status = sapply(gene_of_int$pgr_log, function(x) ifelse(x>=7.3, 1, ifelse(x>=6.1 & x<7.3, "NotSure", 0)))

# her2
gene_of_int$her2_status = 0
gene_of_int$her2_status = sapply(gene_of_int$erbb2_log, function(x) ifelse(x>=11.9, 1, ifelse(x>=11.3 & x<11.9, "NotSure", 0)))

#
table(gene_of_int$er_status)
#   0       1 NotSure 
# 233     781      53 
table(gene_of_int$pr_status)
#   0       1 NotSure 
# 368     573     126 
table(gene_of_int$her2_status)
#   0       1 NotSure 
# 961      84      22 



# Combine all.
tcga_recu_status = cbind(gene_of_int, tcga_recu)

# er+ pr+ her2-   101 recu / out of 533.      LumA
subtype_lum_a = tcga_recu_status$er_status == 1 & tcga_recu_status$pr_status == 1 & tcga_recu_status$her2_status == 0
table(subtype_lum_a)
# FALSE  TRUE 
#   534   533 
table(subtype_lum_a & tcga_recu_status$recu == "Yes")
# FALSE  TRUE 
#   534   101 

# er- pr- her2-   36 recu / out of 173        Basal
subtype_basal = tcga_recu_status$er_status == 0 & tcga_recu_status$pr_status == 0 & tcga_recu_status$her2_status == 0
table(subtype_basal)
# FALSE  TRUE 
#   894   173
table(subtype_basal & tcga_recu_status$recu == "Yes")
# FALSE  TRUE 
#   894    36
# 

# er- pr- her2+   9 recu / out of 39          Her2    Not use, too few samples.
subtype_her2 = tcga_recu_status$er_status == 0 & tcga_recu_status$pr_status == 0 & tcga_recu_status$her2_status == 1
table(subtype_her2)
# FALSE  TRUE 
#  1028    39 
table(subtype_her2 & tcga_recu_status$recu == "Yes")
# FALSE  TRUE 
#  1028     9 


# her2
table(tcga_recu_status$her2_status)
#   0       1 NotSure 
# 961      84      22
table(tcga_recu_status$her2_status == 1 & tcga_recu_status$recu == "Yes")     # 18 / 84
table(tcga_recu_status$her2_status == 1 & tcga_recu_status$er_status == 0 & tcga_recu_status$pr_status == 0 & tcga_recu_status$recu == "Yes")  # 10/41
table(tcga_recu_status$her2_status == 1 & tcga_recu_status$nodePos == "Positive" & tcga_recu_status$recu == "Yes")

table(tcga_recu_status$her2_status == 1 & tcga_recu_status$nodePos == "Negative" & tcga_recu_status$recu == "Yes")     # 18 / 84


# er+ pr+ her2+/-     104 recu / out of 550           LumB
subtype_lum_b = tcga_recu_status$er_status == 1 & tcga_recu_status$pr_status == 1 & tcga_recu_status$her2_status != "NotSure"
table(subtype_lum_b)
# FALSE  TRUE 
#   517   550 
table(subtype_lum_b & tcga_recu_status$recu == "Yes")
# FALSE  TRUE 
#   517   104


# Extract
tcga_lum_a = tcga_recu_status[subtype_lum_a, ]

tcga_lum_b = tcga_recu_status[subtype_lum_b, ]
tcga_basal = tcga_recu_status[subtype_basal, ]

# node
table(tcga_lum_a$nodePos)
# Negative Positive 
#      195      264
table(tcga_lum_a$nodePos == "Negative" & tcga_lum_a$recu == "Yes")  # 25 / 195
table(tcga_lum_a$nodePos == "Positive" & tcga_lum_a$recu == "Yes")  # 63 / 264

table(tcga_lum_b$nodePos)
# Negative Positive 
#      199      276 
table(tcga_lum_b$nodePos == "Negative" & tcga_lum_b$recu == "Yes")  # 25 / 199
table(tcga_lum_b$nodePos == "Positive" & tcga_lum_b$recu == "Yes")  # 65 / 276

table(tcga_basal$nodePos)
# Negative Positive 
#      100       57
table(tcga_basal$nodePos == "Negative" & tcga_basal$recu == "Yes")  # 14 / 100
table(tcga_basal$nodePos == "Positive" & tcga_basal$recu == "Yes")  # 20 / 57


####################           Take LumA and Basal               ##########################
tcga_recu_status_all = cbind(subtype_lum_a = subtype_lum_a, tcga_recu_status)
tcga_recu_status_all = cbind(subtype_basal = subtype_basal, tcga_recu_status_all)
rownames(tcga_recu_status_all) = tcga_recu_status$patient
dim(tcga_recu_status_all)   # 1067*18895.
tcga_recu_status_all[1:3, 1:20]

# Remove outliers seen from PCA.
outliers = c("TCGA-AC-A2QJ", "TCGA-OL-A66P")  # Seen from basal pca.
outlier_row = match(outliers, rownames(tcga_recu_status_all))
tcga_recu_status_all = tcga_recu_status_all[-outlier_row, ]
dim(tcga_recu_status_all)  # 1065*18895

# Remove sample with NA in days.
tcga_recu_status_all = tcga_recu_status_all[!is.na(tcga_recu_status_all$days), ]   # 1064*18895
dim(tcga_recu_status_all)

# add recu1  add recuYes
tcga_recu_status_all = data.frame(recu1 = 0, recuYes = factor("Yes", levels = c("Yes", "No")), tcga_recu_status_all)
tcga_recu_status_all$recu1[tcga_recu_status_all$recu == "Yes"] = 1
table(tcga_recu_status_all$recu1)
tcga_recu_status_all$recuYes[tcga_recu_status_all$recu1 == 0] = factor("No", levels = c("Yes", "No")) 
table(tcga_recu_status_all$recuYes)

# remove microRNA
dim(tcga_recu_status_all)   # 1064*18897
tcga_recu_status_all = tcga_recu_status_all[, !grepl('MIR', colnames(tcga_recu_status_all))]
dim(tcga_recu_status_all)   # 1064*18897  20Jan2020
tcga_recu_status_all[1:2, 1:25]

# Save.
saveRDS(tcga_recu_status_all, "tcga_recu_status_all_outliers_removed_1064_18897_20Jan2020.RDS")


# Read in.
tcga_recu_status_all = readRDS("./tcga_recu_status_all_outliers_removed_1064_18897_20Jan2020.RDS")


######## From below on, use subtype_i to represent subtype, no specific subtype.
########  QC
# Counts.
mat_all = tcga_recu_status_all[, -(1:20)]
dim(mat_all)
mat_all[1:3, 1:3]
recu_status_all = tcga_recu_status_all[, 1:20]
dim(recu_status_all)
recu_status_all[1:3, ]


source("./source_plot_count_pca.R")
plot_count_pca(mat_all, recu_status_all)



