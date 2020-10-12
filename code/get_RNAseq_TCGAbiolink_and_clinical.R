# setwd("/Users/Haipeng/Documents/breast_cancer/RNA-seq/write_paper/data_rnaseq_htseq/")
setwd("/Users/hl/Documents/breast_cancer/write_paper_reDo_3/0.1.1_data__tcga/")
#

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("maftools", version = "3.8")



library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(DT)
library(xlsx)
library(openxlsx)


#
# query = GDCquery(
#   project = "TCGA-BRCA",
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "HTSeq - Counts"
#   # barcode = c("TCGA-BH-A1ET-01A-11R-A137-07", "TCGA-E9-A1NI-01A-11R-A14D-07")
# )
# GDCdownload(query)
# data <- GDCprepare(query)
# saveRDS(data, file = "data_tcga_htseq_sumexp_and_clinical.RDS")
data = readRDS('data_tcga_htseq_sumexp_and_clinical.RDS')            # 2.25.2019


# str(data)
class(data)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"
# 
# mat = assay(data)
# dim(mat)   
# # [1] 56925  1222
# mat[1:3, 1:3]


col_data = colData(data)
dim(col_data)   # 1222 81.
length(col_data$patient)  # 1222
length(unique(col_data$patient))  # 1092 1092 unique patient.
table(col_data$shortLetterCode) 
#  NT   TM   TP   # TP: primary solid tumor. NT: solid tissue normal TM:metastatic
# 113    7 1102   ############# 1102+7=1109 should be used, removing 113 normal. Wrong, TM should not be used.

# keep only TP  and keep only is_ffpe FALSE.
col_data = col_data[(col_data$shortLetterCode %in% "TP") & (col_data$is_ffpe == FALSE), ] 
dim(col_data)  # 1086 81
length(unique(col_data$patient))  # 1081    # still 5 duplicates.
col_data$patient[duplicated(col_data$patient)]
# [1] "TCGA-A7-A26J" "TCGA-A7-A0DB" "TCGA-A7-A13D" "TCGA-A7-A26E" "TCGA-A7-A13E"

#  for same patient with two counts files, remove the one with less counts.
sum(assay(data)[, "TCGA-A7-A13E-01A-11R-A277-07"])
dup_barcode_remove = c("TCGA-A7-A26J-01A-11R-A277-07", "TCGA-A7-A0DB-01A-11R-A277-07",
                       "TCGA-A7-A13D-01A-13R-A277-07", "TCGA-A7-A26E-01A-11R-A277-07",
                       "TCGA-A7-A13E-01A-11R-A277-07")

col_data = col_data[!(col_data$barcode %in% dup_barcode_remove), ]
dim(col_data)   # 1081 81
length(unique(col_data$patient))  # 1081    no more duplicate.

mat_data = assay(data)[, rownames(col_data)] # assay has row genes column sample. col_data has row sample column genes.
dim(mat_data)  # 56925 1081.

saveRDS(mat_data, 'tcga_raw_counts_56925_1081.RDS')

saveRDS(col_data, 'tcag_sumExp_colData_1081_81.RDS')
write.xlsx(as.matrix(col_data), "col_data_clinical_tcga_1081_81.xlsx")

col_data = readRDS("tcag_sumExp_colData_1081_81.RDS")
View(t(as.matrix(col_data)))
View(as.matrix(col_data))
