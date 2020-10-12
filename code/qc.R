setwd('/Users/hl/Documents/breast_cancer/write_paper_reDo_4/0.3_add_ER_node_QC/')


library(ggplot2)
library(DESeq2)
library(ggfortify)
library(ggpubr)
library(grid)
library(gridExtra)

source("./source_read_in_dataset.R")

data_type = "g5"  # GSE96058 SCAN-B.
data_set = read_in_dataset(data_type)
dim(data_set$data_set_mat)  # 3069 18877.
dim(data_set$data_set_recu)  # 3069 24.

#
summary(data_set$data_set_recu)

#       id                   age.at.diagnosis          tumor.size                         lymph.node.status         er.status             pgr.status  
# F1     :   1   age at diagnosis: 69: 122    tumor size: 15: 231   lymph node status: NA          :  96    er status: 0 : 224   pgr status: 0 : 362  
# F10    :   1   age at diagnosis: 67: 113    tumor size: 12: 185   lymph node status: NodeNegative:1878    er status: 1 :2646   pgr status: 1 :2377  
# F100   :   1   age at diagnosis: 66: 110    tumor size: 20: 182   lymph node status: NodePositive:1095    er status: NA: 199   pgr status: NA: 330  
# F1000  :   1   age at diagnosis: 68: 106    tumor size: 11: 169                                                                                     
# F1001  :   1   age at diagnosis: 65: 101    tumor size: 14: 166                                                                                     
# F1002  :   1   age at diagnosis: 63:  95    tumor size: 18: 149                                                                                     
# (Other):3063   (Other)             :2422    (Other)       :1987        

#          her2.status            ki67.status        nhg                     pam50.subtype                  overall.survival.days               overall.survival.event
# her2 status: 0 :2572   ki67 status: 0 : 574   nhg: G1: 454   pam50 subtype: Basal : 325   overall survival days: 1212:  14      overall survival event: 0:2747      
# her2 status: 1 : 392   ki67 status: 1 : 813   nhg: G2:1439   pam50 subtype: Her2  : 307   overall survival days: 2325:  11      overall survival event: 1: 322      
# her2 status: NA: 105   ki67 status: NA:1682   nhg: G3:1115   pam50 subtype: LumA  :1540   overall survival days: 1543:  10                                          
#                                               nhg: NA:  61   pam50 subtype: LumB  : 695   overall survival days: 1640:  10                                          
#                                                              pam50 subtype: Normal: 202   overall survival days: 1900:  10                                          
#                                                                                           overall survival days: 1912:  10                                          
#                                                                                           (Other)                    :3004                                      

#            endocrine.treated           chemo.treated   er_status          pr_status         her2_status        subtype_basal   subtype_lum_a       recu1       
# endocrine treated: 0 : 681    chemo treated: 0 :1807   Length:3069        Length:3069        Length:3069        Mode :logical   Mode :logical   Min.   :0.0000  
# endocrine treated: 1 :2366    chemo treated: 1 :1241   Class :character   Class :character   Class :character   FALSE:2721      FALSE:1232      1st Qu.:0.0000  
# endocrine treated: NA:  22    chemo treated: NA:  21   Mode  :character   Mode  :character   Mode  :character   TRUE :348       TRUE :1837      Median :0.0000  
#                                                                                                                                                 Mean   :0.1049  
#                                                                                                                                                 3rd Qu.:0.0000  
#                                                                                                                                                 Max.   :1.0000  
# 
# recuYes         days        nodePos              years       
# Yes: 322   Min.   :  56   Length:3069        Min.   :0.1533  
# No :2747   1st Qu.:1236   Class :character   1st Qu.:3.3840  
#            Median :1612   Mode  :character   Median :4.4134  
#            Mean   :1606                      Mean   :4.3973  
#            3rd Qu.:2004                      3rd Qu.:5.4867  
#            Max.   :2474                      Max.   :6.7734  


# TCGA
data_type = "tcga"  # GSE96058 SCAN-B.
data_set = read_in_dataset(data_type)
dim(data_set$data_set_mat)  # 1064 18877.
dim(data_set$data_set_recu)  # 1064 20.

#
summary(data_set$data_set_recu)

data_clin = read.csv("../0.2.0_Other_datasets/tcga/1-s2.0-S0092867418302290-mmc1.csv")
data_clin = data_clin[data_clin$type == "BRCA", ]
summary(data_clin)
table(data_clin$PFI, useNA='always')  # progression free interval.
# #N/A    0    1 <NA> 
#    0  952  145    0 
summary(as.numeric(as.character(data_clin$PFI.time)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     0     431     773    1158    1570    8556       1 
