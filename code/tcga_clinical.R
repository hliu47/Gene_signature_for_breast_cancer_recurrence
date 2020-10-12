# setwd("/Users/Haipeng/Documents/breast_cancer/RNA-seq/write_paper/clinical_info/")
setwd("/Users/hl/Documents/breast_cancer/write_paper_reDo_3/0.1.1_data__tcga/")

library(TCGAbiolinks)
library(DT)
library(dplyr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

# Show all project
TCGAbiolinks:::getGDCprojects()$project_id
# project 
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
# $data_categories
# case_count file_count               data_category
# 1       1097       6080     Transcriptome Profiling
# 2       1098       4459       Copy Number Variation
# 3       1044       8648 Simple Nucleotide Variation
# 4       1095       1234             DNA Methylation
# 5       1098       1183                    Clinical
# 6       1098       4604            Sequencing Reads
# 7       1098       5316                 Biospecimen
# 
# $case_count
# [1] 1098
# 
# $file_count
# [1] 31524        31598 (1.17.2020)
# 
# $file_size
# [1] 5.550475e+13




# # indexed clinical data same as col_data from SumExp.
# clinical = GDCquery_clinic(project = "TCGA-BRCA")  # 2.25.2019  # same as col_data from SumExp.
# #
# DT::datatable(
#   clinical,
#   filter = 'top',
#   options = list(scrollX = TRUE, keys = TRUE, pageLength = 10),
#   rownames = FALSE
# )
# 
# summary(clinical$age_at_diagnosis / 365)
# length(unique( clinical$bcr_patient_barcode))  # 1097 unique
# table(clinical$progression_or_recurrence)
# # not reported 
# #         1097


# Raw XML clinical 
# clinical query
# XML 
# query = GDCquery(                 # 2.25.2019
#   project = "TCGA-BRCA",
#   data.category = "Clinical",
#   file.type = "xml"
#   )
# saveRDS(query, "query_tcga_brca_clinical_xml.RDS")          # 2.25.2019
# query = readRDS('query_tcga_brca_clinical_xml.RDS')
# GDCdownload(query)

# xml_clinical_patient = GDCprepare_clinic(query, clinical.info = "patient")
# saveRDS(xml_clinical_patient, 'xml_clinical_patient.RDS')
xml_clinical_patient = readRDS('xml_clinical_patient.RDS')
xml_clinical_patient = apply(xml_clinical_patient, 2, function(x) na_if(x, ""))    # convert "" to NA.
xml_clinical_patient = as.data.frame(xml_clinical_patient)
# Remove duplicate
xml_clinical_patient = xml_clinical_patient[!duplicated(xml_clinical_patient$bcr_patient_barcode), ]
rownames(xml_clinical_patient) = 1:dim(xml_clinical_patient)[1]
# write.csv(xml_clinical_patient, 'xml_clinical_patient.csv')


# View(t(xml_clinical_patient))

# clin_patient = GDCprepare_clinic(query, clinical.info = "patient")

# clin_new_tumor_event = GDCprepare_clinic(query, clinical.info = "new_tumor_event")
# saveRDS(clin_new_tumor_event, "clin_new_tumor_event.RDS")
clin_new_tumor_event = readRDS("clin_new_tumor_event.RDS")
clin_new_tumor_event = apply(clin_new_tumor_event, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_new_tumor_event = as.data.frame(clin_new_tumor_event)


# clin_drug = GDCprepare_clinic(query, clinical.info = "drug")
# saveRDS(clin_drug, "clin_drug.RDS")
clin_drug = readRDS('clin_drug.RDS')
clin_drug = apply(clin_drug, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_drug = as.data.frame(clin_drug)

# clin_follow_up = GDCprepare_clinic(query, clinical.info = "follow_up")
# clin_follow_up = saveRDS(clin_follow_up, "clin_follow_up.RDS")
clin_follow_up = readRDS("clin_follow_up.RDS")
clin_follow_up = apply(clin_follow_up, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_follow_up = as.data.frame(clin_follow_up)

# clin_radiation = GDCprepare_clinic(query, clinical.info = "radiation")
# saveRDS(clin_radiation, "clin_radiation.RDS")
clin_radiation = readRDS("clin_radiation.RDS")
clin_radiation = apply(clin_radiation, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_radiation = as.data.frame(clin_radiation)

# clin_admin = GDCprepare_clinic(query, clinical.info = "admin")
# saveRDS(clin_admin, "clin_admin.RDS")
clin_admin = readRDS("clin_admin.RDS")
clin_admin = apply(clin_admin, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_admin = as.data.frame(clin_admin)

# clin_stage_event = GDCprepare_clinic(query, clinical.info = "stage_event")
# saveRDS(clin_stage_event, "clin_stage_event.RDS")
clin_stage_event = readRDS("clin_stage_event.RDS")
clin_stage_event = apply(clin_stage_event, 2, function(x) na_if(x, ""))    # convert "" to NA.
clin_stage_event = as.data.frame(clin_stage_event)




##########################################
# Using xml_clinical_patient.
# PR ER HER2 NODE.

clin_my = NULL

# PR
clin_my$pr = xml_clinical_patient$breast_carcinoma_progesterone_receptor_status
# clin_my$pr_metastasis = xml_clinical_patient$metastatic_breast_carcinoma_progesterone_receptor_status

# ER
clin_my$er = xml_clinical_patient$breast_carcinoma_estrogen_receptor_status

# HER2
clin_my$her2_ihc = xml_clinical_patient$lab_proc_her2_neu_immunohistochemistry_receptor_status
clin_my$her2_ish = xml_clinical_patient$lab_procedure_her2_neu_in_situ_hybrid_outcome_type

# Node
clin_my$node_pos_by_he = as.numeric(as.character(xml_clinical_patient$number_of_lymphnodes_positive_by_he))
clin_my$node_pos_by_ihc = as.numeric(as.character(xml_clinical_patient$number_of_lymphnodes_positive_by_ihc))

# clin_my$node_exam_count = xml_clinical_patient$lymph_node_examined_count   # may just be how many node examed, not positive ones.

# meatastasis status
clin_my$metastasis_distant_present_ind2 = xml_clinical_patient$distant_metastasis_present_ind2  # direct distance metastasis indicator.
clin_my$metastasis_at_diagnosis_stageIV = xml_clinical_patient$stage_event_pathologic_stage     # stage IV is metastasis.
clin_my$metastasis_at_diagnosis_TNM1 = xml_clinical_patient$stage_event_tnm_categories          # TNM M1 is metastasis.

# time to last followup or death.  
# this needs to be updated, if any, with respect to followup and new tumor event.
clin_my$days_to_last_followup = xml_clinical_patient$days_to_last_followup
clin_my$days_to_death = xml_clinical_patient$days_to_death

# patient barcode
clin_my$patient_barcode = xml_clinical_patient$bcr_patient_barcode

# this needs to be updated, if any, with respect to followup and new tumor event.
clin_my$vital_status = xml_clinical_patient$vital_status 

clin_my$age_days = xml_clinical_patient$days_to_birth
clin_my$age_at_diagnosis = xml_clinical_patient$age_at_initial_pathologic_diagnosis
clin_my$gender = xml_clinical_patient$gender
clin_my$race = xml_clinical_patient$race_list

clin_my$radiation_therapy = xml_clinical_patient$radiation_therapy                 # 80 YES    108 NO.
clin_my$has_drug_info = xml_clinical_patient$has_drugs_information
clin_my$has_radiation_info = xml_clinical_patient$has_radiations_information
clin_my$has_followup_info = xml_clinical_patient$has_follow_ups_information
clin_my$has_new_tumor_event_info = xml_clinical_patient$has_new_tumor_events_information

# history_of_neoadjuvant_treatment   important should remove any YES. 
# only need naive cancer patients without previous neoadjuvent treatment.
clin_my$history_of_neoadjuvant_treatment = xml_clinical_patient$history_of_neoadjuvant_treatment  # show be no for naive cancer.

#
clin_my = as.data.frame(clin_my)

for (i in 1:dim(clin_my)[1]) {
  
  # i = 1
  
  # for node
  if (!is.na(clin_my$node_pos_by_he[i]) & !is.na(clin_my$node_pos_by_ihc[i])) {
    clin_my$node_my[i] = clin_my$node_pos_by_he[i] + clin_my$node_pos_by_ihc[i]

  } else if (!is.na(clin_my$node_pos_by_he[i]) & is.na(clin_my$node_pos_by_ihc[i])) {
    clin_my$node_my[i] = clin_my$node_pos_by_he[i]

  } else if (is.na(clin_my$node_pos_by_he[i]) & !is.na(clin_my$node_pos_by_ihc[i])) {
    clin_my$node_my[i] = clin_my$node_pos_by_ihc[i]

  } else if (is.na(clin_my$node_pos_by_he[i]) & is.na(clin_my$node_pos_by_ihc[i])) {
    clin_my$node_my[i] = NA

  }
  
  
  # for HER2
  temp = c(as.character(clin_my$her2_ihc[i]), as.character(clin_my$her2_ish[i]))
  
  if ("Positive" %in% temp) {
    clin_my$her2_my[i] = "Positive"
  } else if (!("Positive" %in% temp) & "Negative" %in% temp) {
    clin_my$her2_my[i] = "Negative"
  } else {
    clin_my$her2_my[i] = NA
  }
 
  
  
  # Metastasis
  if (grepl('IV', clin_my$metastasis_at_diagnosis_stageIV[i]) | grepl('M1', clin_my$metastasis_at_diagnosis_TNM1[i])) {
    clin_my$metastasis_at_diagnosis_my[i] = "Yes"
  } else if (grepl('M0', clin_my$metastasis_at_diagnosis_TNM1[i])) {
    clin_my$metastasis_at_diagnosis_my[i] = "No"
  } else if (!(clin_my$metastasis_at_diagnosis_stageIV[i] %in% c("", "Stage IV", "Stage X"))) {
    clin_my$metastasis_at_diagnosis_my[i] = "No" 
  } else {
    clin_my$metastasis_at_diagnosis_my[i] = NA
  }
    
    
  
  
  # # time to last followup or death. 
  # # need to be updated with regards to followup and new tumor event.
  # if (clin_my$vital_status[i] == 'Dead') {
  #   clin_my$days_to_death_or_lastfollowup[i] =  clin_my$days_to_death[i]
  # } else if (clin_my$vital_status[i] == 'Alive') {
  #   clin_my$days_to_death_or_lastfollowup[i] = clin_my$days_to_last_followup[i]
  # } else {
  #   clin_my$days_to_death_or_lastfollowup[i] = NA
  # }
  
    
}


                                                                                  # Below updated on 1.17.2020.
table(clin_my$er, useNA = "always")
# Indeterminate      Negative      Positive          <NA> 
#             2           238           808            49 

table(clin_my$pr, useNA = "always")
# Indeterminate      Negative      Positive          <NA> 
#             4           344           699            50 

table(clin_my$her2_my, useNA = "always")
# Negative Positive     <NA> 
#      760      197      140 

table(clin_my$node_my > 0, useNA = "always")
# FALSE  TRUE  <NA> 
#   427   508   162 

table(clin_my$metastasis_distant_present_ind2, useNA = "always")
#    NO  YES <NA> 
#   381   12  704 

table(clin_my$metastasis_at_diagnosis_my, useNA = "always")
#     No  Yes <NA> 
#   1065   22   10 
##########################################


##########################################
# Using clin_new_tumor_event
new_tumor_event_my = data.frame(
  bcr_patient_barcode = clin_new_tumor_event$bcr_patient_barcode,
  days_to_new_tumor_event_after_initial_treatment = clin_new_tumor_event$days_to_new_tumor_event_after_initial_treatment,
  new_neoplasm_event_type = clin_new_tumor_event$new_neoplasm_event_type,
  new_neoplasm_event_occurrence_anatomic_site = clin_new_tumor_event$new_neoplasm_event_occurrence_anatomic_site,
  new_neoplasm_occurrence_anatomic_site_text = clin_new_tumor_event$new_neoplasm_occurrence_anatomic_site_text,
  new_tumor_event_additional_surgery_procedure = clin_new_tumor_event$new_tumor_event_additional_surgery_procedure,
  days_to_new_tumor_event_additional_surgery_procedure = clin_new_tumor_event$days_to_new_tumor_event_additional_surgery_procedure
)

new_tumor_event_my$days_to_recu_or_metastasis_processed = new_tumor_event_my$days_to_new_tumor_event_after_initial_treatment

for(i in 1:dim(new_tumor_event_my)[1]) {
  if (is.na(new_tumor_event_my$days_to_recu_or_metastasis_processed[i])) {
    if (!is.na(new_tumor_event_my$new_neoplasm_event_type[i]) | !is.na(new_tumor_event_my$new_neoplasm_event_occurrence_anatomic_site[i]) | 
        !is.na(new_tumor_event_my$new_neoplasm_occurrence_anatomic_site_text[i]) | ("YES" %in% new_tumor_event_my$new_tumor_event_additional_surgery_procedure[i])) {
      new_tumor_event_my$days_to_recu_or_metastasis_processed[i] = new_tumor_event_my$days_to_new_tumor_event_additional_surgery_procedure[i]
    }
    
  }
}





#########################################
# Using clin_follow_up
follow_up_my = data.frame(
  bcr_patient_barcode = clin_follow_up$bcr_followup_barcode,
  vital_status = clin_follow_up$vital_status,
  days_to_last_followup = as.numeric(as.character(clin_follow_up$days_to_last_followup)),
  days_to_death = as.numeric(as.character(clin_follow_up$days_to_death)),
  days_to_last_known_alive = as.numeric(as.character(clin_follow_up$days_to_last_known_alive)),
  new_tumor_event_after_initial_treatment = clin_follow_up$new_tumor_event_after_initial_treatment,
  days_to_new_tumor_event_after_initial_treatment = as.numeric(as.character(clin_follow_up$days_to_new_tumor_event_after_initial_treatment)),
  additional_surgery_locoregional_procedure = clin_follow_up$additional_surgery_locoregional_procedure,
  days_to_additional_surgery_locoregional_procedure = as.numeric(as.character(clin_follow_up$days_to_additional_surgery_locoregional_procedure)),
  additional_surgery_metastatic_procedure = clin_follow_up$additional_surgery_metastatic_procedure,
  days_to_additional_surgery_metastatic_procedure = as.numeric(as.character(clin_follow_up$days_to_additional_surgery_metastatic_procedure)),
  recu_meta = NA,
  days_to_recu_meta = NA
)

# update days_to_last_followup
for (i in 1:dim(follow_up_my)[1]) {
  if (!is.na(follow_up_my$days_to_last_known_alive[i])) {
    print(i) # only two.
    if (!is.na(follow_up_my$days_to_last_followup[i])) {
       if (as.numeric(as.character(follow_up_my$days_to_last_followup[i])) < as.numeric(as.character(follow_up_my$days_to_last_known_alive[i]))) {
      follow_up_my$days_to_last_followup[i] = follow_up_my$days_to_last_known_alive[i] }
    } else {
      follow_up_my$days_to_last_followup[i] = follow_up_my$days_to_last_known_alive[i]
    }
  }
}




for (i in 1:dim(follow_up_my)[1]) {
  
  # recu_meta
  if ("YES" %in% follow_up_my$new_tumor_event_after_initial_treatment[i] |
      !is.na(follow_up_my$days_to_new_tumor_event_after_initial_treatment[i]) |
      "YES" %in% follow_up_my$additional_surgery_locoregional_procedure[i] |
      "YES" %in% follow_up_my$additional_surgery_metastatic_procedure[i] |
      !is.na(follow_up_my$days_to_additional_surgery_locoregional_procedure[i]) |
      !is.na(follow_up_my$days_to_additional_surgery_metastatic_procedure[i])) {
    follow_up_my$recu_meta[i] = "Yes"
  }
  
  # days
  if (!is.na(follow_up_my$days_to_new_tumor_event_after_initial_treatment[i])) {
    follow_up_my$days_to_recu_meta[i] = follow_up_my$days_to_new_tumor_event_after_initial_treatment[i]
  } else if (!is.na(follow_up_my$days_to_additional_surgery_locoregional_procedure[i])) {
    follow_up_my$days_to_recu_meta[i] = follow_up_my$days_to_additional_surgery_locoregional_procedure[i]
  } else if (!is.na(follow_up_my$days_to_additional_surgery_metastatic_procedure[i])) {
    follow_up_my$days_to_recu_meta[i] = follow_up_my$days_to_additional_surgery_metastatic_procedure[i]
  }
  
  
}


##########################################
# Using clin_radiation
# radiation should be administered after diagnosis. after tissue sample collection.
# breast cancer diagnosis, then colecte tissue sample, then do lab test, then decide treatment plan, chemotherapy, 
# or hormonal therapy or radiotherapy.   
# if previous had radiation or chemo then collect tissue sample, it's not ok, their samples should be removed
# should keep chemo, radiation naive tissue samples for RNA-Seq analysis.
  
#########################################
all_status_my = data.frame(
  stringsAsFactors = F,
  
  # bcr_patient_barcode = unique(xml_clinical_patient$bcr_patient_barcode),
  bcr_patient_barcode = xml_clinical_patient$bcr_patient_barcode,
  
  # final status.
  recu_or_meta_or_dead = NA,
  
  # days to final status, 
  days_to_recu_meta_dead_last_followup = NA, 
  
  # age at first diagonosis.
  age_first_diagnosis = NA,
  
  # these er pr her2 node should be all at the time of first diagnosis, 
  # when RNA-Seq samples were collected.
  # anyone with history_of_neoadjuvent_treatment should be removed, for that would effected gene expression. 
  er = NA,
  pr = NA,
  her2 = NA,
  node = NA,
  
  # history of neoadjuvent treatment.
  history_of_neoadjuvant_treatment = NA,
  
  # need to update with regarding to followup and new tumor event.
  vital_status = NA,
  days_to_death = NA,
  days_to_last_followup = NA,
  
  # dist meta
  distant_metastasis = NA,
  metastasis_at_diagnosis = NA,

  days_to_distant_metastasis = NA,
  
  #
  locoregional_recurrence = NA,
  days_to_locoregional_recurrence = NA,
  
  # from new_tumor_event_my
  days_to_loco_recu_or_dist_meta = NA,
  
  # from follow_up_my
  follow_up_recu_meta = NA,
  
  # these should be treatment after initial diagnosis.
  tamoxifen = NA,  # homonal therapy.
  chemotherapy = NA,
  radiation_therapy = NA
)

# First step, fill it in with clin_my
# Second step, fill in with new_tumor_event_my,
# Third step, update with follow_up_my.

# First step, fill it in with clin_my
for (i in 1:dim(all_status_my)[1]) {
  
  # print(i)
  all_status_my$age_first_diagnosis[i] = as.character(clin_my$age_at_diagnosis[i])
  
  all_status_my$er[i] = as.character(clin_my$er[i])
  all_status_my$pr[i] = as.character(clin_my$pr[i])
  all_status_my$her2[i] = as.character(clin_my$her2_my[i])
  all_status_my$node[i] = clin_my$node_my[i]
  
  all_status_my$history_of_neoadjuvant_treatment[i] = as.character(clin_my$history_of_neoadjuvant_treatment[i])
  
  all_status_my$vital_status[i] = as.character(clin_my$vital_status[i])
  all_status_my$days_to_death[i] = as.numeric(as.character(clin_my$days_to_death[i]))
  all_status_my$days_to_last_followup[i] = as.numeric(as.character(clin_my$days_to_last_followup[i]))
  
  all_status_my$distant_metastasis[i] = as.character(clin_my$metastasis_distant_present_ind2[i])
  all_status_my$metastasis_at_diagnosis[i] = as.character(clin_my$metastasis_at_diagnosis_my[i])

}

# Second step, fill in with new_tumor_event_my
for (i in 1:dim(all_status_my)[1]) {
  temp = as.numeric(as.character(new_tumor_event_my$days_to_recu_or_metastasis_processed[grep(all_status_my$bcr_patient_barcode[i], new_tumor_event_my$bcr_patient_barcode)]))
  temp = temp[!is.na(temp)]
  
  if (length(temp) > 0) {
    temp = min(temp)
    all_status_my$days_to_loco_recu_or_dist_meta[i] = temp
  }
}

# Third step, update with follow_up_my

for (i in 1:dim(all_status_my)[1]) {
  # all_status_my$vital_status[i]
  
  # update vital status.
  temp = as.character(follow_up_my$vital_status[grep(all_status_my$bcr_patient_barcode[i], follow_up_my$bcr_patient_barcode)])
  temp = temp[!is.na(temp)]
  
  if (length(temp) > 0) {
    if('Dead' %in% temp) {
      all_status_my$vital_status[i] = "Dead"
      
      temp2 = as.numeric(as.character(follow_up_my$days_to_death[grep(all_status_my$bcr_patient_barcode[i], follow_up_my$bcr_patient_barcode)]))
      temp2 = temp2[!is.na(temp2)] 
      if(length(temp2) > 0) {
        temp2 = min(temp2)
        all_status_my$days_to_death[i] = temp2
      }
      
    } else {   # alive.
      temp3 = as.numeric(as.character(follow_up_my$days_to_last_followup[grep(all_status_my$bcr_patient_barcode[i], follow_up_my$bcr_patient_barcode)]))
      temp3 = temp3[!is.na(temp3)]
      if (length(temp3) > 0) {
        temp3 = max(temp3)
        
        temp4 = as.numeric(as.character(all_status_my$days_to_last_followup[i]))
        
        if (!is.na(temp4)) {
          if (temp3 >= temp4) {         # only when days_to_last_follow is larger then update.
            all_status_my$days_to_last_followup[i] = temp3 
            # cat('temp3 >= temp4 good updated', as.character(all_status_my$bcr_patient_barcode[i]), "\n")
          } else if (temp3 < temp4){
            cat('temp3 < temp4 not good', as.character(all_status_my$bcr_patient_barcode[i]), "\n")
          }
        } else if (is.na(temp4)) {
          all_status_my$days_to_last_followup[i] = temp3 
          cat('temp4 is NA updated', as.character(all_status_my$bcr_patient_barcode[i]), "\n")
        }
      }
    }
  }
  
  # and add follow_up_recu_meta (from follow_up_my).
  temp = as.character(follow_up_my$recu_meta[grep(all_status_my$bcr_patient_barcode[i], follow_up_my$bcr_patient_barcode)])
  temp = temp[!is.na(temp)]
  
  if (length(temp) > 0) {
    if ('Yes' %in% temp) {
      all_status_my$follow_up_recu_meta[i] = "Yes"
    }
  }
  
  # update days_to_loco_recu_or_dist_meta (from new_tumor_event_my)
  temp = as.numeric(as.character(follow_up_my$days_to_recu_meta[grep(all_status_my$bcr_patient_barcode[i], follow_up_my$bcr_patient_barcode)]))
  temp = temp[!is.na(temp)]
  
  if (length(temp) > 0) {
    temp = min(temp)
    temp2 = all_status_my$days_to_loco_recu_or_dist_meta[i]
    if (!is.na(temp2)) {
      if (temp < temp2) {
        all_status_my$days_to_loco_recu_or_dist_meta[i] = temp
        # cat(as.character(all_status_my$bcr_patient_barcode[i]), 'temp < temp2', '\n')
      }
    } else if (is.na(temp2)) {
      all_status_my$days_to_loco_recu_or_dist_meta[i] = temp
      # cat(as.character(all_status_my$bcr_patient_barcode[i]), 'temp2 is NA updated', '\n')
    }
    
  }
  
}


################    finalize process.
# # final recu meta dead status.
for (i in 1:dim(all_status_my)[1]) {
  if ("Dead" %in% all_status_my$vital_status[i] | 
      "YES" %in% all_status_my$distant_metastasis[i] |
      "Yes" %in% all_status_my$metastasis_at_diagnosis[i] | 
      !is.na(all_status_my$days_to_loco_recu_or_dist_meta[i]) |
      "Yes" %in% all_status_my$follow_up_recu_meta[i]) {
    all_status_my$recu_or_meta_or_dead[i] = "Yes"
  }
  
}


# metastasis at diagnosis is stage IV breast cancer, 
# it can be kept for survival prediction as making their time to event as 0,
# for prediction taking only recu status (recu or not_recu) it can be kept too,
# and we make them as recu cases.

# final timing to event.
for (i in 1:dim(all_status_my)[1]) {
  
  if ("Yes" %in% all_status_my$metastasis_at_diagnosis[i]) {
    all_status_my$days_to_recu_meta_dead_last_followup[i] = 0
    
  } else if (!is.na(all_status_my$days_to_loco_recu_or_dist_meta[i])) {
    all_status_my$days_to_recu_meta_dead_last_followup[i] = all_status_my$days_to_loco_recu_or_dist_meta[i]
    
  } else if ("Dead" %in% all_status_my$vital_status[i]) {
    all_status_my$days_to_recu_meta_dead_last_followup[i] = all_status_my$days_to_death[i]
    
  } else if ("Alive" %in% all_status_my$vital_status[i]) {
    all_status_my$days_to_recu_meta_dead_last_followup[i] = all_status_my$days_to_last_followup[i]
  }
}


dim(all_status_my)
write.csv(all_status_my, 'all_status_my_1097_22.csv')


# extract shorter info
all_status_my_short = data.frame(
  patient = as.character(all_status_my$bcr_patient_barcode),
  age = all_status_my$age_first_diagnosis,
  recu = as.character(all_status_my$recu_or_meta_or_dead),
  days = all_status_my$days_to_recu_meta_dead_last_followup,
  # ER = all_status_my$er,
  # PR = all_status_my$pr,
  # HER2 = all_status_my$her2,
  Node = all_status_my$node,
  history_of_neoadjuvant_treatment = all_status_my$history_of_neoadjuvant_treatment
)

# > 0 node is Positive.
all_status_my_short$nodePos = NA
all_status_my_short$nodePos[all_status_my_short$Node > 0] = "Positive"
all_status_my_short$nodePos[all_status_my_short$Node == 0] = "Negative"

# remove that had previous neoadjuvent treatment
all_status_my_short = all_status_my_short[!(all_status_my_short$history_of_neoadjuvant_treatment %in% "Yes"), ]

dim(all_status_my_short)

saveRDS(all_status_my_short, 'all_status_my_short_1084_7_useThis.RDS')
write.csv(all_status_my_short, 'all_status_my_short_1084_7_useThis.csv')





# Biospecimen query.

# query_specimen = GDCquery(              # 2.25.2019
#   project = "TCGA-BRCA",
#   data.category = "Biospecimen",
#   file.type = 'xml'                     # there are images as well.
#   )
# # View(getResults(query_specimen))
# saveRDS(query_specimen, "query_specimen_tcga_brca_biospecimen_xml.RDS")
query_specimen = readRDS('query_specimen_tcga_brca_biospecimen_xml.RDS')
# GDCdownload(query_specimen)

# specimen_sample = GDCprepare_clinic(query_specimen, clinical.info = "sample")
# saveRDS(specimen_sample, "specimen_sample.RDS")
specimen_sample = readRDS("specimen_sample.RDS")

# specimen_bio_patient = GDCprepare_clinic(query_specimen, clinical.info = "bio_patient")
# saveRDS(specimen_bio_patient, "specimen_bio_patient.RDS")
specimen_bio_patient = readRDS("specimen_bio_patient.RDS")

# specimen_portion = GDCprepare_clinic(query_specimen, clinical.info = "portion")
# saveRDS(specimen_portion, "specimen_portion.RDS")
specimen_portion = readRDS("specimen_portion.RDS")

# specimen_protocol = GDCprepare_clinic(query_specimen, clinical.info = "protocol")
# saveRDS(specimen_protocol, "specimen_protocol.RDS")
specimen_protocol = readRDS("specimen_protocol.RDS")

# specimen_aliquot = GDCprepare_clinic(query_specimen, clinical.info = "aliquot")
# saveRDS(specimen_aliquot, "specimen_aliquot.RDS")
specimen_aliquot = readRDS("specimen_aliquot.RDS")

# specimen_analyte = GDCprepare_clinic(query_specimen, clinical.info = "analyte")
# saveRDS(specimen_analyte, "specimen_analyte.RDS")
specimen_analyte = readRDS("specimen_analyte.RDS")

# specimen_slide = GDCprepare_clinic(query_specimen, clinical.info = "slide")
# saveRDS(specimen_slide, "specimen_slide.RDS")
specimen_slide = readRDS("specimen_slide.RDS")

