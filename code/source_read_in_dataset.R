# Read in data

read_in_dataset = function(data_type) {
  if (data_type == "tcga") {
    # tcga
    data_set = readRDS("../0.3_add_ER_node_QC/tcga_recu_status_all_outliers_removed_1064_18897_20Jan2020.RDS")
    data_set_mat = data_set[, -c(1:20)]
    # dim(data_set_mat)  # 1023 36862  # 1064*18877 (20Jan2020)
    data_set_recu = data_set[, 1:20]
  } else if (data_type == 'gh') {
    # gh
    data_set = readRDS("../0.2.1_commen_genes_tcga_gh_others/mat_recu_gh.RDS")
    data_set_mat = data_set[, -c(1:11)]
    data_set_recu = data_set[, 1:11]
  } else if (data_type == 'g3') {
    # g3   NO subtype_basal.
    data_set = readRDS("../0.2.1_commen_genes_tcga_gh_others/mat_recu_g3.RDS")
    data_set_mat = data_set[, -c(1:20)]
    data_set_recu = data_set[, 1:20]
    print(table(data_set_recu$nodePos))
    print(table(data_set_recu$subtype_basal))
    print(table(data_set_recu$subtype_lum_a))
    print(table(data_set_recu$her2_status))
    data_set_recu$years = data_set_recu$days/365.25
  } else if (data_type == 'g5') {
    # g5
    data_set = readRDS("../0.2.1_commen_genes_tcga_gh_others/mat_recu_g5.RDS")
    data_set_mat = data_set[, -c(1:23)]
    data_set_recu = data_set[, 1:23]
    data_set_recu$years = data_set_recu$days/365.25
  }
  
  return(list(data_set_mat=data_set_mat, data_set_recu=data_set_recu))
}