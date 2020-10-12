give_heat_map = function(sub_mat, sub_recu, genes_heat, genes_type) {
  
  # Re-order by recu1.
  sub_recu = sub_recu[order(sub_recu$recu1), ]
  sub_mat = sub_mat[rownames(sub_recu), ]

  heat_mat = sub_mat[, genes_heat]

  heat_mat = t(varianceStabilizingTransformation(t(heat_mat), fitType = 'mean', blind = FALSE))
  heat_mat = t(heat_mat)

  anno_col = data.frame(recu=sub_recu$recuYes, 
                        node = sub_recu$nodePos)
  rownames(anno_col) = rownames(sub_recu)
  
  for (cluster_col in c(FALSE, TRUE)) {
    
    save_name = paste("heatmap_sample", dim(heat_mat)[2], "gene", dim(heat_mat)[1], 
                      genes_type, 'cluster_col', cluster_col, '.pdf', sep = "_")
    pheatmap(
      heat_mat, 
      kmeans_k = NA, cluster_rows = TRUE, cluster_cols = cluster_col, 
      clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2',
      annotation_col = anno_col,  # col is sample.
      filename = save_name, fontsize = 3)
    
  }

}


