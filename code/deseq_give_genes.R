deseq_give_genes = function(mat_input, recu_input, resample_n, deseq_p, deseq_f) {
  sub_mat=mat_input
  sub_recu=recu_input
  
  ##############   DESeq2.
  dds = DESeqDataSetFromMatrix(
    countData = t(sub_mat),      # input raw count.
    colData = sub_recu,
    design = ~ recuYes
  )
  set.seed(123)
  dds = DESeq(dds, parallel = TRUE, fitType = 'mean')
  
  res1 = results(dds, contrast = c("recuYes", "No", "Yes"))
  
  res_plot = data.frame(l2fc=res1$log2FoldChange, padj=res1$padj)
  p_cutoff = deseq_p
  f_cutoff = deseq_f
  
  cat("p_cutoff : ", p_cutoff, '\n')
  cat("f_cutoff : ", f_cutoff, '\n')
  
  res_plot$sig_status = "notSig"
  res_plot$sig_status[abs(res_plot$l2fc) > f_cutoff & res_plot$padj < p_cutoff] = "Sig"
  res_plot$padj_neg_log10 = -log10(res_plot$padj)
  
  vol_plot = ggplot(res_plot, aes(l2fc, padj_neg_log10, colour=sig_status)) +
    geom_point(size=2) +
    geom_vline(xintercept = c(f_cutoff, -f_cutoff), linetype=2) +
    geom_hline(yintercept = -log10(p_cutoff), linetype=2) +
    theme_bw(base_size = 22) +
    theme(legend.position = "none") +
    labs(x=bquote(~Log[2]~ "fold change"), y=bquote(~-Log[10]~italic(Padj))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
    scale_color_manual(values = c("darkgrey", "red"))
  
  
  resSig <- subset(res1, padj < p_cutoff & abs(log2FoldChange) > f_cutoff)  # padj can be changed.
  
  genes_deseq = rownames(resSig)
  
  # ggsave(paste0("volcano_selected_", length(genes_deseq), "_resample_", resample_n, ".pdf"), vol_plot, height = 10, width = 10)
  # saveRDS(vol_plot, paste0("volcano_selected_", length(genes_deseq), "_resample_", resample_n, ".RDS"))
  
  # Heatmap for deseq genes.
  # genes_deseq = c("HMGA2","SOCS2","BMP6","MN1")
  # give_heat_map(sub_mat, sub_recu, genes_heat = genes_deseq, genes_type = 'deseq')
  
  return(genes_deseq)
}
