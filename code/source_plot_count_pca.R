##### Source file
##### Plot count and PCA.



plot_count_pca = function(mat_all, recu_status_all) {
  
  
  # dim(mat_all)
  # mat_all[1:3, 1:3]
  # dim(recu_status_all)
  # recu_status_all[1:3, 1:3]
  
  
  subtypes = colnames(recu_status_all)[grep("subtype", colnames(recu_status_all))]
  
  
  for (node_status in c("Negative", "Positive")) {
    
    print(node_status)
    
    for (subtype_i in subtypes) {
      print(subtype_i)
      # subtype_i = subtypes[1]
      
      # mat
      sample_select = 
        recu_status_all[, subtype_i] & recu_status_all$nodePos == node_status & !is.na(recu_status_all$nodePos)
      
      mat_i = mat_all[sample_select, ]
      
      # recu
      recu_i = recu_status_all[sample_select, ]
      recu_i$recu1 = 0
      recu_i$recu1[recu_i$recu == 'Yes'] = 1
      
      total_count = data.frame(sample = rownames(mat_i), totalCount = rowSums(mat_i), recu1 = recu_i$recu1)

      # count plot
      total_count = total_count[order(total_count$totalCount), ]
      
      total_count$reorder = with(data = total_count, reorder(sample, totalCount,  function(x) x))
      
      fill_type = 'recu1' # er recu Type er recu node pr her2
      total_count$recu1 = as.factor(total_count$recu1)
      total_count$reorder2 = seq(1, dim(total_count)[1])
      
      pl = ggplot(total_count, aes(reorder2, totalCount, fill = total_count[, fill_type])) +  
        geom_bar(stat = 'identity') +
        geom_hline(yintercept =  20000000, color = 'green', size = 1.5, linetype = 'dashed') +  # 20 Million
        
        scale_x_continuous(breaks = pretty(seq(0, 1000, 21), 21)) +
        scale_y_continuous(breaks = pretty(seq(0, 130000000), 14), 
                           labels = pretty(seq(0, 130), 14)) +
        
        scale_fill_discrete(name = "Recurrence", labels = c("No", "Yes")) +
        
        theme(axis.text.x = element_text(angle = 90, hjust = -0.3)) +
        # grids(linetype = 'dashed') +
        border(color = 'black', size = 2, linetype = 'solid') +
        theme_bw()  +
        # theme_light() +
        # theme_pubr() +
        theme(legend.position = 'right') 
      
      #
      pl = ggpar(
        pl,
        xlab = "Sample",
        ylab = 'Sequencing counts (million)',
        font.x = c(25, "plain", 'black'),
        font.y = c(25, 'plain', 'black'),
        font.tickslab = c(20, 'plain', 'black'),
        font.legend = c(20, 'plain', 'black'),
        
        # legend.title = "Recurrence",
        ticks = TRUE,
        tickslab = TRUE
      )
      
      save_name = paste("Plot_count", subtype_i, "node", node_status,'.pdf', sep ='_')
      ggsave(pl, file = save_name, width = 15, height = 10)
      saveRDS(pl, paste0(save_name, ".RDS"))
      
      
      # PCA
      # vst.
      vsd = vst(t(mat_i), blind = FALSE)   # vst needs row gene colomn sample.
      vsd = t(vsd)
      recu_i$recu1 = as.factor(recu_i$recu1)
      
      color_by = 'recu1' 
      
      prin_comp = prcomp(vsd)
      data = recu_i
      
      colour = color_by
      
      pl2 = autoplot(
        prin_comp,
        data = data,
        colour = colour,
        label = F
      )
      
      pl2 = pl2 +
        scale_color_discrete(name = 'Recurrence', labels = c("No", "Yes")) +
        theme(plot.margin = margin(1.1, 1.1, 1.1, 1.1, 'cm')) +
        # theme(legend.position = ifelse(plot_legend, "right", 'none')) 
        # geom_point() +
        # geom_text(aes(label=patient), nudge_x = 0, nudge_y = 0) +
        theme(legend.position = "right") 
      
      
      pl2 = ggpar(
        pl2,
        title = subtype_i,
        # xlab = "Sample",
        # ylab = 'Sequencing counts (million)',
        font.x = c(15, "plain", 'black'),
        font.y = c(15, 'plain', 'black'),
        font.tickslab = c(12, 'plain', 'black'),
        font.legend = c(15, 'plain', 'black')
        
        # legend.title = "Recurrence",
        # ticks = TRUE,
        # tickslab = TRUE
      )
      pl2
      
      save_name = paste("plot_pca", subtype_i,"node", node_status, ".pdf", sep = "_")
      ggsave(pl2, file = save_name, width = 15, height = 10, dpi = 400, scale = 0.8)
      saveRDS(pl2, paste0(save_name, '.RDS'))
      
    }
    
    
  }
  
}

