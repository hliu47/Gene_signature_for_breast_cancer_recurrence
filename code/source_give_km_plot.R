give_km_plot = function(surv_data, plot_legend) {
  km_plot = ggsurvplot(
    surv_data$surv_fit,
    data = predicted_class,
    
    risk.table = TRUE,          # show risk table.
    risk.table.col = "strata",  # Risk table color by groups
    
    pval = TRUE,             # show p-value of log-rank test.
    xlim = c(0,10),        # present narrower X axis, but not affect survival estimates.
    break.time.by = 1,       # break X axis in time intervals by 500.
    surv.scale = "percent",
    
    palette = c('#239B56', '#A04000'),
    
    # xlab = "Year",
    xlab = paste0("Year ", 
                  "HR:", round(surv_data$cox_hr,2), 
                  " 95CI:", "[", round(surv_data$cox_hr_95l,2), ",", round(surv_data$cox_hr_95u,2), "]",
                  " Kappa:", round(surv_data$kappa,2),
                  " Acc:", round(surv_data$acc, 2),
                  " Sens:", round(surv_data$sens, 2),
                  " Spec:", round(surv_data$spec, 2)),
    ylab = "Recurrence-free survival",
    
    legend = ifelse(plot_legend, "right", "none"),
    legend.labs = c(paste0('No-recurrence, predicted'), paste0('Recurrence, predicted')),
    font.legend = c(12, 'bold', 'black'),
    
    ggtheme = theme_light(),
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE, # show bars instead of names in text annotations in legend of risk table
    
    risk.table.fontsize = 4,
    font.x = c(11, "plain", 'black'),
    font.y = c(11, "plain", 'black'),
    font.tickslab = c(8, 'plain', "black")
  )
  
  
  km_plot$plot = km_plot$plot +
    geom_segment(aes(x = 5, y = 0, xend = 5, yend = 1), linetype = "dashed", size = 0.4) +
    geom_segment(aes(x = 10, y = 0, xend = 10, yend = 1), linetype = "dashed", size = 0.4) +
    # theme(axis.text.x=element_text(size=rel(0.5), angle=90)) + 
    # theme(text = element_text(size=20)) +
    theme(plot.margin = margin(0.7,0.7,0,0.7, "cm"))
  
  km_plot$table = ggpubr::ggpar(
    km_plot$table, 
    legend = 'none',
    font.title = c(14, 'plain', 'black'),
    font.x = c(11, "plain", "black"),
    font.y = c(11, "plain", 'black'),
    font.xtickslab = c(8, 'plain', 'black')
  )
  
  km_plot$table = km_plot$table + 
    theme(plot.margin = margin(0,0.7,0.7,0.7, "cm"))
  
  return(km_plot)
}