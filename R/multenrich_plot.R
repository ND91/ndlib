#' multenrich_plot
#'
#' Produces a boxplot overlayed by a heatmap-like plot for multiple enrichment analyses providing a clear overview of the enrichment scores and p-values when having performed multiple enrichment analyses (i.e. different contrasts).
#' It is advisable to ensure that the number of pathways to be plotted is less than 100 or else the list will become very long.
#' The concept and implementation was derived from Martin Garrido Rodriguez-Cordoba's code (https://mgrcbioinfo.github.io/my_GSEA_plot/) and was modified by me.
#' @param enrich_list List of enrichment analyses, wherein each enrichment analysis requires columns containing the pathway name, the enrichment scores (make this 1 or -1 if no particular enrichment score is provided), and a p-value.
#' @param title A title for the plot
#' @param pathway_col The column name of the pathway columns
#' @param es_col The column name of the enrichment scores. Make this 1 or -1 if no particular enrichment score is provided.
#' @param pval_col The column name of the p-values columns. I advise not to use the adjusted p-values as they tend to flatten out (especially Benjamini-Hochberg adjusted p-values).
#'
#' @author Martin Garrido Rodriguez-Cordoba, modified by Andrew Y.F. Li Yim
#'
#' @return A ggplot object containing of heatmap-like plot of the enriched pathways indicating direction, size of the enrichment and the p-value for all provided enrichment analyses.
#' @keywords enrichment analysis, plot, heatmap
#' @export
#' @import ggplot2

multenrich_plot <- function(enrich_list, title, pathway_col, es_col, pval_col){
  plot_groups <- unlist(lapply(enrich_list, nrow))
  plot_df <- data.frame(do.call(rbind, enrich_list))
  plot_df <- data.frame(Pathway = plot_df[,which(colnames(plot_df) == pathway_col)],
                        Enrichment = plot_df[,which(colnames(plot_df) == es_col)],
                        pval = plot_df[,which(colnames(plot_df) == pval_col)],
                        Comparison = rep(names(plot_groups), plot_groups))

  plot_df$Comparison <- factor(plot_df$Comparison, levels = unique(plot_df$Comparison))

  plot_df$Status <- rep("Upregulated", nrow(plot_df))
  plot_df$Status[plot_df$Enrichment < 0] <- "Downregulated"
  plot_df$Status <- factor(plot_df$Status, levels = c("Upregulated", "Downregulated"))

  #Plotting
  ggplot(data = plot_df, mapping = aes(x = Comparison,
                                       y =  Pathway,
                                       colour = Status,
                                       shape = Status,
                                       fill = Enrichment,
                                       size = -log10(pval))) +
    geom_point() +
    scale_shape_manual(values = c(24, 25)) +
    scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
    scale_color_manual(values = c("coral3","deepskyblue3")) +
    ggtitle(title) +
    theme_bw() +
    guides(size = guide_legend(override.aes = list(shape=17))) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y =  element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"))
}
