#' svd_correlator
#'
#' Correlates the principal components with a given confounder. In short, the right singular matrix is correlated column per column with the confounder.
#' @param evs Either the output of SVD: A list containing three outputs: left singular matrix u, vector d and right singular matrix v, or a vector containing an eigen vector matrix obtained from SVD.
#' @param confounder A vector containing the confounder of interest. Order and length must be the same as the number of columns in the right singular matrix v.
#' @param alpha A value depicting the alpha threshold used for defining significance. If NULL, then no significance will be provided. Defaults to 0.05.
#'
#' @author Andrew Y.F. Li Yim
#'
#' @return A scatterplot and data.frame depicting the correlation of the principal components with the confounders.
#' @keywords PCA, SVD, correlation
#' @export
#' @import ggplot2

svd_correlator <- function(evs, confounder, alpha_threshold = 0.05, padj.method = "BH", title){
  if(is.null(evs)) stop("No evs provided")
  if(is.list(evs)){
    evs <- evs$v
  }
  if(ncol(evs) != length(confounder)) stop("Number of columns in the v matrix is not the same as the number of confounders")
  if(missing(title)) title <- paste0("Correlation: Confounder per PC")

  if(!is.numeric(confounder)){
    confounder <- factor(confounder)
  }

  cor.pval <- apply(X = evs, MARGIN = 2, FUN = function(i){
    fit <- lm(i~confounder)
    pval <- lmp(fit)

    return(c(Correlation = summary(fit)$adj.r.squared, P = pval))
  })

  cor.pval.df <- data.frame(t(cor.pval), padj = p.adjust(p = data.frame(t(cor.pval))$P, method = padj.method), PC = 1:ncol(cor.pval))

  if(!is.null(alpha_threshold)){
    cor.pval.df$Significant <- cor.pval.df$padj < alpha_threshold
    plotgraph <- ggplot(cor.pval.df, aes(x = PC, y = Correlation, ymax = 1, col = Significant)) +
      geom_point(size = 3) +
      theme_bw() +
      ggtitle(title) +
      ylab(expression(R^2)) +
      scale_fill_discrete(drop = FALSE) +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = "bottom")
  } else{
    plotgraph <- ggplot(cor.pval.df, aes(x = PC, y = Correlation, ymax = 1)) +
      geom_point(size = 3) +
      theme_bw() +
      ggtitle(title) +
      ylab(expression(R^2)) +
      scale_fill_discrete(drop = FALSE) +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = "bottom")
  }

  return(list(cor.pval.df, plotgraph))
}
