#' volcano_plot
#'
#' This function was written since the volcanoplot is an easy way to globally assess the distribution of the p values relative to the effect sizes. This specific function has been tailored to my need of working with methylation data. 
#' @param effect_sizes A numeric vector of effect sizes
#' @param pvals A numeric vector of (adjusted) p values)
#' @param identifiers A character vector of identifiers
#' @param int_effect_threshold A vertical dashed line to indicate results above an effect size. Defaults to NULL.
#' @param pval_threshold A horizontal dashed line to indicate the statistically significant results. Defaults to 0.05
#' @param effect_limit If you know there is an upper limit to the effect size. Defaults to finding the limits automatically
#' @param top_names A numeric value representing whether you want to plot the names of the top (based on significance) X hits. Defaults to NULL.
#' @param title Title of the plot.
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A volcano plot with the effect size on the X axis and the -log10 adjusted p-value on the Y axis. 
#' @keywords volcanoplot
#' @export
#' @examples 
#' #Load data
#' require(minfiData)
#' baseDir <- system.file("extdata", package = "minfiData")
#' targets <- read.metharray.sheet(base = baseDir)
#' RGset <- read.metharray.exp(targets = targets, recursive = T)
#' Mset <- preprocessIllumina(rgSet = RGset, bg.correct = T, normalize = "controls", reference = 2)
#' Rset <- ratioConvert(Mset)
#' GMset <- mapToGenome(Rset)
#' 
#' Beta <- getBeta(GMset)
#' design <- model.matrix(~targets$Sample_Group)
#' 
#' #Linear regression
#' lfit <- lmFit(Beta, design)
#' lfit <- eBayes(lfit)
#' top_genes <- toptable(fit = lfit, coef = 2, number = Inf)
#' 
#' #Plot the volcano_plot
#' 

volcano_plot <- function(effect_sizes, pvals, identifiers, int_effect_threshold = NULL, pval_threshold = 0.05, effect_limit = NULL, top_names = NULL, title){

  #Sanitize the data
  if(is.null(effect_sizes)) stop("No vector of effect sizes defined")
  else{effect_sizes <- as.numeric(effect_sizes)}
  if(is.null(pvals)) stop("No vector of p-values provided")
  else{pvals <- as.numeric(pvals)}
  if(is.null(identifiers)) stop("No identifiers provided")
  else{identifiers <- as.character(identifiers)}
  if(!is.null(top_names) & !is.numeric(top_names)) stop("The top_names must be a numeric value") 
  
  if(length(effect_sizes) != length(pvals)) stop("Number of effect sizes and the number of p-values does not correspond")
  if(length(effect_sizes) != length(identifiers)) stop("Number of effect sizes and the number of identifiers does not correspond")
  
  #Generate data frame for ggplot
  volplot <- data.frame(effect_sizes, pvals, identifiers)
  volplot$threshold <- (volplot$pvals < pval_threshold)*1
  volplot$threshold[which(volplot$threshold == 1)] <- "Significant"
  volplot$threshold[which(volplot$threshold == 0)] <- "Non significant"
  if(!is.null(int_effect_threshold)){
    volplot$threshold[abs(volplot$effect_sizes) >= abs(int_effect_threshold)] <- "Interesting"
  }
  volplot$threshold <- factor(volplot$threshold, levels = c("Significant", "Non significant", "Interesting"))  
  
  x_lim <- max(abs(effect_sizes))
  
  drawplot <- ggplot(volplot, aes(x = effect_sizes, y = -log10(pvals), label = identifiers, color = threshold))
  if(!is.null(top_names)){
    drawplot <- drawplot + geom_label_repel(data=head(volplot[order(volplot$pvals, decreasing= F),], n=top_names), show.legend = F)
  }
  drawplot <- drawplot + geom_point(alpha = 0.4, size = 1)
  drawplot <- drawplot + geom_hline(yintercept = -log10(pval_threshold), linetype = "longdash")
  if(!is.null(int_effect_threshold)){
    drawplot <- drawplot + geom_vline(xintercept = c(-int_effect_threshold, int_effect_threshold), linetype = "longdash")
  } 
  drawplot <- drawplot + theme_bw()
  drawplot <- drawplot + xlim(c(-x_lim, x_lim))
  drawplot <- drawplot + ylab("-log10(P)")
  drawplot <- drawplot + xlab("Mean effect size")
  drawplot <- drawplot + scale_color_brewer(palette = "Dark2")
  drawplot <- drawplot + guides(colour = guide_legend(override.aes = list(size = 10)))
  drawplot <- drawplot + theme(axis.text = element_text(size = 12), 
                               axis.title = element_text(size = 14),
                               legend.title = element_blank(),
                               legend.position = "bottom")
  drawplot
}