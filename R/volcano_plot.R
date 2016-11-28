#' volcano_plot
#'
#' This function was written since the volcanoplot is an easy way to globally assess the distribution of the p values relative to the effect sizes. This specific function has been tailored to my need of working with methylation data. 
#' @param effect_sizes A numeric vector of effect sizes.
#' @param pvals A numeric vector of p values.
#' @param significance A boolean vector indicating the significant hits (T = significant, F = non-significant).
#' @param identifiers A character vector of identifiers.
#' @param int_effect_threshold A vertical dashed line to indicate results above an effect size. Defaults to NULL.
#' @param effect_limit Upper limit of the plot. Defaults to finding the limits automatically.
#' @param top_names An integer representing the names of the top (based on significance) X hits. Defaults to NULL.
#' @param title Title of the plot.
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A volcano plot with the effect size on the X axis and the -log10 p-value on the Y axis. Adjusted p-values are indicated with a different color.
#' @keywords volcanoplot
#' @export
#' @import ggplot2 
#' @import ggrepel
#' @seealso \pkg{\link{minfi}}
#' @seealso \pkg{\link{minfiData}}
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

volcano_plot <- function(effect_sizes, pvals, significance, identifiers, int_effect_threshold = NULL, effect_limit = NULL, top_names = NULL, title = NULL, x_lab = NULL, y_lab = NULL){

  #Sanitize the data
  if(is.null(effect_sizes)) stop("No vector of effect sizes defined")
  else{effect_sizes <- as.numeric(effect_sizes)}
  if(is.null(pvals)) stop("No vector of p-values provided")
  else{pvals <- as.numeric(pvals)}
  if(is.null(significance)) stop("No vector of significance provided")
  else{significance <- as.logical(significance)}
  if(is.null(identifiers)) stop("No vector of identifiers provided")
  else{identifiers <- as.character(identifiers)}
  if(!is.null(top_names) & !is.numeric(top_names)) stop("The top_names must be a numeric value") 
  
  if(length(effect_sizes) != length(pvals)) stop("Size of the effect size vector and the p-value vector do not correspond")
  if(length(effect_sizes) != length(identifiers)) stop("Size of the effect size vector and the identifier vector do not correspond")
  if(length(effect_sizes) != length(significance)) stop("Size of the effect size vector and the significance vector do not correspond")
  
  #which(is.na(significance))
  
  #Generate data frame for ggplot
  volplot <- data.frame(effect_sizes, pvals, significance, identifiers)
  volplot$threshold[which(volplot$significance == T)] <- "significant"
  volplot$threshold[which(volplot$significance == F)] <- "non significant"
  if(!is.null(int_effect_threshold)){
    volplot$threshold[abs(volplot$effect_sizes) >= abs(int_effect_threshold) & volplot$threshold == "significant"] <- "sig. interesting"
  }
  volplot$threshold <- factor(volplot$threshold, levels = c("sig. interesting", "significant", "non significant"))  
  
  #Plotting range thresholds
  x_lim <- max(abs(effect_sizes))
  y_lim <- c(0, range(-log10(volplot$pvals))[2]+1)

  #Colors
  cbColors <- c("#009E73", "#0072B2", "#D55E00")
  
  #Plotting with ggplot2
  drawplot <- ggplot(volplot, aes(x = effect_sizes, y = -log10(pvals), label = identifiers, color = threshold))
  drawplot <- drawplot + geom_point(alpha = 0.4, size = 1)
  if(!is.null(int_effect_threshold)){
    drawplot <- drawplot + geom_vline(xintercept = c(-int_effect_threshold, int_effect_threshold), linetype = "longdash")
  } 
  if(!is.null(top_names)){
    drawplot <- drawplot + geom_label_repel(data=head(volplot[order(volplot$pvals, decreasing= F),], n=top_names), show.legend = F)
  }
  drawplot <- drawplot + theme_bw()
  drawplot <- drawplot + xlim(c(-x_lim, x_lim))
  if(!is.null(y_lab)){
    drawplot <- drawplot + ylab(y_lab)
  } else{
    drawplot <- drawplot + ylab("-log10(P)")
  }
  if(!is.null(x_lab)){
    drawplot <- drawplot + xlab(x_lab)
  } else{
    drawplot <- drawplot + xlab("Mean effect size")
  }
  if(!is.null(title)){
    drawplot <- drawplot + ggtitle(title)
  } 
  
  drawplot <- drawplot + scale_color_manual(values = c(cbColors),
                                            drop = F)
  drawplot <- drawplot + guides(colour = guide_legend(override.aes = list(size = 10)))
  drawplot <- drawplot + theme(axis.text = element_text(size = 17), 
                               axis.title = element_text(size = 17, face = "bold"),
                               plot.title = element_text(size = 17, face = "bold"),
                               legend.title = element_blank(),
                               legend.text = element_text(size = 17),
                               legend.position = "bottom")
  drawplot
}