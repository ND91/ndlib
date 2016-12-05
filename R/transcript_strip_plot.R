#' transcript_strip_plot 
#'
#' Produces a a stripplot of a specific transcript given a count matrix. Ensure that the count matrix has the identifier as the row names
#' @param id Identifier as provided in the count matrix
#' @param counts Matrix of counts with the identifiers as rownames
#' @param factor_interest A vector based on which you want to stratify your plot
#' @param gg.plot Do you want to plot using ggplot2 (TRUE) or base R plotting (FALSE). Defaults to TRUE.
#' @param enlarged Do you want an enlarged plot on the right side (TRUE or FALSE). Defaults to "auto" whereby the standard deviation will be calculated of the sample. If the sd is smaller than 0.015 the enlarged plot will be generated for visibility purposes. 
#' @param type Choose the type of plot you want: "boxplot", "SE", "CI" or simply "mean". Defaults to "boxplot".
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A stripplot of a specific gene
#' @keywords stripplot, boxplot, expression
#' @export
#' @importFrom Hmisc smean.cl.normal 
#' @importFrom Hmisc smean.cl.boot 
#' @importFrom gridExtra grid.arrange 
#' @importFrom grid textGrob
#' @import ggplot2
#' @examples 
#' Beta <- matrix(c(0.15,0.2,0.17,0.76,0.8,0.65,0.22,0.23,0.24,0.5,0.51,0.52), nrow = 2,ncol = 6, byrow = T)
#' rownames(Beta) <- c("cg123456", "cg789010")
#' Pheno <- c(rep("pheno1", 3), rep("pheno2", 3))
#' cpg_strip_plot("cg123456", Beta, Pheno, type = "SE")

transcript_strip_plot <- function(id, counts, factor_interest, title, gg.plot = T, enlarged = "auto", type = "boxplot"){

  #Convert the counts to a matrix, else the dataframe calling function won't work properly
  counts <- as.matrix(counts)
  factor_interest <- as.factor(factor_interest)
  
  if(is.null(id)) stop("No ID defined")
  if(is.null(counts)) stop("Count matrix provided")
  if(is.null(rownames(counts))) stop("Count matrix does not have any IDs")
  if(!(id %in% rownames(counts))) stop("ID does not exist in given count matrix")
  if(ncol(counts) != length(factor_interest)) stop("Length of factor of interest does not equal number of samples") 
  if(is.null(factor_interest)) stop("No factor of interest provided") 
  if(missing(title)) title <- id
  
  counts.df <- data.frame(count = counts[id,], Group = factor_interest)
  #Not sure of declaring enlarged as a boolean with an extra value (Troolean?) is a good thing. Similarly, not sure if overwriting the initial value afterwards to a Boolean is a good thing either. 
  if(gg.plot == T){
    #Overall plot
    overall_plot <- ggplot(data = counts.df, aes(x = Group, y = count)) + 
      theme_bw() + 
      xlab("") + 
      scale_shape_manual(values=(1:nlevels(counts.df$Group))%%10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 17), 
            axis.title = element_text(size = 17, face = "bold"),
            legend.title = element_text(size = 17, face = "bold"),
            legend.text = element_text(size = 17),
            plot.title = element_text(face = "bold"))
    
    if(type == "boxplot"){
      overall_plot <- overall_plot + 
        geom_boxplot(aes(fill = Group), outlier.colour = NA)
    } else if(type == "SE"){ #Note that stat_summary calls functions from the Hmisc package
      #The following sections call the stat_summary function, which can call functions from the Hmisc package. I highly recommend you to read the following threads to understand what is happening:
      # https://stackoverflow.com/questions/19258460/standard-error-bars-using-stat-summary
      # http://r.789695.n4.nabble.com/ggplot-stat-summary-mean-cl-boot-td4021218.html
      overall_plot <- overall_plot + 
        stat_summary(fun.data = mean_se, geom = "crossbar",  aes(color = Group))
      
    } else if(type == "CI"){ 
      para <- T
      
      if(para == T){
        #mean_cl_normal: Calls Hmisc::smean.cl.normal and is meant for normal distributions
        stat.function <- stat_summary(fun.data = mean_cl_normal, geom = "crossbar",  aes(color = Group))
      } else{
        #mean_cl_boot: Calls Hmisc::smean.cl.boot and is meant for non-parametric distributions, CIs are calculated through bootstraps, which we set to 5000
        B <- 5000
        conf.int <- 0.95
        stat.function <- stat_summary(fun.data = mean_cl_boot, fun.args = list(B = B, conf.int = conf.int), geom = "crossbar",  aes(color = Group))
      }
      overall_plot <- overall_plot + 
        stat.function
    } else if(type == "mean"){
      overall_plot <- overall_plot + 
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "crossbar",  aes(color = Group))
    } else stop("Incorrect value entered for the \"type\" argument! Can only be \"boxplot\", \"SE\", \"CI\" or \"mean\"")
    
    if(enlarged == "auto"){
      if(sd(counts[id,]) < 0.015){ #A standard deviation of 0.015 is chosen somewhat arbitrarily
        message("Standard deviation of sample is smaller than 0.015. Enlarged plot will be generated")
        enlarged <- T
      } else{
        message("Standard deviation of sample is larger than 0.015. Enlarged plot is not necessary")
        enlarged <- F
      }
    } 
    if(enlarged == T){
      large_plot <- overall_plot +
        geom_jitter(size = 2, aes(shape = Group)) +
        ggtitle("Enlarged") + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 17), 
              axis.title = element_blank())
      
      overall_plot <- overall_plot + 
        ylab("Counts")
      
      grid.arrange(overall_plot, large_plot, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 17, fontface = "bold")))
    } else if(enlarged == F){
      overall_plot <- overall_plot + 
        geom_jitter(size = 2, aes(shape = Group)) + 
        ggtitle(title)
      overall_plot
    } else stop("Incorrect value entered for the \"enlarged\" argument!")
  } else{
    set.seed(1)
    boxplot(beta~Group, data = counts.df, main = id, ylab = "Beta", main = title)
    stripchart(beta~Group, data = counts.df, vertical = T, method = "jitter", pch = 16, add = T, col = "blue")
  }
}