#' cpg_strip_plot
#'
#' Produces a boxplot overlayed by a stripplot of a specific CpG given a matrix of beta-values. Ensure that the Beta-matrix has the CpG number as the row names
#' @param cpg_num CpG number (ID)
#' @param betas Matrix of Beta with the CpG numbers as rownames
#' @param factor_interest A vector based on which you want to stratify your plot
#' @param enlarged Do you want an enlarged plot on the right side (TRUE or FALSE). Defaults to "auto" whereby the standard deviation will be calculated of the sample. If the sd is smaller than 0.015 the enlarged plot will be generated for visibility purposes.
#' @param type Choose the type of plot you want: "boxplot", "SE", "CI" or simply "mean". Defaults to "boxplot".
#' @param colors A vector containing the the colors to be utilized for the plot with the contrast of interest as names (default: uses ggplot's default coloring scheme).
#' @param legend Boolean depicting whether a legend should be included (default: TRUE)
#' @param y_lab A specific name for the y labels (default: "Beta")
#'
#' @author Andrew Y.F. Li Yim
#'
#' @return A boxplot overlayed by a stripplot of a specific CpG
#' @keywords cpg, stripplot, boxplot, methylation
#' @export
#' @importFrom Hmisc smean.cl.normal
#' @importFrom Hmisc smean.cl.boot
#' @importFrom ggpubr ggarrange
#' @importFrom ggpubr annotate_figure
#' @import ggplot2
#' @examples
#' Beta <- matrix(c(0.15,0.2,0.17,0.76,0.8,0.65,0.22,0.23,0.24,0.5,0.51,0.52), nrow = 2,ncol = 6, byrow = T)
#' rownames(Beta) <- c("cg123456", "cg789010")
#' Pheno <- c(rep("pheno1", 3), rep("pheno2", 3))
#' cpg_strip_plot("cg123456", Beta, Pheno, type = "SE")

cpg_strip_plot <- function(cpg_num, betas, factor_interest, title, enlarged = "auto", type = c("boxplot", "CI", "SE", "mean"), colors = NULL, legend = TRUE, y_lab = NULL){
  #Convert the betas to a matrix, else the dataframe calling function won't work properly
  betas <- as.matrix(betas)
  factor_interest <- as.factor(factor_interest)
  type <- match.arg(type)

  if(is.null(cpg_num)) stop("No CpG number defined")
  if(is.null(betas)) stop("No Beta matrix provided")
  if(is.null(rownames(betas))) stop("Betas do not have any CpG rownames")
  if(!(cpg_num %in% rownames(betas))) stop("CpG number does not exist in given Beta matrix")
  if(ncol(betas) != length(factor_interest)) stop("Length of factor of interest does not equal number of samples")
  if(is.null(factor_interest)) stop("No factor of interest provided")
  if(missing(title)) title <- cpg_num
  if(is.null(y_lab)) y_lab <- "Beta"
  if(!is.null(colors)){
    if(!all(factor_interest %in% names(colors))) stop("Factor of interest cannot be found in the name of \"colors\", \"colors\" must be a named vector")
  }

  beta.df <- data.frame(beta = betas[cpg_num,], Cohort = factor_interest)
  #Not sure of declaring enlarged as a boolean with an extra value (Troolean?) is a good thing. Similarly, not sure if overwriting the initial value afterwards to a Boolean is a good thing either.
  #Overall plot
  overall_plot <- ggplot(data = beta.df, aes(x = Cohort, y = beta)) +
    theme_bw() +
    xlab("") +
    ylab(y_lab) +
    scale_shape_manual(values=(1:nlevels(beta.df$Cohort))%%10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          plot.title = element_text(face = "bold"))

  if(!is.null(colors)){
    overall_plot <- overall_plot + scale_fill_manual(values = colors)
  }

  if(legend == F){
    overall_plot <- overall_plot + theme(legend.position = "none")
  } else{
    overall_plot <- overall_plot + theme(axis.text.x = element_blank())
  }

  if(type == "boxplot"){
    overall_plot <- overall_plot +
      geom_boxplot(aes(fill = Cohort), outlier.colour = NA)
  } else if(type == "SE"){ #Note that stat_summary calls functions from the Hmisc package
    #The following sections call the stat_summary function, which can call functions from the Hmisc package. I highly recommend you to read the following threads to understand what is happening:
    # https://stackoverflow.com/questions/19258460/standard-error-bars-using-stat-summary
    # http://r.789695.n4.nabble.com/ggplot-stat-summary-mean-cl-boot-td4021218.html
    overall_plot <- overall_plot +
      stat_summary(fun.data = mean_se, geom = "crossbar",  aes(color = Cohort))

  } else if(type == "CI"){
    para <- T

    if(para == T){
      #mean_cl_normal: Calls Hmisc::smean.cl.normal and is meant for normal distributions
      stat.function <- stat_summary(fun.data = mean_cl_normal, geom = "crossbar",  aes(color = Cohort))
    } else{
      #mean_cl_boot: Calls Hmisc::smean.cl.boot and is meant for non-parametric distributions, CIs are calculated through bootstraps, which we set to 5000
      B <- 5000
      conf.int <- 0.95
      stat.function <- stat_summary(fun.data = mean_cl_boot, fun.args = list(B = B, conf.int = conf.int), geom = "crossbar",  aes(color = Cohort))
    }
    overall_plot <- overall_plot +
      stat.function
  } else if(type == "mean"){
    overall_plot <- overall_plot +
      stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "crossbar",  aes(color = Cohort))
  }

  if(enlarged == "auto"){
    if(sd(betas[cpg_num,]) < 0.015){ #A standard deviation of 0.015 is chosen somewhat arbitrarily
      message("Standard deviation of sample is smaller than 0.015. Enlarged plot will be generated")
      enlarged <- T
    } else{
      message("Standard deviation of sample is larger than 0.015. Enlarged plot is not necessary")
      enlarged <- F
    }
  }
  if(enlarged == T){
    large_plot <- overall_plot +
      geom_jitter(size = 2, aes(shape = Cohort)) +
      ggtitle("Enlarged") +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 12),
            axis.title = element_blank())

    overall_plot <- overall_plot +
      coord_cartesian(ylim=c(0,1)) +
      ylab("Beta")

    if(legend == F){
      large_plot <- large_plot + theme(legend.position = "none")
      comb_plot <- ggarrange(overall_plot, large_plot, ncol = 2)
    } else{
      large_plot <- large_plot + theme(axis.text.x = element_blank())
      comb_plot <- ggarrange(overall_plot, large_plot, ncol = 2, common.legend = T, legend = "bottom")
    }

    comb_plot <- annotate_figure(comb_plot, text_grob(title, size = 17, face = "bold"))
    return(comb_plot)
  } else if(enlarged == F){
    overall_plot <- overall_plot +
      geom_jitter(size = 2, aes(shape = Cohort)) +
      coord_cartesian(ylim=c(0,1)) +
      ggtitle(title)
    overall_plot
  } else stop("Incorrect value entered for the \"enlarged\" argument!")
}
