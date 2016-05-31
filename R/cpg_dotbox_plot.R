#' cpg_dotbox_plot 
#'
#' Produces a boxplot overlayed by a stripplot of a specific CpG given a matrix of beta-values. Ensure that the Beta-matrix has the CpG number as the row names
#' @param cpg_num CpG number (ID)
#' @param betas Matrix of Beta with the CpG numbers as rownames
#' @param factor_interest A vector based on which you want to stratify your plot
#' @param gg.plot Do you want to plot using ggplot2 (TRUE) or base R plotting (FALSE). Defaults to TRUE.
#' @param enlarged Do you want an enlarged plot on the right side (TRUE or FALSE). Defaults to "auto" whereby the standard deviation will be calculated of the sample. If the sd is smaller than 0.015 the enlarged plot will be generated for visibility purposes. 
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A boxplot overlayed by a stripplot of a specific CpG
#' @keywords cpg, jitter, boxplot, methylation
#' @export
#' @examples 
#' Beta <- matrix(c(0.15,0.2,0.17,0.76,0.8,0.65,0.22,0.23,0.24,0.5,0.51,0.52), nrow = 2,ncol = 6, byrow = T)
#' rownames(Beta) <- c("cg123456", "cg789010")
#' Pheno <- c(rep("pheno1", 3), rep("pheno2", 3))
#' cpg_dotbox_plot("cg123456", Beta, Pheno)


#TODO: Give it confidence intervals options. 
cpg_dotbox_plot <- function(cpg_num, betas, factor_interest, title, gg.plot = T, enlarged = "auto"){

  #Convert the betas to a matrix, else the dataframe calling function won't work properly
  betas <- as.matrix(betas)
  factor_interest <- as.factor(factor_interest)
  
  if(is.null(cpg_num)) stop("No CpG number defined")
  if(is.null(betas)) stop("No Beta matrix provided")
  if(is.null(rownames(betas))) stop("Betas do not have any CpG rownames")
  if(!(cpg_num %in% rownames(betas))) stop("CpG number does not exist in given Beta matrix")
  if(ncol(betas) != length(factor_interest)) stop("Length of factor of interest does not equal number of samples") 
  if(is.null(factor_interest)) stop("No factor of interest provided") 
  if(missing(title)) title <- cpg_num
  
  #Not sure of declaring enlarged as a boolean with an extra value (Troolean?) is a good thing. Additionally, not sure if overwriting the initial value afterwards to a Boolean is a good thing either. 
  beta.df <- data.frame(beta = betas[cpg_num,], Cohort = factor_interest)
  if(gg.plot == T){
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
      set.seed(1)
      small_plot <- ggplot(data = beta.df, aes(x = Cohort, y = beta)) + 
        ylim(c(0,1)) + 
        geom_boxplot(aes(fill = Cohort), outlier.colour = NA) + 
        #geom_jitter(size = 1.5, aes(shape = Cohort)) +
        theme_bw() + 
        xlab("") + 
        ylab("Beta") +
        scale_shape_manual(values=(1:nlevels(beta.df$Cohort))%%10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 17), 
              axis.title = element_text(size = 17, face = "bold"),
              legend.title = element_text(size = 17, face = "bold"),
              legend.text = element_text(size = 17))
      
      set.seed(1)
      large_plot <- ggplot(data = beta.df, aes(x = Cohort, y = beta)) + 
        geom_boxplot(aes(fill = Cohort), outlier.colour = NA) + 
        geom_jitter(size = 2.5, aes(shape = Cohort)) + 
        theme_bw() + 
        ggtitle("Enlarged") + 
        theme(legend.position = "none") + 
        xlab("") + 
        scale_shape_manual(values=(1:nlevels(beta.df$Cohort))%%10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 17), 
              axis.title = element_blank())
      grid.arrange(small_plot, large_plot, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 17, fontface = "bold")))
    } else if(enlarged == F){
      set.seed(1)
      small_plot <- ggplot(data = beta.df, aes(x = Cohort, y = beta)) + 
        ylim(c(0,1)) + 
        geom_boxplot(aes(fill = Cohort), outlier.colour = NA) + 
        geom_jitter(size = 2, aes(shape = Cohort)) +
        theme_bw() + 
        xlab("") + 
        scale_shape_manual(values=(1:nlevels(beta.df$Cohort))%%10) +
        ggtitle(title) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 17), 
              axis.title = element_text(size = 17, face = "bold"),
              legend.title = element_text(size = 17, face = "bold"),
              legend.text = element_text(size = 17),
              plot.title = element_text(face = "bold"))
      small_plot
    } else stop("Incorrect value entered for \"enlarged\"!")
  } else{
    set.seed(1)
    boxplot(beta~Cohort, data = beta.df, main = cpg_num, ylim = c(0,1), ylab = "Beta", main = title)
    stripchart(beta~Cohort, data = beta.df, vertical = T, method = "jitter", pch = 16, add = T, col = "blue")
  }
}