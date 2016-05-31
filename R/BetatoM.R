#' BetatoM
#'
#' Beta values represent the level of methylation on a scale of 0 and 100, where 0 means 0% methylation and 100 means 100% methylation. Despite the interpretability of Beta values, Du et al. 2010 argue that Beta values are more heteroscedastic than M-values, making them less reliable for downstream statistical analysis. This function converts Beta effect sizes into mean M effect sizes for methylation analysis. 
#' @param betas A numeric vector containing the Beta values
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return M values
#' @keywords Beta, M-values
#' @export
#' @examples 
#' Beta <- c(0.74125143,0.93553086,0.76974692,0.83446278,0.71556639,0.05195812,0.05815754,0.33628677,0.73010305,0.40927036)
#' Mvals <- BetatoM(Beta)

BetatoM <- function(betas){
  if(is.null(betas)) stop("No Beta values defined")
  else if(is.character(betas)) stop("Provided Betas are non-numeric")
  betas <- as.numeric(betas)
  
  log2(betas/(1-betas))
}