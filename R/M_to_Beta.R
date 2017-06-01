#' M_to_Beta
#'
#' Beta values represent the level of methylation on a scale of 0 and 100, where 0 means 0% methylation and 100 means 100% methylation. Despite the interpretability of Beta values, Du et al. 2010 argue that Beta values are more heteroscedastic than M-values, making them less reliable for downstream statistical analysis. This function converts mean M values into their respective Beta values for methylation analysis. 
#' @param mvals A numeric vector containing the M-values
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return Beta values 
#' @keywords Beta, M-values
#' @export
#' @examples 
#' Mvals <- c(1.5184121,3.8591045,1.7411637,2.3336920,1.3309937,-4.1895298,-4.0174476,-0.9808680,1.4356915,-0.5294438)
#' Betas <- MtoBeta(Mvals)

MtoBeta <- function(mvals){
  if(is.null(mvals)) stop("No M values defined")
  else if(is.character(mvals)) stop("Provided M values are non-numeric")
  mvals <- as.numeric(mvals)
  
  (2^(mvals))/(2^(mvals)+1)
}