#' minfi_to_GEO
#'
#' Provide a GEO submission for Illumina methylation data according to the template provided with GPL16304
#' @param values A matrix or dataframe containing the methylation values (i.e. Beta or M-values)
#' @param detP A matrix containing the detection P values per probe
#' @param unmeth A matrix containing the unmethylated values per probe with the appropriate row names
#' @param meth A matrix containing the methylated values per probe with the appropriate column names
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A list of two dataframes: 1) The unmethylated, methylated and detection P-value per sample per probe, and 2) The degree of methylation (Beta, M-value) and the detection P-value per sample per probe.
#' @keywords cpg, methylation, 450k, EPIC, GEO
#' @export
#' @examples 
#' 
#' #"I generally provide GEO the processed beta/M-values, and the unprocessed methylated and unmethylated channels as unprocessed data. The detection P-values should be the same regardless of processing"
#' 
#' #Load data
#' require(minfiData)
#' baseDir <- system.file("extdata", package = "minfiData")
#' targets <- read.metharray.sheet(base = baseDir)
#' RGset <- read.metharray.exp(targets = targets, recursive = T)
#' 
#' #Preprocess it raw to get the raw methylated and unmethylated channels
#' Mset.raw <- preprocessRaw(RGset)
#' meth.raw <- getMeth(Mset.raw)
#' unmeth.raw <- getUnmeth(Mset.raw)
#' detP.raw <- detectionP(RGset)
#' 
#' #Preprocess it using your normalization tool of interest
#' Mset <- preprocessIllumina(rgSet = RGset, bg.correct = T, normalize = "controls", reference = 2)
#' Rset <- ratioConvert(Mset)
#' GMset <- mapToGenome(Rset)
#' 
#' #Beta values
#' Beta.processed <- getBeta(GMset)
#' 
#' geo.submission <- minfi_to_GEO(Beta.processed, detP.raw, unmeth.raw, meth.raw)
#' 

minfi_to_GEO <- function(values, detP, unmeth, meth){
  #Sanitization
  if(is.null(values)) stop("No Beta or M-values provided")
  if(is.null(detP)) stop("No detection P-values provided")
  else if(min(detP) < 0 | max(detP) > 1) stop("detP is smaller than 0, or larger than 1, which is not possible for a P-value")
  if(ncol(detP) != ncol(values)) stop("The matrices containing the detection P value and the Beta/M-values contain unequal samples")
  if(nrow(detP) != nrow(values)) stop("The matrices containing the detection P value and the Beta/M-values contain unequal CpGs")
  if(is.null(unmeth)) stop("No unmethylated values provided")
  if(is.null(meth)) stop("No methylated values provided")
  if(ncol(unmeth) != ncol(meth)) stop("The matrices containing the methylated and unmethylated channel contain unequal samples")
  else if(nrow(unmeth) != nrow(meth)) stop("The matrices containing the methylated and unmethylated channel contain unequal CpGs")
  else if(ncol(unmeth) != ncol(detP)) stop("The matrices containing the methylated and unmethylated channel, and the detection P values contain unequal samples")
  else if(nrow(unmeth) != nrow(detP)) stop("The matrices containing the methylated and unmethylated channel, and the detection P values contain unequal CpGs")
  if("FALSE" %in% names(table(rownames(unmeth) == rownames(meth)))) stop("The row names of the methylated and unmethylated channel do not correspond")
  if("FALSE" %in% names(table(colnames(unmeth) == colnames(meth)))) stop("The column names of the methylated and unmethylated channel do not correspond")
  
  #Creation of the Beta/M-value, detP matrix
  rows.combined <- nrow(values)
  cols.combined <- ncol(values) + ncol(detP)
  values.detP <- matrix(NA, nrow = rows.combined, ncol = cols.combined)
  values.detP[, seq(1, cols.combined, 2)] <- values
  values.detP[, seq(2, cols.combined, 2)] <- detP
  
  colnames(values.detP) <- 1:ncol(values.detP)
  colnames(values.detP)[seq(1, cols.combined, 2)] <- colnames(values)
  colnames(values.detP)[seq(2, cols.combined, 2)] <- paste(colnames(values), "Detection Pval")
  rownames(values.detP) <- rownames(values)
  
  #Creation of the resulting unmethylated, methylated, detP matrix
  rows.combined <- nrow(meth)
  cols.combined <- ncol(meth) + ncol(unmeth) + ncol(detP)
  
  unmeth.meth.detP <- matrix(NA, nrow = rows.combined, ncol = cols.combined)
  unmeth.meth.detP[, seq(1, cols.combined, 3)] <- unmeth
  unmeth.meth.detP[, seq(2, cols.combined, 3)] <- meth
  unmeth.meth.detP[, seq(3, cols.combined, 3)] <- detP
  
  colnames(unmeth.meth.detP) <- 1:ncol(unmeth.meth.detP)
  colnames(unmeth.meth.detP)[seq(1, cols.combined, 3)] <- paste(colnames(meth), "Unmethylated Signal")
  colnames(unmeth.meth.detP)[seq(2, cols.combined, 3)] <- paste(colnames(meth), "Methylated Signal")
  colnames(unmeth.meth.detP)[seq(3, cols.combined, 3)] <- paste(colnames(meth), "Detection Pval")
  rownames(unmeth.meth.detP) <- rownames(meth)
  
  return(list(unmeth.meth.detP = unmeth.meth.detP, values.detP = values.detP))
}