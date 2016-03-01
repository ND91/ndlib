#' gene_cpg_finder
#'
#' Find the probes associated to a given gene symbol (UCSC) using Regex
#' @param gene_symbol A (list of) UCSC gene symbol(s)
#' @param annotation.gr GenomicRanges object containing the annotation/features of the GMset (obtained through minfi::getAnnotation())
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A (list of) GenomicRanges object(s) containing the CpG annotations associated to the gene of interest
#' @keywords cpg, methylation
#' @export
#' @examples 
#' #Load data
#' require(minfiData)
#' baseDir <- system.file("extdata", package = "minfiData")
#' targets <- read.450k.sheet(base = baseDir)
#' RGset <- read.450k.exp(targets = targets, recursive = T)
#' Mset <- preprocessIllumina(rgSet = RGset, bg.correct = T, normalize = "controls", reference = 2)
#' Rset <- ratioConvert(Mset)
#' GMset <- mapToGenome(Rset)
#' 
#' annotation.gr <- makeGRangesFromDataFrame(getAnnotation(GMset), keep.extra.columns = T, start.field = "pos", end.field = "pos")
#' 
#' gene_cpg_finder(gene_symbol = "TNF", annotation.gr = annotation.gr)

#Extract the CpGs that are associated to the gene_symbol variable using a grep function. 
#The "(^|;)*" and "($|;)*" are used to find the gene specifically generically named gene symbols. For instance, simply grepping with TNF will also yield TNFSF4, TNFa, TNFb etc.

gene_cpg_finder <- function(gene_symbol, annotation.gr){
  if(is.null(gene_symbol)) stop("No UCSC gene symbol provided")
  if(is.null(annotation.gr)) stop("No annotation.gr provided")
  else if(!class(annotation.gr)[1] == "GRanges") stop("annotation.gr is not in GenomicRanges format")
  
  #If two or more genes are provided, yield a list of GenomicRanges objects
  if(length(gene_symbol) == 1){
    cpgs <- annotation.gr[grep(paste0("(^|;)*", gene_symbol, "($|;)*"), annotation.gr$UCSC_RefGene_Name), ]
  } else {
    cpgs <- sapply(gene_symbol, function(x) annotation.gr[grep(paste0("(^|;)*", x, "($|;)*"), annotation.gr$UCSC_RefGene_Name), ])
  }
  
  return(cpg = cpgs)
}