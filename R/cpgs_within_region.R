#' cpgs_within_region
#'
#' Produces a dotplot of a specific (region of) CpG(s) given a matrix of beta-values superposed onto the RefSeq genomic location. 
#' @param chr Chromosome (example "chr6")
#' @param annotation.gr GenomicRanges object containing the annotation/features of the GMset (obtained through minfi::getAnnotation())
#' @param start_region The start of the region you want to plot
#' @param end_region The end of the region you want to plot. If the region is only 1 bp, simply enter the same value as the start_dmr
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A region of CpG(s) within the range of interest
#' @keywords cpg, dmr, methylation
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
#' Beta <- getBeta(GMset)
#' annotation.gr <- makeGRangesFromDataFrame(getAnnotation(GMset), keep.extra.columns = T, start.field = "pos", end.field = "pos")
#' 
#' #Find DMRs
#' design <- model.matrix(~pData(Mset)$Sample_Group)
#' bumps <- bumphunter(object = GMset, design = design, coef = 2, cutoff = 0.2, B = 0, type = "Beta")
#' dmr <- bumps$table[1,]
#'
#' #Find DMPs within DMRs
#' cpg_within_region(chr = "chr6", annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$start)
#'

cpgs_within_region <- function(chr, annotation.gr, start_region, end_region){
  
  #Input check
  if(is.null(chr)) stop("No chr provided")
  if(!is.character(chr)) stop("chr is not a character, please input it as \"chr[1-9]\"")
  if(is.null(start_region)) stop("No start given")
  if(is.null(end_region)) stop("No end given")
  
  #Check Region
  dmr.hg19.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_region & start(annotation.gr) <= end_region),]
  
  return(cpgs = dmr.hg19.annotation)
}