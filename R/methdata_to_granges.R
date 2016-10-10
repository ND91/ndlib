#' methdata_to_granges
#'
#' Produces a dotplot of a specific (region of) CpG(s) given a matrix of beta-values superposed onto the RefSeq genomic location. 
#' @param methdata Methylation data (betas or M-values)
#' @param annotation.gr Annotation data in GenomicRanges format
#' @param chr The chromosome as character (ex. "chr6")
#' @param start_region The start of the region as numeric
#' @param end_region The end of the region as numeric
#' @param flanks The size of the around the region of interest. Defaults to 0
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A region of the annotation.gr and betas given the coordinates in GenomicRanges format
#' @keywords methylation, dmr, cpg, granges
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
#' annotation.gr <- makeGRangesFromDataFrame(getAnnotation(GMset), keep.extra.columns = T, start.field = "pos", end.field = "pos")
#' 
#' methdata_to_granges(methdata = Beta, annotation.gr = annotation.gr, chr = "chr6", start_region = 70672841, end_region = 70672878)

methdata_to_granges <- function(methdata, annotation.gr, chr, start_region, end_region, flanks){
  #Input check
  if(is.null(methdata)) stop("No methdata (Betas or M-values) provided")
  if(is.null(annotation.gr)) stop("No annotation (in GenomicRanges format) provided")
  if(is.null(chr)){
    stop("No chr provided")
  } else if(!is.character(chr)) stop("chr is not a character, please input it as \"chr[1-9]\"")
  if(is.null(start_region)){
    stop("No start position provided")
  } else if(!is.numeric(start_region)) stop("start_region is not a numeric")
  if(is.null(end_region)){
    stop("No end position provided")
  } else if(!is.numeric(end_region)) stop("end_region is not a numeric")
  
  #In case the wrong direction has been provided
  if(start_region > end_region){
    x <- end_region
    end_region <- start_region
    start_region <- x
  }
  
  #Option to find CpGs in flanking regions
  if(missing(flanks)){
    flanks <- 0
  }
  start_region <- start_region + flanks
  end_region <- end_region + flanks
  
  GMset.beta <- as.data.frame(methdata)
  
  #Find all the hits in the annotation.gr
  region.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_region & start(annotation.gr) <= end_region),]
  if(length(region.annotation)==0) stop("Cannot find any CpGs within the given coordinates!")
  
  region.meth <- data.frame(seqnames = as.character(seqnames(region.annotation)), pos = start(region.annotation), GMset.beta[names(region.annotation),])
  colnames(region.meth) <- c("seqnames", "pos", colnames(GMset.beta))
  region.meth.gr <- makeGRangesFromDataFrame(df = region.meth, keep.extra.columns = T, start.field = "pos", end.field = "pos")
  
  return(list(region.annotation.gr = region.annotation, region.meth.gr = region.meth.gr))
}

