#' dmr_genome_plot
#'
#' Produces a dotplot of a specific (region of) CpG(s) given a matrix of beta-values superposed onto the RefSeq genomic location. 
#' @param name Name of the DMR/DMP
#' @param chr Chromosome (example "chr6")
#' @param betas Matrix of Beta with the CpG numbers as rownames
#' @param annotation.gr GenomicRanges object containing the annotation/features of the GMset (obtained through minfi::getAnnotation())
#' @param factor_interest A vector based on which you want to stratify your plot
#' @param start_dmr The start of the region you want to plot
#' @param end_dmr The end of the region you want to plot. If the region is only 1 bp, simply enter the same value as the start_dmr
#' @param flanks How large should the flanking region be? Defaults to the same size as the DMR.
#' @param SNP Do you want SNP information sourced from UCSC? Defaults to FALSE, because calling UCSC is slow.
#' @param Reg Do you want Regulatory information sourced from UCSC? Defaults to FALSE, because calling UCSC is slow.
#' @param dmps Do you want to indicate them individually? Defaults to FALSE
#' @param genome_version Which genome build to use? Defaults to "hg19"
#' @param diff_symbol Do you want the different groups to be represented by different symbols (TRUE)? Defaults to FALSE
#' @param dotsize The pixel size for plotting symbols. Defaults to 0.7 (standard of Gviz)
#' @param plotgrid Do you want to plot a grid (TRUE)? Defaults to TRUE
#' @param confint Do you want to plot the 95\% confidence intervals (TRUE)? Defaults to TRUE
#' @param highlight Do you want to highlight the region of interest (TRUE)? Defaults to TRUE
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A dotplot of a specific (region of) CpG(s) of beta-values superposed onto the RefSeq genome. Additionally, it returns a genomic ranges object with the CpGs within the chosen region.
#' @keywords cpg, dmr, methylation
#' @export
#' @import Gviz
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
#' annotation.gr <- makeGRangesFromDataFrame(getAnnotation(GMset), keep.extra.columns = T, start.field = "pos", end.field = "pos")
#' 
#' #Find DMRs
#' design <- model.matrix(~pData(Mset)$Sample_Group)
#' bumps <- bumphunter(object = GMset, design = design, coef = 2, cutoff = 0.2, B = 0, type = "Beta")
#' dmr <- bumps$table[1,]
#'
#' #Plot DMP
#' dmr_genomeplot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$start)
#'
#' #Plot DMR (without SNPs)
#' dmr_genomeplot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$end)
#' 
#' #Plot DMR (with SNPs)
#' dmr_genome_plot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$end, SNP = T) 

dmr_genome_plot <- function(name, chr, betas, annotation.gr, factor_interest, start_dmr, end_dmr, flanks, SNP = F, Reg = F, dmps = F, genome_version = "hg19", diff_symbol = T, dotsize = 0.7, plotgrid = T, confint = T, highlight = T){
  
  #Input check
  if(is.null(chr)) stop("No chr provided")
  if(!is.character(chr)) stop("chr is not a character, please input it as \"chr[1-9]\"")
  if(is.null(betas)) stop("No Beta matrix provided")
  if(is.null(annotation.gr)) stop("No annotation.gr provided")
  else if(!class(annotation.gr)[1] == "GRanges") stop("annotation.gr is not in GenomicRanges format")
  if(is.null(factor_interest)) stop("No factor of interest provided")
  if(ncol(betas) != length(factor_interest)) stop("Length of factor of interest does not equal number of samples")
  if(is.null(start_dmr)) stop("No start given")
  if(is.null(end_dmr)) stop("No end given")
  if(missing(flanks)){ #If no flanks are given, automate the flanks to be the same size as the DMR. If a DMP is provided, set it to 100 (arbitrarily chosen)
    region_size <- end_dmr-start_dmr
    if(region_size == 0){
      flanks <- 100
    } else{
      flanks <- region_size
    }
  }
  
  #Check DMR
  dmr.hg19.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_dmr & start(annotation.gr) <= end_dmr),]
  #print(dmr.hg19.annotation)
  
  #Add flanks
  start_dmr_flanks <- start_dmr-flanks
  end_dmr_flanks <- end_dmr+flanks
  
#   #Betas
#   GMset.beta <- as.data.frame(betas)
#   dmr.flanks.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_dmr_flanks & start(annotation.gr) <= end_dmr_flanks),]
#   dmr.meth <- data.frame(seqnames = as.character(seqnames(dmr.flanks.annotation)), pos = start(dmr.flanks.annotation), GMset.beta[names(dmr.flanks.annotation),])
#   colnames(dmr.meth) <- c("seqnames", "pos", colnames(GMset.beta))
#   dmr.meth.gr <- makeGRangesFromDataFrame(df = dmr.meth, keep.extra.columns = T, start.field = "pos", end.field = "pos")

  dmr.meth.gr <- methdata_to_granges(methdata = betas, annotation.gr = annotation.gr, chr = chr, start_region = start_dmr_flanks, end_region = end_dmr_flanks)$region.meth.gr
  
  #############
  # Plot:Gviz #
  #############
  
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
  grtrack <- UcscTrack(track = "RefSeq Genes", 
                       table = "refGene", 
                       trackType = "GeneRegionTrack", 
                       chromosome = chr, 
                       genome = genome_version, 
                       rstart = "exonStarts", 
                       rends = "exonEnds", 
                       symbol = "name2", 
                       strand = "strand", 
                       shape = "arrow", 
                       transcriptAnnotation = "symbol", 
                       name = "RefSeq", 
                       showId = T, 
                       from = start_dmr_flanks, 
                       to = end_dmr_flanks, 
                       fill = "#8282d2",
                       size = 0.4,
                       rot.title = 0,
                       fontsize = 10)
  
  factor_interest <- as.factor(factor_interest)
  
  #Grid
  if(plotgrid == T){
    type_plot <- c("a", "p", "g")
  } else{
    type_plot <- c("a", "p")
  }

  #Confidence intervals (95%)
  if(confint == T){
    type_plot <- c(type_plot, "confint")
  }
    
  #Change symbols if necessary
  if(diff_symbol == T){
    pchp <- as.numeric(as.factor(levels(factor_interest)))
  } else{
    pchp <- 20
  }
  #Change size symbols
  cexp <- dotsize
  
  dtrack.meth <- DataTrack(range = dmr.meth.gr, 
                           name = "Methylation", 
                           ylim = c(0,1), 
                           groups = factor_interest,
                           lty = sort(as.numeric(as.factor(levels(factor_interest)))),
                           pch = pchp,
                           cex = cexp,
                           type = type_plot, 
                           legend = TRUE,
                           fontsize = 10)
 
  ###################
  # Base track list #
  ###################
  tracklist <- c(itrack, gtrack, grtrack, dtrack.meth)
  
  ########################
  # Accessory track list #
  ########################
  
  #DMP track: Useful indication where the individual CpGs lie within a DMR when the number of CpGs within a DMR become large
  if(dmps){
    dmp.hg19.annotation <- dmr.hg19.annotation
    
    #Add the G and direction of each CpG
    start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"]) <- start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"])-1
    end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"]) <- end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"])+1
    dmptrack <- AnnotationTrack(dmp.hg19.annotation, name = "DMP", col = dmp.hg19.annotation$red_end)
  }
  
  #SNP track
  if(SNP){
    snptrack <- UcscTrack(genome = genome_version, 
                          chromosome = chr, 
                          track = "Common SNPs(141)", 
                          from = start_dmr_flanks, 
                          to = end_dmr_flanks, 
                          trackType = "AnnotationTrack", 
                          start = "chromStart", 
                          end = "chromEnd", 
                          gene = "name", 
                          feature = "func", 
                          strand = "strand", 
                          shape = "box",
                          stacking = "dense",
                          fill = "black", 
                          name = "SNPs",
                          fontsize = 10)
    
    tracklist <- c(tracklist, snptrack)
  }
  #Regulatory tracks
  #CpG Islands
  if(Reg){
    cgitrack <- UcscTrack(genome = genome_version, 
                          chromosome = chr,
                          track = "cpgIslandExt", 
                          from = start_dmr_flanks-1000, 
                          to = end_dmr_flanks+1000,
                          trackType = "AnnotationTrack", 
                          start = "chromStart",
                          end = "chromEnd", 
                          id = "name", 
                          shape = "box",
                          fill = "#006400",
                          name = "CGI",
                          size = 0.4,
                          rot.title = 0,
                          fontsize = 10)
    
#     chromhmm <- UcscTrack(genome = genome_version, 
#                           chromosome = chr,
#                           track = "Broad ChromHMM", 
#                           table = "wgEncodeBroadHmmGm12878HMM",
#                           from = start_dmr_flanks-1000, 
#                           to = end_dmr_flanks+1000,
#                           trackType = "AnnotationTrack", 
#                           start = "chromStart",
#                           end = "chromEnd", 
#                           feature = "name", 
#                           shape = "box",
#                           fill = "#960000",
#                           name = "ChromHMM",
#                           featureAnnotation = "feature",
#                           size = 0.4,
#                           rot.title = 0,
#                           fontsize = 8)
    
    #Not very relevant: No blood cells
    #Transcription Factor binding sites
    tfbstrack <- UcscTrack(genome = genome_version, 
                           chromosome = chr,
                           track = "Txn Factor ChIP", 
                           table = "wgEncodeRegTfbsClusteredV3", 
                           from = start_dmr_flanks-1000, 
                           to = end_dmr_flanks+1000,
                           trackType = "AnnotationTrack", 
                           start = "chromStart",
                           end = "chromEnd", 
                           feature = "name", 
                           shape = "box",
                           fill = "#960000",
                           name = "TFBS",
                           featureAnnotation = "feature",
                           size = 0.4,
                           rot.title = 0,
                           fontsize = 10)
    
    #Does not always work
#     #H3K27Ac
#     h3k27actrack <- UcscTrack(genome = genome_version, 
#                               chromosome = chr,
#                               track = "Broad Histone", 
#                               table = "wgEncodeBroadHistoneGm12878H3k27acStdSig",
#                               from = start_dmr_flanks-1000, 
#                               to = end_dmr_flanks+1000,
#                               trackType = "DataTrack", 
#                               start = "start", 
#                               end = "end", 
#                               data = "score", 
#                               type = "hist",
#                               #window = -1,
#                               #windowSize = 50,
#                               fill.histogram = "black",
#                               name = "H3K27Ac", 
#                               ylim = c(0,100), 
#                               size = 1.3,
#                               fontsize = 8)
    #H3K4Me1
#     h3k4me1track <- UcscTrack(genome = genome_version, 
#                               chromosome = chr,
#                               track = "Broad Histone", 
#                               table = "wgEncodeBroadHistoneGm12878H3k04me1StdSigV2",
#                               from = start_dmr_flanks-1000, 
#                               to = end_dmr_flanks+1000,
#                               trackType = "DataTrack", 
#                               start = "start", 
#                               end = "end", 
#                               data = "score", 
#                               type = "hist",
#                               window = -1,
#                               windowSize = 50,
#                               fill.histogram = "black",
#                               name = "H3K4Me1", 
#                               ylim = c(0,100),
#                               size = 1.3,
#                               fontsize = 8)
#     #H3K4Me3
#     h3k4me3track <- UcscTrack(genome = genome_version, 
#                               chromosome = chr,
#                               track = "Broad Histone", 
#                               table = "wgEncodeBroadHistoneGm12878H3k04me3StdSigV2",
#                               from = start_dmr_flanks-1000, 
#                               to = end_dmr_flanks+1000,
#                               trackType = "DataTrack", 
#                               start = "start", 
#                               end = "end", 
#                               data = "score", 
#                               type = "hist",
#                               window = -1,
#                               windowSize = 10,
#                               fill.histogram = "black",
#                               name = "H3K4Me3", 
#                               ylim = c(0,100),
#                               size = 1.3,
#                               fontsize = 8)
#     #DNaseI
#     dnasetrack <- UcscTrack(genome = genome_version, 
#                             chromosome = chr,
#                             track = "DNase Clusters", 
#                             table = "wgEncodeRegDnaseClusteredV3",
#                             from = start_dmr_flanks-1000, 
#                             to = end_dmr_flanks+1000,
#                             trackType = "DataTrack", 
#                             start = "start", 
#                             end = "end", 
#                             data = "score", 
#                             type = "hist",
#                             window = -1,
#                             windowSize = 10,
#                             fill.histogram = "black",
#                             name = "DNaseI",
#                             ylim = c(0,1000),
#                             size = 1.3,
#                             fontsize = 8)
    
    #tracklist <- c(tracklist, cgitrack, h3k4me1track, h3k4me3track, dnasetrack, tfbstrack)
    tracklist <- c(tracklist, cgitrack, tfbstrack)
  }
  
  #Highlight track: Useful indication of where the region of interest was found
  if(end_dmr_flanks-start_dmr_flanks < 100){
    offset <- 1
  } else{
    offset <- 5
  }
  
  if(highlight == T){
    alltrack <- HighlightTrack(trackList = tracklist,
                               start = start_dmr-offset,
                               end = end_dmr+offset,
                               chromosome = chr)  
  } else{
    alltrack <- tracklist
  }
  
  plotTracks(trackList = alltrack, 
             from = start_dmr_flanks, 
             to = end_dmr_flanks, 
             main = name,
             cex.title = 1.5,
             cex.axis = 1.5,
             cex.legend = 1.5, 
             add53 = T,
             add35 = T,
             littleTicks = F,
             #fontsize = 12, 
             fontcolor = "black",
             background.title = "darkgray") 
  return(dmr_annotation = dmr.hg19.annotation)
}