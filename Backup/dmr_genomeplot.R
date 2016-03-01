#' dmr_genomeplot
#'
#' Produces a dotplot of a specific (region of) CpG(s) given a matrix of beta-values superposed onto the RefSeq genomic location. 
#' @param name Name of the DMR/DMP
#' @param chr Chromosome (example "chr6")
#' @param betas Matrix of Beta with the CpG numbers as rownames
#' @param factor_interest A vector based on which you want to stratify your plot
#' @param start_dmr The start of the region you want to plot
#' @param end_dmr The end of the region you want to plot. If the region is only 1 bp, simply enter the same value as the start_dmr
#' @param flanks How large should the flanking region be? Defaults to the same size as the DMR.
#' @param simplified Do you simply want the position/region overlayed on the refseq genome, or do you want SNP information as well? Defaults to TRUE because it is faster.
#' @param dmps Do you want to indicate them individually? Defaults to FALSE
#' @param genome_version Which genome build to use? Defaults to "hg19"
#' @param diff_symbol Do you want the different groups to be represented by different symbols (TRUE)? Defaults to FALSE
#' @param cex The pixel size for plotting symbols. Defaults to 0.7 (standard of Gviz)
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A dotplot of a specific (region of) CpG(s) of beta-values superposed onto the RefSeq genome. Additionally, it returns a genomic ranges object with the CpGs within the chosen region.
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
#' #Plot DMP
#' dmr_genomeplot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$start)
#'
#' #Plot DMR (without SNPs)
#' dmr_genomeplot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$end)
#' 
#' #Plot DMR (with SNPs)
#' dmr_genomeplot(name = "DMR_1", chr = "chr6", betas = Beta, factor_interest = pData(Mset)$Sample_Group, annotation.gr = annotation.gr, start_dmr = dmr$start, end_dmr = dmr$end, simplified = F) 

dmr_genomeplot <- function(name, chr, betas, annotation.gr, factor_interest, start_dmr, end_dmr, flanks, simplified = T, dmps = F, genome_version = "hg19", diff_symbol = F, dotsize = 0.7){

  #Input check
  if(is.null(chr)) stop("No chr provided")
  if(!is.character(chr)) stop("chr is not a character, please input it as \"chr[1-9]\"")
  if(is.null(betas)) stop("No Beta matrix provided")
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
  if(is.null(name)){
    if(end_dmr == start_dmr){
      name <- "DMP"
    } else{
      name <- "DMR"
    }
  }
  
  #Check DMR
  dmr.hg19.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_dmr & start(annotation.gr) <= end_dmr),]
  #print(dmr.hg19.annotation)
  
  #DMR track: Useful indication of where the DMR was found
  if(start_dmr != end_dmr){
    dmrtrack <- AnnotationTrack(start = start_dmr, end = end_dmr, name = "DMR", chromosome = chr, genome = genome_version)
  } else{
    dmrtrack <- AnnotationTrack(start = start_dmr, end = end_dmr, name = "DMP", chromosome = chr, genome = genome_version)
  }
  
  #DMP track: Useful indication where the individual CpGs lie within a DMR when the number of CpGs within a DMR become large
  if(dmps){
    dmp.hg19.annotation <- dmr.hg19.annotation
    
    #Add the G and direction of each CpG
    start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"]) <- start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"])-1
    end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"]) <- end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"])+1
    dmptrack <- AnnotationTrack(dmp.hg19.annotation, name = "DMP", col = dmp.hg19.annotation$red_end)
  }
  
  #Add flanks
  start_dmr_flanks <- start_dmr-flanks
  end_dmr_flanks <- end_dmr+flanks
  dmr.flanks.annotation <- annotation.gr[which(as.character(seqnames(annotation.gr)) == chr & start(annotation.gr) >= start_dmr_flanks & start(annotation.gr) <= end_dmr_flanks),]
  
  #Betas
  GMset.beta <- as.data.frame(betas)
  dmr.450k <- data.frame(seqnames = as.character(seqnames(dmr.flanks.annotation)), pos = start(dmr.flanks.annotation), GMset.beta[names(dmr.flanks.annotation),])
  colnames(dmr.450k) <- c("seqnames", "pos", colnames(GMset.beta))
  dmr.450k.gr <- makeGRangesFromDataFrame(df = dmr.450k, keep.extra.columns = T, start.field = "pos", end.field = "pos")

  #############
  # Plot:Gviz #
  #############
  
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
  grtrack <- UcscTrack(track = "RefSeq Genes", 
                       rot.title = 0,
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
                       name = "RefSeq Genes", 
                       showId = T, 
                       from = start_dmr_flanks, 
                       to = end_dmr_flanks, 
                       fill = "#8282d2")
  
  if(simplified == F){
    cpgtrack <- UcscTrack(genome = genome_version, 
                          chromosome = chr,
                          track = "cpgIslandExt", 
                          from = start_dmr_flanks, 
                          to = end_dmr_flanks,
                          trackType = "AnnotationTrack", 
                          start = "chromStart",
                          end = "chromEnd", 
                          id = "name", 
                          shape = "box",
                          fill = "#006400",
                          name = "CpG Islands")
    txntrack <- UcscTrack(genome = genome_version, 
                          chromosome = chr,
                          track = "Txn Factor ChIP", 
                          table = "wgEncodeRegTfbsClusteredV3", 
                          from = start_dmr_flanks, 
                          to = end_dmr_flanks,
                          trackType = "AnnotationTrack", 
                          start = "chromStart",
                          end = "chromEnd", 
                          feature = "name", 
                          shape = "box",
                          fill = "#960000",
                          name = "TFBS",
                          featureAnnotation = "feature")
    h3k27actrack <- UcscTrack(genome = genome_version, 
                              chromosome = chr,
                              track = "Layered H3K27Ac", 
                              table = "wgEncodeBroadHistoneGm12878H3k27acStdSig",
                              from = start_dmr_flanks, 
                              to = end_dmr_flanks,
                              trackType = "DataTrack", 
                              start = "start", 
                              end = "end", 
                              data = "score", 
                              type = "hist",
                              window = -1,
                              windowSize = 50,
                              fill.histogram = "black",
                              name = "H3K27Ac", 
                              ylim = c(0,100))
    h3k4me1track <- UcscTrack(genome = genome_version, 
                              chromosome = chr,
                              track = "Layered H3K4Me1", 
                              table = "wgEncodeBroadHistoneGm12878H3k4me1StdSig",
                              from = start_dmr_flanks, 
                              to = end_dmr_flanks,
                              trackType = "DataTrack", 
                              start = "start", 
                              end = "end", 
                              data = "score", 
                              type = "hist",
                              window = -1,
                              windowSize = 50,
                              fill.histogram = "black",
                              name = "H3K4Me1", 
                              ylim = c(0,100))
    h3k4me3track <- UcscTrack(genome = genome_version, 
                              chromosome = chr,
                              track = "Layered H3K4Me3", 
                              table = "wgEncodeBroadHistoneGm12878H3k4me3StdSig",
                              from = start_dmr_flanks, 
                              to = end_dmr_flanks,
                              trackType = "DataTrack", 
                              start = "start", 
                              end = "end", 
                              data = "score", 
                              type = "hist",
                              window = -1,
                              windowSize = 50,
                              fill.histogram = "black",
                              name = "H3K4Me3", 
                              ylim = c(0,100))
    dnasetrack <- UcscTrack(genome = genome_version, 
                            chromosome = chr,
                            track = "DNase Clusters", 
                            table = "wgEncodeRegDnaseClusteredV3",
                            from = start_dmr_flanks, 
                            to = end_dmr_flanks,
                            trackType = "AnnotationTrack", 
                            start = "chromStart", 
                            end = "chromEnd", 
                            data = "score", 
                            shape = "box",
                            name = "DNase")
    mnasetrack <- UcscTrack(genome = genome_version, 
                            chromosome = chr,
                            track = "DNase Clusters", 
                            table = "wgEncodeRegDnaseClusteredV3",
                            from = start_dmr_flanks, 
                            to = end_dmr_flanks,
                            trackType = "AnnotationTrack", 
                            start = "chromStart", 
                            end = "chromEnd", 
                            data = "score", 
                            shape = "box",
                            name = "DNase")
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
                          name = "SNPs")
  }
  
  factor_interest <- as.factor(factor_interest)
  type_plot <- c("a", "p", "g")
  #type_plot <- c("a", "p", "g", "confint")
  
  #Change symbols if necessary
  if(diff_symbol == T){
    pchp <- as.numeric(as.factor(levels(factor_interest)))
  } else{
    pchp <- 20
  }
  #Change size symbols
  cexp <- dotsize
  
  dtrack.450k <- DataTrack(range = dmr.450k.gr, 
                           name = "450K", 
                           ylim = c(0,1), 
                           groups = factor_interest,
                           lty = as.numeric(as.factor(levels(factor_interest))),
                           pch = pchp,
                           cex = cexp,
                           type = type_plot, 
                           legend = TRUE)
  
  if(simplified == F){
    if(dmps){
      plotTracks(list(itrack, 
                      gtrack, 
                      grtrack, 
                      dtrack.450k, 
                      dmrtrack, 
                      dmptrack, 
                      snptrack), 
                 from = start_dmr_flanks, 
                 to = end_dmr_flanks, 
                 main = name,
                 cex.title = 0.9,
                 cex.axis = 1,
                 add53 = T,
                 add35 = T,
                 littleTicks = F,
                 fontsize = 12, 
                 fontcolor = "black",
                 background.title = "darkgray") 
    } else{
      plotTracks(list(itrack, 
                      gtrack, 
                      grtrack, 
                      dtrack.450k, 
                      dmrtrack,
                      snptrack), 
                 from = start_dmr_flanks, 
                 to = end_dmr_flanks, 
                 main = name,
                 cex.title = 0.9,
                 cex.axis = 1,
                 add53 = T,
                 add35 = T,
                 littleTicks = F,
                 fontsize = 12, 
                 fontcolor = "black",
                 background.title = "darkgray") 
    }
  } else{
    if(dmps){
      plotTracks(list(itrack, 
                      gtrack, 
                      grtrack, 
                      dtrack.450k,
                      dmrtrack, 
                      dmptrack), 
                 from = start_dmr_flanks, 
                 to = end_dmr_flanks, 
                 main = name,
                 cex.title = 0.9,
                 cex.axis = 1,
                 add53 = T,
                 add35 = T,
                 littleTicks = F,
                 fontsize = 12, 
                 fontcolor = "black",
                 background.title = "darkgray")  
    } else{
      plotTracks(list(itrack, 
                      gtrack, 
                      grtrack, 
                      dtrack.450k,
                      dmrtrack), 
                 from = start_dmr_flanks, 
                 to = end_dmr_flanks, 
                 main = name,
                 cex.title = 1.2,
                 cex.axis = 1,
                 add53 = T,
                 add35 = T,
                 littleTicks = F,
                 fontsize = 12, 
                 fontcolor = "black",
                 background.title = "darkgray")  
    }
  }
  
  return(dmr_annotation = dmr.hg19.annotation)
}