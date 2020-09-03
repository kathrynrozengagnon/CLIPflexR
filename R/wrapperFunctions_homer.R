#' Wrapper function for Homer's makeTagDirectory and findpeaks
#'
#' Wrapper function for Homer's makeTagDirectory and findpeaks
#'
#'
#' @docType methods
#' @name homer_peaks
#' @rdname homer_peaks
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqs File to process.
#' @param maketagdir Path to makeTagDirectory from Homer toolkit
#' @param findpeaks Path to findpeaks from Homer toolkit
#' @param format Input format
#' @param createSingleTagsTSV Create a single tags.tsv file for all "chromosomes"
#' @param tagdir Name of Tag Directory
#' @param style Type of Homer peak calling
#' @param foldEnrichmentOverLocal fold enrichment over local tag count
#' @param localSize region to check for local tag enrichment
#' @param strand find peaks using tags on both strands or separate
#' @param minDist minimum distance between peaks
#' @param size Peak size
#' @param genomeSize genome size
#' @param fragLength Approximate fragment length,
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
#' FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
#' bam <- bowtie_align(FqFile_QFColStpClipped,myIndex)
#' bed <- bamtobed(bam)
#' homer_peaks(bed)
#' @export
homer_peaks <- function(fileTofqs,maketagdir="makeTagDirectory",
                        findpeaks="findpeaks",
                        format="bed",
                        createSingleTagsTSV=TRUE,
                        tagdir=file.path(dirname(fileTofqs),
                                         gsub("\\.bed","",make.names(basename(fileTofqs)))),
                        style="factor",
                        foldEnrichmentOverLocal=2,
                        localSize=10000,
                        strand="seperate",
                        minDist=50,
                        size=10,
                        fragLength=10,
                        genomeSize=NULL, 
                        stderr=paste0(getwd(),"homer_stats_stderr"),
                        stdout=paste0(getwd(),"homer_stats_stdout"),
                        useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                        additionalArgumements=NULL,verbose=FALSE){
  cmd <- maketagdir
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqs))stop("File does not exist")
  baseNAME <- make.names(basename(fileTofqs))
  tagDir <- tagdir
  
  if(file.exists(fileTofqs) & !dir.exists(tagDir)){
    
    dir.create(tagDir,showWarnings = TRUE,recursive = TRUE)
    formatTagIn <- paste0("-format ",format)
    
    args <- c(
      paste0(tagDir),
      paste0(fileTofqs),
      ifelse(createSingleTagsTSV,"-single",NULL),
      formatTagIn
    )
    if(verbose){      
      message("makeTagDirectory command is ",cmd)
      message("makeTagDirectory arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    
    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
  }
  
  cmd <- findpeaks
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(dir.exists(tagDir) & !file.exists(file.path(tagDir,"peaks.txt"))){
    
    # style="factor",
    # foldEnrichmentOverLocal=2,
    # localSize=10000,
    # strand="seperate",
    # minDist=50,
    # size=10,
    # fragLength=10,
    # 
    args <- c(
      paste0(tagDir),
      paste0("-o auto"),
      paste0("-style ",style),
      paste0("-L ",foldEnrichmentOverLocal),
      paste0("-localSize ",localSize),
      paste0("-strand ",strand),      
      paste0("-minDist ",minDist),
      paste0("-size ",size),
      paste0("-fragLength ",fragLength),
      paste0("-gsize ", genomeSize))
    if(verbose){      
      message("findPeaks command is ",cmd)
      message("findPeaks arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    
    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
    
  }
  return(file.path(tagDir,"peaks.txt"))
}
