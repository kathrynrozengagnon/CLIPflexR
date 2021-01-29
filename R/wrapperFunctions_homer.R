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
#' @param style Type of Homer peak calling, factor | histone | groseq | tss | dnase | super | mCs
#' @param foldEnrichmentOverLocal fold enrichment over local tag count
#' @param localSize region to check for local tag enrichment
#' @param strand find peaks using tags on both | separate strands
#' @param minDist minimum distance between peaks
#' @param size Peak size
#' @param genomeSize genome size, NULL (default) = 2E9, applicable for human or mouse), or set for your genome
#' @param fragLength Approximate fragment length
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen.
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @return Path to unzipped file
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
#' FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
#' bam <- bowtie_align(FqFile_QFColStpClipped,myIndex)
#' bed <- bamtobed(bam)
#' homer_peaks(bed)
#' @export
homer_peaks <- function(fileTofqs, maketagdir = "makeTagDirectory", findpeaks = "findpeaks", 
                        format = "bed", createSingleTagsTSV = TRUE, tagdir = file.path(dirname(fileTofqs), 
                                                                                       gsub("\\.bed", "", make.names(basename(fileTofqs)))), 
                        style = "factor", foldEnrichmentOverLocal = 2, localSize = 10000, 
                        strand = "seperate", minDist = 50, size = 10, fragLength = 10, genomeSize =  NULL, stderr = file.path(dirname(fileTofqs),paste0(basename(fileTofqs), "_homer_stderr.txt")), 
                        stdout = file.path(dirname(fileTofqs),paste0(basename(fileTofqs), "_homer_stderr.txt")), useClipRConda = ifelse(is.null(getOption("CLIPflexR.condaEnv")), 
                                                                                                                                        FALSE, TRUE), additionalArgumements = NULL, verbose = FALSE, writelog  = T) {
  cmd <- maketagdir
  if (useClipRConda) 
    cmd <- file.path(getOption("CLIPflexR.condaEnv"), "bin", 
                     cmd)
  if (!file.exists(fileTofqs)) 
    stop("File does not exist")
  baseNAME <- make.names(basename(fileTofqs))
  tagDir <- tagdir
  if (file.exists(fileTofqs) & !dir.exists(tagDir)) {
    dir.create(tagDir, showWarnings = TRUE, recursive = TRUE)
    formatTagIn <- paste0("-format ", format)
    args <- c(paste0(tagDir), paste0(fileTofqs), ifelse(createSingleTagsTSV, 
                                                        "-single", NULL), formatTagIn)
    if (verbose) {
      message("makeTagDirectory command is ", cmd)
      message("makeTagDirectory arguments are ", paste0(args, 
                                                        sep = " ", collapse = " "))
    }
    if(writelog) {
      system2(cmd, args, stdout = stdout, stderr = stderr)}  else {
        system2(cmd, args, stdout = "", stderr = "")
      }
  }
  cmd <- findpeaks
  if (useClipRConda) 
    cmd <- file.path(getOption("CLIPflexR.condaEnv"), "bin", cmd)
  if (dir.exists(tagDir) & !file.exists(file.path(tagDir, "peaks.txt")) & is.null(genomeSize)){
    args <- c(paste0(tagDir), paste0("-o auto"), paste0("-style ", style), paste0("-L ", foldEnrichmentOverLocal), 
              paste0("-localSize ", localSize), paste0("-strand ", strand), paste0("-minDist ", minDist),
              paste0("-size ", size), paste0("-fragLength ", fragLength))
  } 
  else if (dir.exists(tagDir) & !file.exists(file.path(tagDir,"peaks.txt")) & !is.null(genomeSize)) {
    args <- c(paste0(tagDir), paste0("-o auto"), paste0("-style ", style), paste0("-L ", foldEnrichmentOverLocal), 
              paste0("-localSize ", localSize), paste0("-strand ", strand), paste0("-minDist ", minDist), 
              paste0("-size ", size), paste0("-fragLength ", fragLength), paste0("-gsize ", genomeSize))
  }
  if (verbose) {
    message("findPeaks command is ", cmd)
    message("findPeaks arguments are ", paste0(args, sep = " ", collapse = " "))
  }
  if(writelog) {
    system2(cmd, args, stdout = stdout, stderr = stderr)}  else {
      system2(cmd, args, stdout = "", stderr = "")
    }
  return(file.path(tagDir, "peaks.txt"))
}
