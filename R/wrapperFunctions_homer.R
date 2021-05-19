#' Wrapper function for Homer's makeTagDirectory and findPeaks
#'
#' Wrapper function for Homer's makeTagDirectory and findPeaks
#'
#'
#' @docType methods
#' @name homer_peaks
#' @rdname homer_peaks
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqs path to file to process (BED).
#' @param maketagdir path to makeTagDirectory from Homer toolkit.
#' @param findpeaks path to findpeaks from Homer toolkit.
#' @param format input format, "bed" (default).
#' @param createSingleTagsTSV create a single tags.tsv file for all "chromosomes".
#' @param tagdir name of Tag Directory.
#' @param style type of Homer peak calling, "factor", "histone", "groseq", "tss", "dnase", "super" or "mCs".
#' @param foldEnrichmentOverLocal fold enrichment over local tag count.
#' @param localSize region to check for local tag enrichment.
#' @param strand find peaks using tags on "both" or "separate" (default) strands.
#' @param minDist minimum distance between peaks.
#' @param size peak size.
#' @param genomeSize genome size, default is NULL (2E9, applicable for human or mouse), set integer to specify for your genome.
#' @param fragLength approximate fragment length.
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @return path to unzipped file
#' @examples
#' \dontrun{
#' testFasta <- system.file("extdata/hg19Small.fa",package = "CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package = "CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite = TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
#' FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStpClipped,myIndex, overwrite = TRUE))
#' bed <- bamtobed(bam)
#' homer_peaks(bed)
#' }
#' @export
homer_peaks <- function(fileTofqs, maketagdir = "makeTagDirectory", findpeaks = "findPeaks", 
                        format = "bed", createSingleTagsTSV = TRUE, tagdir = file.path(dirname(fileTofqs), 
                                                                                       gsub("\\.bed", "", make.names(basename(fileTofqs)))), 
                        style = "factor", foldEnrichmentOverLocal = 2, localSize = 10000, 
                        strand = "seperate", minDist = 50, size = 10, fragLength = 10, genomeSize =  NULL, stderr = file.path(dirname(fileTofqs),paste0(basename(fileTofqs), "_homer_stderr.txt")), 
                        stdout = file.path(dirname(fileTofqs),paste0(basename(fileTofqs), "_homer_stderr.txt")), useClipRConda = ifelse(is.null(getOption("CLIPflexR.condaEnv")), 
                                                                                                                                        FALSE, TRUE), additional_Args = NULL, verbose = FALSE, writelog  = T) {
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
