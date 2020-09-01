#' make count matrix from bed
#'
#' make count matrix from bed
#'
#'
#' @docType methods
#' @name countFromBed
#' @rdname countFromBed
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param bed mapped reads, bed files.
#' @param GR count reads over these genomic ranges. 
#' @param notStranded Stranded or not stranded TRUE or FALSE.
#' @param interFeature Count reads mapping to multiple features TRUE or FALSE.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' bowtie2_index(testFasta)
#' @import rtracklayer GenomicRanges GenomicAlignments
#' @return counting matrix
#' @export
countFromBed <- function(Bed,GR,notStranded=TRUE,interFeature=FALSE){
  reads <- import.bed(Bed)
  fk <- summarizeOverlaps(GR,reads,ignore.strand = notStranded,inter.feature=interFeature)
  assay(fk)
}