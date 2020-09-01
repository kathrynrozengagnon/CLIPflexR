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
#' @param bed Bed file
#' @param GR GR need add description
#' @param notStranded Stranded or not stranded
#' @param interFeature interfeature # need add description
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' bowtie2_index(testFasta)
#' @import rtracklayer GenomicRanges GenomicAlignments
#' @return counting matrix
#' @export
countFromBed <- function(Bed,GR,notStranded=TRUE,interFeature=FALSE){
  # require(rtracklayer)
  # require(GenomicRanges)
  # require(GenomicAlignments)
  reads <- import.bed(Bed)
  fk <- summarizeOverlaps(GR,reads,ignore.strand = notStranded,inter.feature=interFeature)
  assay(fk)
}
