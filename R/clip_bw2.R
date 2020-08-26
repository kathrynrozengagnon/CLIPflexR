#' Build bigwigs
#'
#' @rdname ClIP_bw2
#' @param sort_bam sorted bam file
#' @param res_dir result directory
#' @param normalized Normalised to total reads.
#' 
#' @docType methods
#' @import  Rsamtools
#' @return bigwig
#' @export
ClIP_bw2 <- function(sort_bam,res_dir=NULL,normalized=TRUE){
  if(dir.exists(res_dir)){dir.create(res_dir)}else{print(paste0("folder is exist: ",res_dir))}
  samID <- gsub("_sort.bam","",basename(sort_bam))
  bw_out <- file.path(res_dir,paste0(samID,".bigwig"))
  if(isTRUE(normalized)){
    allChromosomeStats <- idxstatsBam(sort_bam)
    mapped <- sum(allChromosomeStats[,"mapped"])
    toRun <- coverage(BamFile(sort_bam),weight = (10^6)/mapped)
    export.bw(toRun,con=bw_out)
  }else{
    toRun <- coverage(BamFile(sort_bam))
    export.bw(toRun,con=bw_out)}
}