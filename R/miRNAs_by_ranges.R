#' Process reads for miRNA counting
#'
#' Process reads for miRNA counting
#'
#'
#' @docType methods
#' @name revmap_process
#' @rdname revmap_process
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fasta FASTA file to process
#' @param readlength length requirements for FASTA
#' @param linkers contaminating seqences to remove
#' @param outfa PATH to output FASTA
#' @examples
#' bam <- system.file("extdata/example_bams/KRG121817A_Aag2_Ago1_Exp1A_rm5_rm3.bam",package="CLIPflexR")
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped,paste0(FqFile_QF,".gz"))
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5)
#' bam <- bowtie_align(FqFile_QFColStripped,myIndex)
#' 
#' @return Path to FASTA file
#' @import reticulate GenomicAlignments Biostrings
#' @export
#' 
revmap_process <- function(input, linkers = NULL, length_max = NULL, length_min = NULL) { 
  require(Biostrings)
  require(IRanges)
  require(magrittr)
  fa <- readDNAStringSet(input, format = "fasta", nrec = -1L)
  if (!is.null(length_max)) {
    fa <- fa[width(fa) <= length_max]}
  if (!is.null(length_min)) {
    fa <- fa[width(fa) >= length_min]}
  if (!is.null(linkers)) {
    fa2  <- vector("list", length(linkers))
    for(x in 1:length(linkers)){
      fa2[[x]] <- vmatchPattern(pattern=DNAString(linkers[[x]]),subject=fa) %>% unlist()}
    fa2 <-  unlist(IRangesList(fa2))
    fa2 <- fa2@NAMES
    fa <- fa[!names(fa) %in% fa2]  }
  outname = paste(input, sep = "")
  outname = gsub(".fa$", "_processed.fa", outname)
  writeXStringSet(fa, outname)
  return(outname)
}

#' @docType methods
#' @name chimera_process
#' @rdname chimera_process
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas FASTA files to process.
#' @param miRNA_ranges knownMiRNAs FASTA file.
#' @param linkers contaminating sequences to remove 
#' @param length_min minimum length required during fasta processing
#' @param length_max maximum length required during fasta processing 
#' @param genomeIndex Full genome index 
#' @return Path 
#' @import GenomicAlignments BiocParallel stringr rtracklayer GenomicRanges
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @export
#' 
Ranges_count <- function(fastas,miRNA_ranges,genomeIndex,linkers = NULL, length_max = NULL, length_min = NULL, bpparam=NULL,verbose=TRUE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(!is.null(linkers) | !is.null(length_max) | !is.null(length_min)) { 
    if (verbose) message("Processing fastas..",appendLF = FALSE)
    if (!is.null(linkers)) message("removing linkers... ", linkers)
    if (!is.null(length_max)) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min)) message ("getting sequences longer than ", as.character(length_min), " nt")
    pro_fastas <- vector("list",length = length(fastas)) 
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }} else {
      pro_fastas <-  fastas
    }
  if(verbose) message("Mapping reads to genome..",appendLF = FALSE)
  mappedBams <- bplapply(pro_fastas, bowtie_align, index = genomeIndex, maxMismatches=0,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
  Mybeds <- bplapply(mappedBams, bamtobed,BPPARAM=bpparam)
  if(verbose) message("done") 
  if(verbose) message("Importing ranges to search..",appendLF = FALSE)
  qranges <- rtracklayer::import(miRNA_ranges)
  if(verbose) message("Counting mapped reads..",appendLF = FALSE)
  kks <- bplapply(Mybeds,countFromBed,GR=qranges,notStranded=FALSE, BPPARAM=bpparam)
  kksMat <- do.call(cbind,kks)
  colnames(kksMat) <- basename(unlist(Mybeds))
  kksMat <- as.data.frame(kksMat)
  outname <- unique(dirname(unlist(Mybeds)))
  outname <- paste0(outname, "/miRNA_count_matrix.txt")
  write.table(kksMat, outname,  col.names = TRUE,  row.names = FALSE,  sep =  "\t", quote = F)
  return(outname)
  if (verbose) message("done")
  }



