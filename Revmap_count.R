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
#' unbam(bam)
#' @return Path to FASTA file
#' @import reticulate GenomicAlignments Biostrings
#' @export
revmap_process <- function(input, linkers = NULL, length_max = NULL, length_min = NULL) { 
  require(Biostrings)
  require(IRanges)
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
  outname = gsub(".fa", "_processed.fa", outname)
  writeXStringSet(fa, outname)
  return(outname)
   }
  
  
  #' Merge read sequence with bedfile
  #'
  #'Merge read sequence with bedfile
  #'
  #'
  #' @docType methods
  #' @name chimera_joinread
  #' @rdname chimera_joinread
  #'
  #' @author Kathryn Rozen-Gagnon
  #'
  #' @param file File to process.
  #' @param outFile Name of output file.
  #' @param filtDup Output index name
  #' @examples
  #' bams <- system.file("extdata/example_bams/hg19Small.fa",package="CLIPflexR")
  #' genomeIndex <- bowtie2_index(system.file("extdata/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa",package="CLIPflexR"))
  #' knownMiRNAs <- system.file("extdata/aae_miRNAs_mature_fixed.fa",package="CLIPflexR")
  #' exclude <- decompress(testFQ,overwrite=TRUE)
  #' FqFile_QF <- fastq_quality_filter(FqFile)
  #' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
  #' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
  #' FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
  #' bam <- bowtie_align(FqFile_QFColStpClipped,myIndex)
  #' bamtobed(bam)
  #' @return Path 
  #' @import GenomicAlignments
  #' @importMethodsFrom rtracklayer export.bed export.bw mcols
  #' @export
  #' 
  
  
  
  
  
  #' Merge read sequence with bedfile
  #'
  #'Merge read sequence with bedfile
  #'
  #'
  #' @docType methods
  #' @name chimera_process
  #' @rdname chimera_process
  #'
  #' @author Kathryn Rozen-Gagnon
  #'
  #' @param fastas FASTA files to reverse map.
  #' @param knownMiRNAs knownMiRNAs FASTA file.
  #' @param linkers contaminating sequences to remove 
  #' @return Path 
  #' @import GenomicAlignments BiocParallel stringr Biostrings
  #' @importMethodsFrom rtracklayer export.bed export.bw mcols
  #' @export
  revmap_count <- function(fastas, knownMiRNAs, bpparam=NULL,verbose=TRUE, linkers = NULL, length_max = NULL, length_min = NULL, removedups =FALSE){
    if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
    if(verbose) message("Processing fastas..",appendLF = FALSE)
    pro_fastas <- vector("list",length = length(fastas)) 
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }
    if (verbose) {
    if (!is.null(linkers)) message("removing linkers... ", linkers)
    if (!is.null(length_max)) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min)) message ("getting sequences longer than ", as.character(length_min), " nt")  }
    if(verbose) message("Creating indices from FASTA files..",appendLF = FALSE)
    indicies <- bplapply(pro_fastas,bowtie2_index,BPPARAM=bpparam)
    if(verbose) message("done")
    if(verbose) message("Mapping miRNAs to processed reads..",appendLF = FALSE)
    revBams <- bplapply(indicies,
                        function(x,knownMiRNAs){
                          bowtie_align(knownMiRNAs,index=x, bam = paste0(x, ".bam"))
                        },knownMiRNAs=knownMiRNAs,BPPARAM=bpparam)
    if(verbose) message("done")
    if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
    beds <- bplapply(revBams, bamtobed,BPPARAM=bpparam)
  if (removedups) { message("deduplicating...") 
      dedup <- lapply(beds,read.delim,  header = F)
      for (i in 1:length(dedup)) {
      dedup[[i]]$miRNAnum  <- ifelse(!grepl("miR|let|iab", dedup[[i]]$V4), str_sub(dedup[[i]]$V4, 2, 7), NA)
      dedup[[i]]$miRNAnum <- ifelse(grepl("miR|let|iab", dedup[[i]]$V4), as.numeric(apply(dedup[[i]],1,function(x) regmatches(x["name"],regexpr("[0-9]+",x["name"])))), paste(dedup[[i]]$miRNAnum))
      dedup[[i]] <- dedup[[i]][order(dedup[[i]]$V1, dedup[[i]]$V2, dedup[[i]]$miRNAnum,dedup[[i]]$V4),]
      dedup[[i]] <- dedup[[i]][!duplicated(dedup[[i]]$V1),] 
      dedup[[i]]$miRNAnum <- NULL
      } 
    names(dedup) <- gsub(".bed", "_deduplicated.bed", beds)
    for (i in seq_along(dedup)) {
      write.table(dedup[i],names(dedup)[i], col.names =  FALSE, row.names= FALSE, sep = "\t", quote = F)
      } 
    return(names(dedup))[i]} 
    else { return(beds)}
    if(verbose) message("done")  }
    
   
  