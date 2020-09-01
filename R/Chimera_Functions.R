#' Turns a BAM to a FASTA file
#'
#' Turns a BAM to a FASTA file
#'
#'
#' @docType methods
#' @name unbam
#' @rdname unbam
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param bam BAM file to turn to FASTA
#' @param outfa PATH to output FASTA
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' myIndex <- bowtie2_index(testFasta)
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
unbam <- function(bam,outfa=NULL){
  require(GenomicAlignments)
  require(Rsamtools)
  require(Biostrings)
  inBAM <- scanBam(bam,param=ScanBamParam(what = c("qname","seq"),flag = scanBamFlag(isUnmappedQuery = TRUE)))
  if(is.null(outfa)) outfa <- gsub(".bam","_unmapped.fa",bam)
  toWrite <- inBAM[[1]]$seq
  names(toWrite) <- inBAM[[1]]$qname
  writeXStringSet(toWrite,filepath = outfa,append = FALSE)
  return(outfa)
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
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
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
#' @param bams BAM files to process.
#' @param knownMiRNAs knownMiRNAs FASTA file.
#' @param genomeIndex Full genome index.
#' @param exclude Names of miRNAs or sequences to remove, character vector. 
#' @param length_min minimum length required during fasta processing, integer.
#' @param length_max maximum length required during fasta processing, integer. 
#' @return Path 
#' @import GenomicAlignments BiocParallel purrr stringr
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @export
chimera_Process <- function(bams,knownMiRNAs,genomeIndex,exclude, bpparam=NULL,verbose=TRUE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(verbose) message("Extracting unmapped reads to FASTA..",appendLF = FALSE)
  fastas <- bplapply(bams,unbam,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Creating indicies from FASTA files..",appendLF = FALSE)
  indicies <- bplapply(fastas,bowtie2_index,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Mapping miRNAs to unmapped reads..",appendLF = FALSE)
 # newBams <- bplapply(indicies,
 #                      function(x,knownMiRNAs){
 #                        bowtie_align2(knownMiRNAs,index=x, sam = paste0(x, ".bam"))
 #                     },knownMiRNAs=knownMiRNAs,BPPARAM=bpparam)
 newBams <- bplapply(indicies,
 function(x,knownMiRNAs){
   bowtie_align(knownMiRNAs,index=x, bam = paste0(x, ".bam"), maxMismatches = 0, report_k = 1000000)
 },knownMiRNAs=knownMiRNAs,BPPARAM=bpparam)
  
  if(verbose) message("done")
  if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
  beds <- bplapply(newBams, bamtobed,BPPARAM=bpparam)
  if(verbose) message("done")  
  
  beds <- beds[file.size(unlist(beds))!=0]
  chimera <- lapply(beds, read.delim, header = FALSE, sep = "")
  names(chimera) <- beds 
  col.names <- c("rowname", "start", "stop", "name", "score", "strand")
  chimera <- lapply(chimera, setNames, col.names)
  fanams <- gsub(".bed", ".fa", beds)
  
  fa <- fastas[fastas %in% fanams]
  fasta <- bplapply(fa, readDNAStringSet, format = "fasta", nrec = -1L,BPPARAM=bpparam)
  names(fasta) <- fa 
  fasta <- bplapply(fasta, function(x) tibble::rownames_to_column(as.data.frame(x)),BPPARAM=bpparam)
  #check for duplicate readnames
  if (verbose) message("Checking read names..", appendLF = FALSE)
  checkreads <- lapply(fasta, function(x) duplicated(x$rowname))
  checkreads <- unlist(lapply(checkreads, function(x) length(which(x==TRUE))))
  idx <- which(checkreads > 0)
  if (verbose & isTRUE(length(idx)==0)) message("read names ok") else message(paste0("Warning, the following sample(s) have duplicate read names! \n", names(idx)))
  #merge together read sequence and bed by rowname (read name)
  chimera <- purrr::map2(fasta,chimera, ~merge(.x,.y, by = "rowname"))
  
  #write files to chimera inputs: bed with read name, start, stop, srand, read name, and read sequence
  # setwd(Dir)
  fastaTxts <- vector("list",length = length(chimera))
  for (x in 1:length(chimera)) {
    write.table(chimera[[x]], file=paste0(names(chimera)[x],".txt"), sep="\t", quote = FALSE)
    fastaTxts[[x]] <- paste0(names(chimera)[x],".txt")
  }
 names(fastaTxts) <- names(chimera)

  chimerafiles <- vector("list",length = length(fastaTxts))
  for (i in 1:length(fastaTxts)) {
    chimerafiles[[i]] <- chimeraProcess(fastaTxts[[i]], exclude)
  }
  
  rffiles <- vector("list",length = length(chimerafiles)) 

  for (i in 1:length(chimerafiles)) {
    rffiles[[i]] <- reformat(chimerafiles[[i]])
  }
  
  
  
  remappedBams <- bplapply(rffiles, bowtie_align, index =genomeIndex,BPPARAM=bpparam)
  
  myBeds <- bplapply(remappedBams, bamtobed,BPPARAM=bpparam)
  
  
  return(myBeds)
  
}

chimeraProcess <- function(input, exclude) {
  BR1 <- read.delim(input, header=T)
  BR1 <- BR1[!BR1$name %in% exclude,] 
  BR1<- BR1[BR1$strand=="+",]
    BR1$miRNAnum <- ifelse(!grepl("miR|let|iab", BR1$name), stringr::str_sub(BR1$name, 2, 7), NA)
    BR1$miRNAnum <- ifelse(grepl("miR|let|iab", BR1$name), as.numeric(apply(BR1,1,function(x) regmatches(x["name"],regexpr("[0-9]+",x["name"])))), paste(BR1$miRNAnum))
    BR1 <- BR1[order(BR1$rowname, BR1$start, BR1$miRNAnum, BR1$name),]
    BR1 <- BR1[!duplicated(BR1$rowname),] 
    BR1$miRNAnum <- NULL
  BR1$length <- sapply(as.character(BR1$x), nchar)
  BR1$ups.seq <- mapply(substr, x=BR1$x, start=0, stop=BR1$start)
  BR1$dns.seq <- mapply(substr, x=BR1$x, start=BR1$stop+1, stop=BR1$length)
  outname = paste(input, sep = "")
  outname = gsub(".fa.txt", "_chimera.txt", outname)
  write.table(BR1, outname, quote=F, sep="\t")
  return(outname)
}

reformat <- function(input) {
  BR1 <- read.delim(input, header=T)
  BR1$ID <- paste(BR1$rowname, BR1$name, sep = ";") 
  BR1 <- BR1[c("ID","dns.seq")]
  BR1$dns.seq <- as.character(BR1$dns.seq)
  BR1 <- BR1[nchar(BR1$dns.seq)>=18,] 
  outname = paste(input, '.fa', sep = "")
  BR2 <- DNAStringSet(BR1$dns.seq, use.names = TRUE)
  names(BR2) <- BR1$ID
  writeXStringSet(BR2, outname, format ="fasta")  
  return(outname)
}
