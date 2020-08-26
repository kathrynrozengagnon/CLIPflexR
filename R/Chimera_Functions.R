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
#' @param BAMs BAM files to process.
#' @param knownMiRNAs knownMiRNAs FASTA file.
#' @param genomeIndex Full genome index
#' @return Path 
#' @import GenomicAlignments
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @export
chimera_Process <- function(bams,knownMiRNAs,genomeIndex,bpparam=NULL,verbose=TRUE){
  if(!is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(verbose) message("Extracting unmapped reads to FASTA..",appendLF = FALSE)
  fastas <- bplappy(bams,unbam,bpparam=bpparam)
  if(verbose) message("done")
  if(verbose) message("Creating indicies from FASTA files..",appendLF = FALSE)
  indicies <- bplappy(fastas,bowtie2_index,bpparam=bpparam)
  if(verbose) message("done")
  if(verbose) message("Mapping miRNAs to unmapped reads..",appendLF = FALSE)
  newBams <- bplappy(indicies,
                      function(x,knownMiRNAs){
                        bowtie_align(knownMiRNAs,index=x)
                      },knownMiRNAs=knownMiRNAs,bpparam=bpparam)
  if(verbose) message("done")
  if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
  beds <- bplapply(newBams, bamtobed,bpparam=bpparam)
  if(verbose) message("done")  
  
  
  # Not sure here...
  chimera <- bplapply(beds, function(x) {
    if (!file.size(x) == 0) {
      read.delim(x, header = FALSE, sep = "")
    }
  },bpparam=bpparam)
  names(chimera) <- beds
  t <- lapply(chimera, nrow) 
  t <- unlist(t) 
  t <- names(t) #get names of beds with entried and read these in 
  
  
  test <- lapply(t, read.delim, header = FALSE, sep = "")
  names(test) <- t 
  col.names <- c("rowname", "start", "stop", "name", "score", "strand")
  test <- lapply(test, setNames, col.names)
  t <- gsub(".bed", ".fa", t) #so t has names of all fastas to read in 
 
  fa <- fastas
  
  fasta <- bplapply(fa, readDNAStringSet, format = "fasta", nrec = -1L,bpparam=bpparam)
  names(fasta) <- fa 
  fasta <- fasta[names(fasta) %in% t]
  
  # fasta <- lapply(fasta, as.data.frame)
  fasta <- bplapplylapply(fasta, function(x) tibble::rownames_to_column(as.data.frame(x)),bpparam=bpparam)
  
  
  #merge together read sequence and bed by rowname (read name)
  chimera <- map2(fasta, test, ~merge(.x,.y, by = "rowname"))
  # names(chimera) <- gsub("/rugpfs/fs0/rice_lab/scratch/krozen/AaegL5_mapped/unmapped_for_chimera/unmapped_", "/rugpfs/fs0/rice_lab/scratch/krozen/AaegL5_mapped/unmapped_for_chimera/known_novel_sRNA_revmapped/unmapped_", names(chimera))
  names(chimera) <- names(chimera)
  
  #write files to chimera inputs: bed with read name, start, stop, srand, read name, and read sequence
  # setwd(Dir)
  fastaTxts <- vector("list",length = length(chimera))
  for (x in names(chimera)) {
    write.table(chimera[[x]], file=paste0(x,".txt"), sep="\t", quote = FALSE)
    fastaTxts[[x]] <- paste0(x,".txt")
  }
 
  ####
  # Eh?
  ##
  # write.table(test, file=".txt", sep="\t")
  
  
  chimeraProcess <- function(input) {
    
    BR1 <- read.delim(input, header=T)
    BR1 <- BR1[!grepl("AAAAAAAAAAAAAAAAAAAAAAAAA", BR1$name),] #this was high frequency in my libraries removed
    BR1<- BR1[BR1$strand=="+",]
    BR1 <- BR1[duplicated(BR1$rowname)==F,] #come back to this, to remove reads with more than one miRNA mapped 
    BR1$length <- sapply(as.character(BR1$x), nchar)
    BR1$ups.seq <- mapply(substr, x=BR1$x, start=0, stop=BR1$start)
    BR1$dns.seq <- mapply(substr, x=BR1$x, start=BR1$stop, stop=BR1$length)
    outname = paste(input, sep = "")
    outname = gsub(".fa.txt", "_chimera.txt", outname)
    write.table(BR1, outname, quote=F, sep="\t")
    return(outname)
  }
  files <- vector("list",length = length(fastaTxts))
  for (i in 1:length(fastaTxts)) {
    files[[i]] <- chimeraProcess(fastaTxts[i])
  }
  rffiles <- vector("list",length = length(files))  
  for (i in 1:length(files)) {
    rffiles[[i]] <- reformat(files[i])
  }
  
  
  
  remappedBams <- bplapply(rffiles, bowtie_align, index =genomeIndex,bpparam=bpparam)
  
  myBeds <- bplapply(remappedBams, bamtobed,bpparam=bpparam)
  
  
  return(myBeds)
  
}

chimeraProcess <- function(input) {
  
  BR1 <- read.delim(input, header=T)
  BR1 <- BR1[!grepl("AAAAAAAAAAAAAAAAAAAAAAAAA", BR1$name),] #this was high frequency in my libraries removed
  BR1<- BR1[BR1$strand=="+",]
  BR1 <- BR1[duplicated(BR1$rowname)==F,] #come back to this, to remove reads with more than one miRNA mapped 
  BR1$length <- sapply(as.character(BR1$x), nchar)
  BR1$ups.seq <- mapply(substr, x=BR1$x, start=0, stop=BR1$start)
  BR1$dns.seq <- mapply(substr, x=BR1$x, start=BR1$stop, stop=BR1$length)
  outname = paste(input, sep = "")
  outname = gsub(".fa.txt", "_chimera.txt", outname)
  write.table(BR1, outname, quote=F, sep="\t")
  return(outname)
}

reformat <- function(input) {
  BR1 <- read.delim(input, header=T)
  BR1$V8 <- paste(BR1$rowname, BR1$name, sep = ";")
  BR1 <- BR1[c(11,10 )]
  BR1$dns.seq <- as.character(BR1$dns.seq)
  BR1$length <- nchar(BR1$dns.seq)
  BR1 <- BR1[BR1$length >= 18,]
  BR1 <- BR1[1:2]
  outname = paste(input, '.fa', sep = "")
  BR2 <- DNAStringSet(BR1$dns.seq, use.names = TRUE)
  names(BR2) <- BR1$V8
  writeXStringSet(BR2, outname, format ="fasta") 
  return(outname)
}


