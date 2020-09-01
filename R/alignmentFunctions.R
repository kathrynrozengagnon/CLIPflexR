
#' Make index for Rbowtie2
#'
#' Make index for Rbowtie2
#'
#'
#' @docType methods
#' @name bowtie2_index
#' @rdname bowtie2_index
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param genomeFasta File to process.
#' @param outIndex Output index name
#' @param overwrite Overwrite if directory exists.
#' @param threads Number of threads to use.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' bowtie2_index(testFasta)
#' @importFrom Rbowtie2 bowtie2_build
#' @return Path to index 
#' @export
bowtie2_index <- function(genomeFasta,
                         outIndex=gsub("\\.fa","",genomeFasta),
                         overwrite=TRUE,threads=1
) {
  bowtieArgs <- paste0("--threads ",threads)
  if(!dir.exists(outIndex)){
    suppressMessages(
    bowtie2_build(references=genomeFasta,
                  bt2Index=outIndex,
                  overwrite = overwrite,bowtieArgs))
  }
  return(outIndex)
}


#' Alignment using Rbowtie2
#'
#' Alignment using Rbowtie2
#'
#'
#' @docType methods
#' @name bowtie_align
#' @rdname bowtie_align
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fq File to process.
#' @param index Index name
#' @param bam Output bam name
#' @param format Format of reads (fastq or fasta)
#' @param maxMismatches max mismatches
#' @param seedSubString length of seed substrings
#' @param threads Number of threads to use
#' @param report_k number of 
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
#' bowtie_align(FqFile_QFColStripped,myIndex)
#' @importFrom Rsamtools asBam indexBam sortBam
#' @importFrom Rbowtie2 bowtie2
#' @return Path to BAM 
#' @export
bowtie_align <- function(fq,index,
                         bam=file.path(dirname(fq),
                                       paste0("Sorted_",basename(fq),".bam")),
                         format="fasta",maxMismatches=1,seedSubString=18,threads=1,report_k=NULL
) {

    if(format == "fasta"){
      optionFormat <- "-f"
    }else{
      optionFormat <- ""
    }
  
  if(file_ext(fq) == "gz"){
      R.utils::gunzip(fq,
                    destname=gsub("\\.gz$","",fq))
    fq <- gsub("\\.gz$","",fq)
    
  }
  if(is.null(report_k)){
    bowtieArgs <- paste0(optionFormat,
                         " -N ",maxMismatches,
                         " -L ",seedSubString,
                         " --threads ",threads)
    }else{
      report_k = as.integer(report_k)
      bowtieArgs <- paste0(optionFormat,
                           " -N ",maxMismatches,
                           " -L ",seedSubString,
                           " --threads ",threads,
                           " -k ",report_k)
      }
  if(!file.exists(bam)){
    suppressMessages(bowtie2(bt2Index = index,
            samOutput = gsub("\\.bam$",".sam",
                             bam),
            seq1 = fq,
            bowtieArgs))
    asBam(gsub("\\.bam$",".sam",
               bam),
          gsub("\\.bam$",".temp",
               bam))
    sortBam(gsub("\\.bam$",".temp.bam",
                 bam),
            destination = gsub("\\.bam$","",
                               bam))
    unlink(gsub("\\.bam$",".temp.bam",
                bam)) # remove temp.bam
    indexBam(bam)
    
  }
  return(bam)
  
}
