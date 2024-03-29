% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\docType{methods}
\name{CLIP_bw2}
\alias{CLIP_bw2}
\title{Build bigwigs}
\usage{
CLIP_bw2(
  sort_bam,
  res_dir = dirname(sort_bam),
  normalized = TRUE,
  stranded = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{sort_bam}{path to sorted BAM file.}

\item{res_dir}{result directory, default is directory where BAM file is located, to specify other/create new directory, enter path as character string.}

\item{normalized}{normalized to total reads.}

\item{stranded}{Whether to make separate bigwigs for each strand, TRUE or FALSE (default).}

\item{verbose}{print messages, TRUE or FALSE (default).}
}
\value{
bigwig.
}
\description{
Build bigwigs
}
\examples{
testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
FqFile <- decompress(testFQ,overwrite = TRUE)
FqFile_FF <- ctk_fastqFilter(FqFile,qsFilter = "mean:0-29:20",verbose=TRUE)
FqFile_clipped <- fastx_clipper(FqFile_FF,length=20)
FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose = TRUE)
FqFile_ColStrip <- ctk_stripBarcode(FqFile_Col,linkerlength=5, inputFormat="fastq") 
##map reads to genome
mapped <- suppressWarnings(bowtie_align(FqFile_ColStrip,myIndex, 
mode="genome_map", inputFormat="fastq"))
wig <- CLIP_bw2(mapped)
}
