% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\docType{methods}
\name{extract_unmapped}
\alias{extract_unmapped}
\title{Extracts unmapped reads from BAM and writes them to a FASTA file}
\usage{
extract_unmapped(bam, outfa = NULL)
}
\arguments{
\item{bam}{path to BAM file.}

\item{outfa}{path to output FASTA (default is same directory as input BAM).}
}
\value{
path to FASTA file of unmapped reads.
}
\description{
Extracts unmapped reads from BAM and writes them to a FASTA file
}
\examples{
testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
myIndex <- suppressWarnings(bowtie2_index(testFasta, overwrite=TRUE))
testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20")
FqFile <- decompress(FqFile_FF,overwrite=TRUE)
FqFile_clipped <- fastx_clipper(FqFile,length=20)
FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
FqFile_Col <- ctk_fastq2collapse(FqFile_QF)
FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5, inputFormat="fastq")
bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, 
overwrite=TRUE, inputFormat="fastq"))
extract_unmapped(bam)
}
\author{
Kathryn Rozen-Gagnon
}
