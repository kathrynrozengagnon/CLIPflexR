% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\docType{methods}
\name{chimera_Process}
\alias{chimera_Process}
\title{Chimera process CLIP data}
\usage{
chimera_Process(
  bams,
  fastas = NULL,
  knownMiRNAs,
  genomeIndex,
  exclude = NULL,
  bpparam = NULL,
  verbose = FALSE,
  removedups = TRUE,
  overwrite = FALSE
)
}
\arguments{
\item{bams}{path to BAM files mapped to the genome, unmapped reads will be extracted for chimera processing. This requires that you mapped your reads to the genome disallowing soft clipping. If you allowed soft clipping, use fasta param and set bams = NULL.}

\item{fastas}{path to reads to process if you want to perform chimera analysis on all reads, i.e. if you mapped your reads to the genome allowing soft clipping, default is NULL.}

\item{knownMiRNAs}{path to FASTA file containing annotated miRNA or other small RNA sequence. Known miRNAs will be prioritized and known miRNA names must be in the format "miR-", "let-", "bantam-", "iab-".}

\item{genomeIndex}{path to genome index.}

\item{exclude}{names of small RNAs to remove, default is NULL, can be set to character vector to specify; must match names in knownMiRNAs file.}

\item{bpparam}{TRUE or FALSE (default).}

\item{verbose}{print messages, TRUE or FALSE (default).}

\item{removedups}{remove multiple small RNAs mapping to the same read, TRUE (default) or FALSE. If TRUE, known miRNAs will be prioritized and known miRNA names must be in the format "miR-", "let-", "bantam-", "iab-".}

\item{overwrite}{overwrite existing output files, TRUE or FALSE (default).}
}
\value{
path to chimera output table.
}
\description{
Chimera process CLIP data
}
\examples{
testFastq <- system.file("extdata/SRR1742056.fastq.gz",package="CLIPflexR")
FqFile <- decompress(testFastq,overwrite=TRUE)
FaFile <- fastx_qtoa(FqFile)
FaFile_clip <- fastx_clipper(FaFile, writelog=FALSE)
myGenome <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
myIndex <- suppressWarnings(bowtie2_index(myGenome, overwrite = TRUE))
myBam <- suppressWarnings(bowtie_align(FaFile_clip,myIndex, overwrite=TRUE))
miRNAs <- system.file("extdata/hsa_mature.fa",package="CLIPflexR")
chimera_bed <- chimera_Process(myBam, knownMiRNAs= miRNAs, genomeIndex = myIndex, 
exclude="hsa-miR-19a-3p", overwrite=TRUE, fastas = NULL)
}
\author{
Kathryn Rozen-Gagnon
}
