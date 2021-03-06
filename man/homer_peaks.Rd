% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapperFunctions_homer.R
\docType{methods}
\name{homer_peaks}
\alias{homer_peaks}
\title{Wrapper function for Homer's makeTagDirectory and findPeaks}
\usage{
homer_peaks(
  fileTofqs,
  maketagdir = "makeTagDirectory",
  findpeaks = "findPeaks",
  format = "bed",
  createSingleTagsTSV = TRUE,
  tagdir = file.path(dirname(fileTofqs), gsub("\\\\.bed", "",
    make.names(basename(fileTofqs)))),
  style = "factor",
  foldEnrichmentOverLocal = 2,
  localSize = 10000,
  strand = "seperate",
  minDist = 50,
  size = 10,
  fragLength = 10,
  genomeSize = NULL,
  stderr = file.path(dirname(fileTofqs), paste0(basename(fileTofqs),
    "_homer_stderr.txt")),
  stdout = file.path(dirname(fileTofqs), paste0(basename(fileTofqs),
    "_homer_stderr.txt")),
  useClipRConda = ifelse(is.null(getOption("CLIPflexR.condaEnv")), FALSE, TRUE),
  additional_Args = NULL,
  verbose = FALSE,
  writelog = T
)
}
\arguments{
\item{fileTofqs}{path to file to process (BED).}

\item{maketagdir}{path to makeTagDirectory from Homer toolkit.}

\item{findpeaks}{path to findpeaks from Homer toolkit.}

\item{format}{input format, "bed" (default).}

\item{createSingleTagsTSV}{create a single tags.tsv file for all "chromosomes".}

\item{tagdir}{name of Tag Directory.}

\item{style}{type of Homer peak calling, "factor", "histone", "groseq", "tss", "dnase", "super" or "mCs".}

\item{foldEnrichmentOverLocal}{fold enrichment over local tag count.}

\item{localSize}{region to check for local tag enrichment.}

\item{strand}{find peaks using tags on "both" or "separate" (default) strands.}

\item{minDist}{minimum distance between peaks.}

\item{size}{peak size.}

\item{fragLength}{approximate fragment length.}

\item{genomeSize}{genome size, default is NULL (2E9, applicable for human or mouse), set integer to specify for your genome.}

\item{stderr}{path to stderr file.}

\item{stdout}{path to stdout file.}

\item{useClipRConda}{use conda environment installed by Herper, TRUE (default) or FALSE.}

\item{additional_Args}{additional arguments to be passed to system call.}

\item{verbose}{print messages, TRUE or FALSE (default).}

\item{writelog}{write stderr/stdout logs, TRUE (default) or FALSE.}
}
\value{
path to unzipped file
}
\description{
Wrapper function for Homer's makeTagDirectory and findPeaks
}
\examples{
\dontrun{
testFasta <- system.file("extdata/hg19Small.fa",package = "CLIPflexR")
myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package = "CLIPflexR")
FqFile <- decompress(testFQ,overwrite = TRUE)
FqFile_QF <- fastq_quality_filter(FqFile)
FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
bam <- suppressWarnings(bowtie_align(FqFile_QFColStpClipped,myIndex, overwrite = TRUE))
bed <- bamtobed(bam)
homer_peaks(bed)
}
}
\author{
Kathryn Rozen-Gagnon
}
