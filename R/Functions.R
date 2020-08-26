file_path_sans_ext <- function (x, compression = FALSE) 
{
  if (compression) 
    x <- sub("[.](gz|bz2|xz)$", "", x)
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
}

file_ext <- function (x) 
{
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}


#' Convert Bam to Bed
#'
#' Convert Bam to Bed
#'
#'
#' @docType methods
#' @name bamtobed
#' @rdname bamtobed
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
bamtobed <- function(file,
                     outFile=gsub("\\.bam",".bed",file),
                     filtDup=FALSE){
  if(!file.exists(file))stop("No BAM file found at", file)
  if(!file.exists(outFile)){
  temp <- readGAlignments(file,
                          param=ScanBamParam(what = "qname"))
  names(temp) <- mcols(temp)$qname
  temp <- granges(temp)
  
  if(filtDup) temp <- temp[!duplicated(temp),]
  export.bed(temp,
             con = outFile)
}
  return(outFile)
}

#' Wrapper function for bzip2
#'
#' Wrapper function for bzip2
#'
#'
#' @docType methods
#' @name decompress
#' @rdname decompress
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToDecompress File to bzip.
#' @param outDir Output directory
#' @param keep keep (don't delete) input files.
#' @param overwrite overwrite existing output files.
#' @return Path to unzipped file.
#' @importFrom R.utils bunzip2 gunzip
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' decompress(testFQ,overwrite=TRUE)
#' @export
decompress <- function(fileToDecompress,
                       outDir=file_path_sans_ext(fileToDecompress),
                       keep=TRUE,
                       overwrite=FALSE){
  fileExt <- file_ext(fileToDecompress)
  cmd <- switch(fileExt,gz=gunzip,bz=bunzip2,bz2=bunzip2)
  
  if(!file.exists(fileToDecompress))stop("File does not exist")
  
  fileWithoutExtension <- file_path_sans_ext(fileToDecompress)
  if(file.exists(fileToDecompress) & !file.exists(fileWithoutExtension)){
    cmd(fileToDecompress,overwrite=overwrite,remove=!keep)
  }
  return(fileWithoutExtension)
}

#' Wrapper function for Compress
#'
#' Wrapper function for Compress
#'
#'
#' @docType methods
#' @name compress
#' @rdname compress
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToCompress File to bzip.
#' @param keep keep (don't delete) input files.
#' @param overwrite overwrite existing output files.
#' @param methodCompress Method to use for compression
#' @return Path to unzipped file.
#' @importFrom R.utils bunzip2 gunzip
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' decom <- decompress(testFQ,overwrite=TRUE)
#' com <- compress(decom,overwrite=TRUE)
#' @export
compress <- function(fileToCompress,
                       keep=TRUE,
                       methodCompress="gunzip",
                       overwrite=FALSE){
  fileExt <- switch(methodCompress,gunzip="gz",bunzip2="bz2")
  cmd <- switch(fileExt,gz=gunzip,bz=bunzip2,bz2=bunzip2)
  
  if(!file.exists(fileToCompress))stop("File does not exist")
  
  fileCompressed <- paste0(fileToCompress,".",fileExt)
  if(file.exists(fileToCompress) & !file.exists(fileCompressed)){
    cmd(fileToCompress,overwrite=overwrite,remove=!keep)
  }
  return(fileCompressed)
}

readinpeaks <- function(peaks,verbose=FALSE,sep="\t",header=FALSE){
if(verbose) message("Reading in peaks....",appendLF = FALSE)
if(is(peaks,"GRanges")){
  peaks <- peaks
  peaks$originalPeak <-  paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks))
}else if(is(peaks,"character")){
  if(!file.exists(peaks))stop(paste0("The file ",peaks," does not exist"))
  peaks <- read.table(peaks,sep=sep,header=header)
  if(!header){
    if(ncol(peaks) == 3)    colnames(peaks) <- c("seqnames","start","end")
    if(ncol(peaks) >= 6)    colnames(peaks)[1:6] <- c("seqnames","start","end","name","score","strand")
  }
  peaks <- makeGRangesFromDataFrame(peaks,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field="seqnames",
                                    start.field="start",
                                    end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)
  peaks$originalPeak <-  paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks))
  
}
  if(verbose) message("...done")
  if(verbose) message("Read in ",length(peaks)," peaks")
  
  return(peaks)
}


#' Retrieve sequences from under peaks/regions
#'
#' Retrieve sequences from under peaks/regions
#'
#'
#' @docType methods
#' @name fetchSequencesForClIP
#' @rdname fetchSequencesForClIP
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param peaks ClIP peaks to process.
#' @param fasta Genome fasta file.
#' @param resize Size of window in bp for resizing around peak center.
#' @param add5 Bp to add to 5' of resized peak.
#' @param add3 Bp to add to 3' of resized peak.
#' @param verbose Print messages.
#' @param bedHeader TRUE if peak file contains column headers.
#' @param bedSep Seperator in BED file.
#' @return Path to unzipped file.
#' @importFrom R.utils bunzip2 gunzip
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' decom <- decompress(testFQ,overwrite=TRUE)
#' com <- compress(decom,overwrite=TRUE)
#' @import Rsamtools
#' @export
fetchSequencesForClIP <- function(peaks,resize=NULL,fasta,add5=0,add3=0,verbose=FALSE,bedHeader=TRUE,bedSep="\t"){
  
  peaks <- readinpeaks(peaks,verbose=verbose,header = bedHeader,sep=bedSep)

  swd <- FaFile(fasta)
  if(!file.exists(index(swd)))indexFa(swd)
  fastaLens <- seqlengths(swd)
  start(peaks) <- start(peaks)-add5
  end(peaks) <- end(peaks)+add3
  sees <- seqlevels(peaks)
  seqlengths(peaks) <- fastaLens[match(sees,names(fastaLens),incomparables = 0)]
  
  
  
  Boundaries <- GRanges(seqlevels(peaks),IRanges(1,seqlengths(peaks)))
  
  if(!is.null(resize)){
    resizePeaks <- resize(peaks, resize, fix="center", use.names=TRUE)
    
  }else{
    resizePeaks <- peaks
  }
  names(resizePeaks) <- paste0(seqnames(resizePeaks),":",start(resizePeaks),"-",end(resizePeaks))
  #temp <- subsetByOverlaps(resizePeaks,Boundaries,type=c("within"))
  resizePeaks <- trim(resizePeaks)##added this line cause off error when extending going ut of bounds
  temp <- findOverlaps(resizePeaks,Boundaries,type=c("within"))
  temp2 <- Rsamtools::getSeq(FaFile(fasta),resizePeaks[temp@from])
  names(temp2) <- names(resizePeaks[temp@from])
  temp2
}


#' Convert Bam to Bed
#'
#' Convert Bam to Bed
#'
#'
#' @docType methods
#' @name annotatePeaksWithPatterns
#' @rdname annotatePeaksWithPatterns
#'
#' @author Kathryn Rozen-Gagnon
#' @param peaks ClIP peaks to process.
#' @param fasta Genome fasta file.
#' @param patterns Patterns to scan in ClIP peaks.
#' @param resize Size of window in bp for resizing around peak center.
#' @param add5 Bp to add to 5' of resized peak.
#' @param add3 Bp to add to 3' of resized peak.
#' @param verbose Print messages.
#' @param checkReverse Check the reverse strand for pattern.
#' @param bedHeader TRUE if peak file contains column headers.
#' @param bedSep Seperator in BED file.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="clipR")
#' @return Path 
#' @import GenomicAlignments GenomicFeatures TxDb.Mmusculus.UCSC.mm10.knownGene GenomicRanges
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet readDNAStringSet
#' @importMethodsFrom Biostrings as.data.frame coverage duplicated end getSeq match reverseComplement start strsplit substring vmatchPattern
#' @export
annotatePeaksWithPatterns  <- function(peaks,fasta,patterns,resize=64,add5=0,add3=0,verbose=FALSE,checkReverse=TRUE,bedHeader=TRUE,bedSep="\t"){
  
  peaks <- readinpeaks(peaks,verbose=verbose,header = bedHeader,sep=bedSep)
  
  
  if(!file.exists(fasta))stop(paste0("The file ",fasta," does not exist"))
  if(verbose) message("Indexing FASTA...",appendLF = FALSE)
  indexFa(fasta)
  if(verbose) message("...done")
  
  
  if(verbose) message("Aligning seqlevels across peaks and FASTA...",appendLF = FALSE)
  fastaLens <- FaFile(fasta) %>% seqlengths
  sees <- seqlevels(peaks)
  seqlengths(peaks) <- fastaLens[match(sees,names(fastaLens),incomparables = 0)]
  if(verbose) message("...done")
  
  if(verbose) message("Removing invalid peaks which are outside contig boundaries...",appendLF = FALSE)
  Boundaries <- GRanges(seqlevels(peaks),IRanges(1,seqlengths(peaks)))
  resizePeaks <- resize(peaks, resize, fix="center", use.names=TRUE)
  resizePeaks <- trim(resizePeaks)
  names(resizePeaks) <- paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks)) #rename accoridng to original (not extended) peakID to be able
  temp <- findOverlaps(resizePeaks,Boundaries,type=c("within"))
  validPeaks <- resizePeaks[temp@from]
  # validPeaks$originalPeak <- names(validPeaks)
  
  
  ##
  # Here we make a name for extended peaks...so we know its extended coordinates even if we change this GRanges!
  ##
  validPeaks$extendedPeak <- paste0(seqnames(validPeaks),":",start(validPeaks),"_",end(validPeaks), "_", strand(validPeaks))
  if(verbose) message("...done")
  if(verbose) message("Removed in ",length(validPeaks)-length(resizePeaks)," peaks")
  
  if(verbose) message("Retrieving sequence from under peaks....",appendLF = FALSE)
  myRes <- getSeq(Rsamtools::FaFile(fasta),validPeaks)
  
  ##
  ## We set myRes (the reverse complement) names by extended peak names.. 
  ##
  names(myRes) <- validPeaks$extendedPeak
  
  validPeaks$seqUnderExtendedPeak <- myRes
  if(verbose) message("...done")
  
  if(verbose) message("Retrieving patterns to search for....",appendLF = FALSE)
  if(class(patterns) == "DNAStringSet"){
    pattern <- patterns
  }else if(class(patterns) == "character"){
    motifSeq <- DNAStringSet(patterns)
    pattern <- as.character(motifSeq)
  }else{
    motifSeq <- readDNAStringSet(patterns)
    pattern <- as.character(motifSeq)
  }
  if(verbose) message("...done")
  if(verbose) message("Read in ",length(pattern)," patterns")
  
  
  # rm(fixedInRegion)
  if(verbose) message("Searching for ",length(pattern)," patterns in peaks....")
  mergePeaksAndSites <- patternCallList <- list()
  for(i in 1:length(pattern)){
    if(verbose) message("Search for ",pattern[i])
    peaks_Sites <- validPeaks
    fixedInRegion <- vmatchPattern(pattern=pattern[i],subject=myRes) %>% unlist()
    
    
    startOfPmatch <- resize(fixedInRegion,width=1,fix="start") %>% unname %>% as.character
    
    ##
    ## Pos is made from names from original sequence we scanned in..which came from extened peak coordnates. 
    ##  or..Names contains extended position and strand information because myRes did in its name.
    ##
    pos <- matrix(unlist(strsplit(gsub(".*:","",names(fixedInRegion)),"_")),ncol=3,byrow =T)
    seq <- gsub(":.*","",names(fixedInRegion))
    
    ###
    ## Here we lose strand information from our GRanges!! 
    ## sitesGRpos are from Pos strand based on name... but they dont have strand information in them !! 
    ## That why we need to reverse complement.
    ## We could add back strand information somewhere here but doesnt matter as we take care of it ourselves
    
    ## I checked a few sites in IGV from original FASTA and all good! And now we can remember why too :)
    ###
    sitesGRpos <- GRanges(seq[pos[,3] =="+"],
                          IRanges(as.numeric(pos[pos[,3] =="+",1])+start(fixedInRegion[pos[,3] =="+"])-1,
                                  as.numeric(pos[pos[,3] =="+",1])+end(fixedInRegion[pos[,3] =="+"])-1),startOfPmatch=startOfPmatch[pos[,3] =="+"])
    sitesGRpos$PeakID <- names(fixedInRegion[pos[,3] =="+"])
    
    sitesGRneg <- GRanges(seq[pos[,3] =="-"],
                          IRanges(as.numeric(pos[pos[,3] =="-",2])-end(fixedInRegion[pos[,3] =="-"])+1,
                                  as.numeric(pos[pos[,3] =="-",2])-start(fixedInRegion[pos[,3] =="-"])+1),startOfPmatch=startOfPmatch[pos[,3] =="-"])
    sitesGRneg$PeakID <- names(fixedInRegion[pos[,3] =="-"])
    
    ## So we need to reverse complement based on
    sitesFA <- c(fetchSequencesForClIP(sitesGRpos,resize=NULL,fasta),reverseComplement(fetchSequencesForClIP(sitesGRneg,resize=NULL,fasta)))
    sitesFAExtended <- c(fetchSequencesForClIP(sitesGRpos,resize=NULL,fasta,add5=add5,add3=add3),
                         reverseComplement(fetchSequencesForClIP(sitesGRneg,resize=NULL,fasta,add5=add5,add3=add3))
    )
    sitesGR <- suppressWarnings(c(sitesGRpos,sitesGRneg))
    
    # getSeq(Rsamtools:::FaFile(fasta),sitesGR)
    sitesGR$Seq <- sitesFA
    sitesGR$SeqExtended <- sitesFAExtended
    sitesGR$siteOfPattern <- granges(sitesGR) %>% unname %>% as.character
    sitesGR$centerOfPattern <-  resize(sitesGR,width=1,fix="center") %>% unname %>% as.character
    dfSitesGR <- as.data.frame(sitesGR)
    # peaks_Sites$extendedRanges <- resize(granges(peaks_Sites), reSize, fix="center", use.names=TRUE)
    # peaks_Sites$peakNames <- paste0(seqnames(peaks_Sites$extendedRanges),":",start(peaks_Sites$extendedRanges),"-",end(peaks_Sites$extendedRanges))
    peaks_SitesDF <- as.data.frame(peaks_Sites,row.names = NULL)
    
    mergePeaksAndSites[[i]] <- merge(peaks_SitesDF,dfSitesGR,by.x="extendedPeak",by.y="PeakID",all=FALSE)
    patternCallList[[i]] <- mergePeaksAndSites[[i]] %>% dplyr::select(Seq,SeqExtended,startOfPmatch,centerOfPattern,siteOfPattern,originalPeak,extendedPeak)
    rm(peaks_SitesDF)
    rm(peaks_Sites)
  }
  #if(verbose) message("....finished search")
  names(mergePeaksAndSites) <- names(patternCallList) <- as.character(pattern)
  return(list(mergePeaksAndSites,patternCallList))
}
