if(getRversion() >= "2.15.1")  utils::globalVariables(c("seq","Seq","SeqExtended","centerOfPattern",
                                                        "siteOfPattern","originalPeak","extendedPeak",
                                                        "baseDir"))


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

readinpeaks <- function(peaks,verbose=FALSE,sep="\t",header=FALSE){
  if(verbose) message("Reading in peaks....",appendLF = FALSE)
  if(is(peaks,"GRanges")){
    peaks <- peaks
    peaks$originalPeak <-  paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks))}
  else if(is(peaks,"character")){
    if(!file.exists(peaks))stop(paste0("The file ",peaks," does not exist"))
    peaks <- read.delim(peaks,sep=sep,header=header)
    if(!header){
      if(ncol(peaks) == 3)    colnames(peaks) <- c("seqnames","start","end")  ##not sure about this case
      if(ncol(peaks) >= 6)    colnames(peaks)[1:6] <- c("seqnames","start","end","name", "score","strand")}
    peaks <- makeGRangesFromDataFrame(peaks,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="seqnames",
                                      start.field="start",
                                      end.field="end",
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE)
    peaks$originalPeak <-  paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks))}
  
   else if(is.data.frame(peaks)) { 
      peaks <- peaks 
    
    peaks <- makeGRangesFromDataFrame(peaks,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="seqnames",
                                      start.field="start",
                                      end.field="end",
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE)
    peaks$originalPeak <-  paste0(seqnames(peaks),":",start(peaks),"_",end(peaks), "_", strand(peaks))}

  if(verbose) message("...done")
  if(verbose) message("Read in ",length(peaks)," peaks")
  
  return(peaks)
}

utils::globalVariables(c("V4"))

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
#' @param file path to file to process (BAM).
#' @param outFile path to output file (BED).
#' @param filtDup output index name.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta, overwrite = TRUE)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
#' FqFile_QFColStpClipped <- fastx_clipper(FqFile_QFColStripped)
#' bam <- bowtie_align(FqFile_QFColStpClipped,myIndex, overwrite = TRUE)
#' bamtobed(bam)
#' @return path to BED. 
#' @import GenomicAlignments
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @importFrom methods is
#' @importFrom stats setNames
#' @export
bamtobed <- function(file,
                     outFile=gsub("\\.bam",".bed",file),
                     filtDup=FALSE){
  if(!file.exists(file))stop("No BAM file found at", file)
  if(file.exists(outFile))stop("Output BED alredy exists at", outFile)
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

#' Decompress
#'
#' Wrapper function for gunzip, bunzip
#'
#'
#' @docType methods
#' @name decompress
#' @rdname decompress
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToDecompress path to file to decompress.
#' @param outDir output directory.
#' @param keep keep (don't delete) input files, TRUE (default) or FALSE.
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).
#' @return path to unzipped file.
#' @importFrom R.utils bunzip2 gunzip
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
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
  if(file.exists(fileToDecompress) & (overwrite | !file.exists(fileWithoutExtension))){
    cmd(fileToDecompress,overwrite=overwrite,remove=!keep)
  }
  return(fileWithoutExtension)
}

#' Compress
#'
#' Wrapper function for gzip, bzip
#'
#'
#' @docType methods
#' @name compress
#' @rdname compress
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToCompress path to file to compress.
#' @param keep keep (don't delete) input files, TRUE (default) or FALSE.
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).
#' @param methodCompress method to use for compression.
#' @return path to zipped file.
#' @importFrom R.utils bunzip2 gunzip gzip bzip2
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' decom <- decompress(testFQ,overwrite=TRUE)
#' com <- compress(decom,overwrite=TRUE)
#' @export
compress <- function(fileToCompress,
                     keep=TRUE,
                     methodCompress="gzip",
                     overwrite=FALSE){
  fileExt <- switch(methodCompress,gzip="gz",bzip2="bz2")
  cmd <- switch(fileExt,gz=gzip,bz=bzip2,bz2=bzip2)
  
  if(!file.exists(fileToCompress))stop("File does not exist")
  
  fileCompressed <- paste0(fileToCompress,".",fileExt)
  if(file.exists(fileToCompress) & (overwrite |!file.exists(fileCompressed))){
    cmd(filename  = fileToCompress, destname= fileCompressed,overwrite=overwrite,remove=!keep)
  }
  return(fileCompressed)
}



#' Retrieve sequences from under peaks/regions
#'
#' Retrieve sequences from under peaks/regions
#'
#'
#' @docType methods
#' @name fetchSequencesForCLIP
#' @rdname fetchSequencesForCLIP
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param peaks CLIP peaks to process; accepts path to peak files (BED file or BED-formatted tab-delimited text files) or R objects (GRanges or BED-formatted data.frame); peak locations must be unique.
#' @param fasta path to genome file (fasta).
#' @param resize size of window in bp for resizing around peak center, default is NULL (set by specifying integer).
#' @param add5 bp to add to 5' of resized peak, default is 0 (set by specifying integer).
#' @param add3 bp to add to 3' of resized peak, default is 0 (set by specifying integer).
#' @param verbose print messages, TRUE or FALSE (default).
#' @param bedHeader if peak file contains column headers, TRUE (default) or FALSE; if TRUE, column names must be "seqnames" (chromosome), "start" (peak start), "end" (peak end), "strand" (peak strand).
#' @param bedSep separator in BED file.
#' @return path to unzipped file.
#' @importFrom R.utils bunzip2 gunzip
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' decom <- decompress(testFQ,overwrite=TRUE)
#' com <- compress(decom,overwrite=TRUE)
#' @import Rsamtools
#' @export
fetchSequencesForCLIP <- function(peaks,resize=NULL,fasta,add5=0,add3=0,verbose=FALSE,bedHeader=TRUE,bedSep="\t"){
  
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


#' annotatePeaksWithPatterns
#'
#' Search for small RNA targets, or any pattern, in peaks
#'
#'
#' @docType methods
#' @name annotatePeaksWithPatterns
#' @rdname annotatePeaksWithPatterns
#'
#' @author Kathryn Rozen-Gagnon
#' @param peaks CLIP peaks to process; accepts path to peak files (BED file or BED-formatted tab-delimited text files) or R objects (GRanges or BED-formatted data.frame); peak locations must be unique.
#' @param fasta path to genome file (fasta).
#' @param patterns patterns to scan in CLIP peaks (character vector, DNAStringSet, or file path to a fasta sequence).
#' @param resize size of window in bp for resizing around peak center.
#' @param add5 bp to add to 5' of resized peak, default is 0 (set by specifying integer).
#' @param add3 bp to add to 3' of resized peak, default is 0 (set by specifying integer).
#' @param verbose print messages, TRUE or FALSE (default).
#' @param checkReverse check the reverse strand for pattern.
#' @param bedHeader if peak file contains column headers, TRUE (default) or FALSE; if TRUE, column names must be "seqnames" (chromosome), "start" (peak start), "end" (peak end), "strand" (peak strand).
#' @param bedSep separator in BED file.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
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
  fastaLens <- seqlengths(FaFile(fasta))
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
    pattern <- as.character(patterns)
  }else if(grepl(".fa", patterns)){
    motifSeq <- readDNAStringSet(patterns)
    pattern <- as.character(motifSeq)
  } else if (class(patterns) == "character"){
    pattern <- patterns }
  if(verbose) message("...done")
  if(verbose) message("Read in ",length(pattern)," patterns")
  

  
  
  # rm(fixedInRegion)
  if(verbose) message("Searching for ",length(pattern)," patterns in peaks....")
  mergePeaksAndSites <- patternCallList <- list()
  for(i in 1:length(pattern)){
    if(verbose) message("Search for ",pattern[i])
    peaks_Sites <- validPeaks
    fixedInRegion <- vmatchPattern(pattern=pattern[i],subject=myRes) %>% unlist()
    if(!isEmpty(fixedInRegion)) {
      
      startOfPmatch <- resize(fixedInRegion,width=1,fix="start") %>% unname %>% as.character
      pos <- matrix(unlist(strsplit(gsub(".*:","",names(fixedInRegion)),"_")),ncol=3,byrow =T)
      seq <- gsub(":.*","",names(fixedInRegion))
      
      sitesGRpos <- GRanges(seq[pos[,3] =="+"],
                            IRanges(as.numeric(pos[pos[,3] =="+",1])+start(fixedInRegion[pos[,3] =="+"])-1,
                                    as.numeric(pos[pos[,3] =="+",1])+end(fixedInRegion[pos[,3] =="+"])-1),startOfPmatch=startOfPmatch[pos[,3] =="+"])
      sitesGRpos$PeakID <- names(fixedInRegion[pos[,3] =="+"])
      
      sitesGRneg <- GRanges(seq[pos[,3] =="-"],
                            IRanges(as.numeric(pos[pos[,3] =="-",2])-end(fixedInRegion[pos[,3] =="-"])+1,
                                    as.numeric(pos[pos[,3] =="-",2])-start(fixedInRegion[pos[,3] =="-"])+1),startOfPmatch=startOfPmatch[pos[,3] =="-"])
      sitesGRneg$PeakID <- names(fixedInRegion[pos[,3] =="-"])
      if(isEmpty(sitesGRpos)) { sitesFA <-  Biostrings::reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta))
      sitesFAExtended <- Biostrings::reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta,add5=add5,add3=add3))
      sitesGR <- sitesGRneg
      
      sitesGR$Seq <- sitesFA
      sitesGR$SeqExtended <- sitesFAExtended
      sitesGR$siteOfPattern <- granges(sitesGR) %>% unname %>% as.character
      sitesGR$centerOfPattern <-  GenomicRanges::resize(sitesGR,width=1,fix="center") %>% unname %>% as.character
      dfSitesGR <- as.data.frame(sitesGR)
      peaks_SitesDF <- as.data.frame(peaks_Sites)
      
      mergePeaksAndSites[[i]] <- merge(peaks_SitesDF,dfSitesGR,by.x="extendedPeak",by.y="PeakID",all=FALSE)
      patternCallList[[i]] <- mergePeaksAndSites[[i]] %>% dplyr::select(Seq,SeqExtended,startOfPmatch,centerOfPattern,siteOfPattern,originalPeak,extendedPeak)
      rm(peaks_SitesDF)
      rm(peaks_Sites)
      } 
      else if(isEmpty(sitesGRneg)) {
        sitesFA <- fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta)
        sitesFAExtended <- fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta,add5=add5,add3=add3)
        sitesGR <- sitesGRpos
        
        sitesGR$Seq <- sitesFA
        sitesGR$SeqExtended <- sitesFAExtended
        sitesGR$siteOfPattern <- granges(sitesGR) %>% unname %>% as.character
        sitesGR$centerOfPattern <-  GenomicRanges::resize(sitesGR,width=1,fix="center") %>% unname %>% as.character
        dfSitesGR <- as.data.frame(sitesGR)
        peaks_SitesDF <- as.data.frame(peaks_Sites)
        
        mergePeaksAndSites[[i]] <- merge(peaks_SitesDF,dfSitesGR,by.x="extendedPeak",by.y="PeakID",all=FALSE)
        patternCallList[[i]] <- mergePeaksAndSites[[i]] %>% dplyr::select(Seq,SeqExtended,startOfPmatch,centerOfPattern,siteOfPattern,originalPeak,extendedPeak)
        rm(peaks_SitesDF)
        rm(peaks_Sites)
      } else {
        sitesFA <- c(fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta),reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta)))
        sitesFAExtended <- c(fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta,add5=add5,add3=add3),
                             reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta,add5=add5,add3=add3))
        )
        sitesGR <- suppressWarnings(c(sitesGRpos,sitesGRneg))
        
        sitesGR$Seq <- sitesFA
        sitesGR$SeqExtended <- sitesFAExtended
        sitesGR$siteOfPattern <- granges(sitesGR) %>% unname %>% as.character
        sitesGR$centerOfPattern <-  GenomicRanges::resize(sitesGR,width=1,fix="center") %>% unname %>% as.character
        dfSitesGR <- as.data.frame(sitesGR)
        peaks_SitesDF <- as.data.frame(peaks_Sites)
        
        mergePeaksAndSites[[i]] <- merge(peaks_SitesDF,dfSitesGR,by.x="extendedPeak",by.y="PeakID",all=FALSE)
        patternCallList[[i]] <- mergePeaksAndSites[[i]] %>% dplyr::select(Seq,SeqExtended,startOfPmatch,centerOfPattern,siteOfPattern,originalPeak,extendedPeak)
        rm(peaks_SitesDF)
        rm(peaks_Sites)}} 
    else { message(pattern[i], "not found")}}
  if(verbose) message("....finished search")
  names(mergePeaksAndSites) <- names(patternCallList) <- as.character(pattern)
  return(list(mergePeaksAndSites,patternCallList))
  
}


#' bowtie2_index
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
#' @param genomeFasta path to file to process.
#' @param outIndex path to output index.
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).
#' @param threads number of threads to use, default is 1.
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' bowtie2_index(testFasta, overwrite = TRUE)
#' @importFrom Rbowtie2 bowtie2_build
#' @return Path to index 
#' @export
bowtie2_index <- function(genomeFasta,
                          outIndex=gsub("\\.fa","",genomeFasta),
                          overwrite=F,threads=1) {
  bowtieArgs <- paste0("--threads ",threads)
    suppressMessages(
      bowtie2_build(references=genomeFasta,
                    bt2Index=outIndex,
                    overwrite = overwrite,bowtieArgs))
  return(outIndex)
}


#' bowtie_align
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
#' @param fq path to file to process (fasta or fastq).
#' @param index bowtie2 index name without extension.
#' @param outbam output bam name, not required if mode is set. If mode is set to "NULL", specify "fq" (default - name of file to process/reads) or "index" (name of index)
#' @param keepSAM keep bowtie2 SAM output file, TRUE (default) or FALSE.
#' @param inputFormat Format of reads, "fastq" or "fasta" (default).
#' @param mode mapping mode, "genome_map" (default - 1 mismatch in seed alignment, seed substring of 18) or "reverse_map" (for mapping short sequences to reads, i.e. miRNAs or chimera analysis, 0 mismatches and a seed substring of 18); set mode=NULL to custom set bowtie2 options using arguments below or additional arguments.
#' @param maxMismatches max mismatches in seed alignment, (default is 0).
#' @param seedSubString length of seed substrings (default is 18).
#' @param threads number of threads to use (default is 1).
#' @param report_k number of alignments to report, default is NULL (the best alignment is reported), set to integer to specify.
#' @param keep_all report all alignments, default is NULL (the best alignment is reported), set to TRUE to report all.
#' @param soft_clip allow soft clipping, default is NULL (no soft clipping), set to TRUE to to soft clip.
#' @param additional_Args any additional mapping arguments not set above, default is NULL, can be set by specifying a single character string with spaces between arguments (see example below); run "Rbowtie2::bowtie2_usage()" to see all options; please note that due to a Rbowtie2 bug, the "--no-1mm-upfront" is not available.
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta, overwrite=TRUE)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_FF <- ctk_fastqFilter(FqFile,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile_FF,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_ColStrip <- ctk_stripBarcode(FqFile_Col,linkerlength=5, inputFormat="fastq")
#' 
#' ##map reads to genome
#' bowtie_align(FqFile_ColStrip,myIndex, inputFormat="fastq", overwrite=TRUE)
#' 
#' ##map reads to genome, custom mode
#' bowtie_align(FqFile_ColStrip,myIndex, mode=NULL, 
#' soft_clip=TRUE, additional_Args="--mp 15 -R 4", inputFormat="fastq", overwrite=TRUE) 
#' 
#' ##map miRNAs to reads
#' miRNAs <- system.file("extdata/hsa_mature.fa",package="CLIPflexR")
#' FaFile_Fa <- fastx_qtoa(FqFile_ColStrip)
#' readIndex <- bowtie2_index(FaFile_Fa, overwrite=TRUE)
#' bowtie_align(miRNAs,readIndex, mode="reverse_map",overwrite=TRUE)

#' @importFrom Rsamtools asBam indexBam sortBam
#' @importFrom Rbowtie2 bowtie2
#' @return path to sorted BAM. 
#' @export
bowtie_align <- function(fq,index,
                          outbam="fq",
                         inputFormat="fasta", keepSAM = T, mode = "genome_map", maxMismatches=0,seedSubString=18,threads=1,report_k=NULL, keep_all=NULL, soft_clip = NULL, additional_Args = NULL, overwrite=FALSE) {
  
  if(inputFormat == "fasta"){
    optionFormat <- "-f"
  }else{
    optionFormat <- ""
  }
  
  if(file_ext(fq) == "gz"){
    R.utils::gunzip(fq,
                    destname=gsub("\\.gz$","",fq))
    fq <- gsub("\\.gz$","",fq)
    
  }
  if(!is.null(mode)){
    if(mode=="genome_map") {
      bowtieArgs <- paste0(optionFormat, " -N 1 -L 18 --threads ",threads)
      bam=file.path(dirname(fq),
                    paste0(gsub("\\.fa|\\.fasta|\\.fastq|\\.fq","",basename(fq)),".bam"))
        suppressMessages(bowtie2(bt2Index = index,
                                 samOutput = gsub("\\.bam$",".sam",
                                                  bam),
                                 seq1 = fq,
                                 bowtieArgs, overwrite = overwrite)) }
    else if(mode=="reverse_map") {
      bowtieArgs <- paste0(optionFormat,  " --score-min C,0,0 -L 18 -a --threads ",threads)
      bam=file.path(dirname(index),
                    paste0(gsub("\\.fa|\\.fasta|\\.fastq|\\.fq","",basename(index)),".bam"))
        suppressMessages(bowtie2(bt2Index = index,
                                 samOutput = gsub("\\.bam$",".sam",
                                                  bam),
                                 seq1 = fq,
                                 bowtieArgs, overwrite = overwrite)) }}
  else if(is.null(mode)){
    bowtieArgs <- paste0(optionFormat, 
                         paste0(" --threads ", threads),
                         ifelse(!is.null(maxMismatches), paste0(" -N ", maxMismatches), ""),
                         ifelse(!is.null(seedSubString), paste0(" -L ", seedSubString), ""),
                         ifelse(!is.null(report_k), paste0(" -k ", as.integer(report_k)), ""),
                         ifelse(!is.null(keep_all), " -a ", ""),
                         ifelse(!is.null(soft_clip), " --local ", ""),
                         ifelse(!is.null(additional_Args), paste0(" ", additional_Args), "")
    )
      if(outbam=="fq") {
        bam=file.path(dirname(fq),
                      paste0(gsub("\\.fa|\\.fasta|\\.fastq|\\.fq","",basename(fq)),".bam"))
        suppressMessages(bowtie2(bt2Index = index,
                                 samOutput = gsub("\\.bam$",".sam",
                                                  bam),
                                 seq1 = fq,
                                 bowtieArgs, overwrite = overwrite)) } else if(outbam=="index"){
                                   bam=file.path(dirname(index), paste0(gsub("\\.fa|\\.fasta|\\.fastq|\\.fq","",basename(index)),".bam"))
                                   suppressMessages(bowtie2(bt2Index = index,
                                                            samOutput = gsub("\\.bam$",".sam",
                                                                             bam),
                                                            seq1 = fq,
                                                            bowtieArgs, overwrite = overwrite))  } 
    }
  
  asBam(gsub("\\.bam$",".sam",
             bam),
        gsub("\\.bam$",".temp",
             bam), overwrite = overwrite)
  sortBam(gsub("\\.bam$",".temp.bam",
               bam),
          destination = gsub("\\.bam$","",
                             bam))
  unlink(gsub("\\.bam$",".temp.bam",
              bam)) # remove temp.bam
  unlink(gsub("\\.bam$",".temp.bam.bai",
              bam))
  indexBam(bam)
  if(!keepSAM) {
    unlink(gsub("\\.bam$",".sam",
                bam))
  }
  
  return(bam)
  
}

#' countFromBed to make count matrix
#'
#' count from mapped beds to make count matrix
#'
#'
#' @docType methods
#' @name countFromBed
#' @rdname countFromBed
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param Bed BED bed files containing mapped reads.
#' @param GR genomic ranges object over which to count. 
#' @param notStranded if strand information should be considered when counting, TRUE (ignore strand) or FALSE (consider strand; default).
#' @param interFeature discard reads mapping to multiple features, TRUE or FALSE (default).
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' bowtie2_index(testFasta, overwrite = TRUE)
#' @import GenomicRanges GenomicAlignments
#' @importMethodsFrom SummarizedExperiment assay
#' @importMethodsFrom rtracklayer import.bed
#' @importMethodsFrom GenomeInfoDb seqlengths "seqlengths<-" seqlevels
#' @return count matrix.
#' @export
countFromBed <- function(Bed,GR,notStranded=TRUE,interFeature=FALSE){
  reads <- import.bed(Bed)
  fk <- summarizeOverlaps(GR,reads,ignore.strand = notStranded,inter.feature=interFeature)
  assay(fk)
}


#' Build bigwigs
#'
#' @rdname CLIP_bw2
#' @param sort_bam path to sorted BAM file.
#' @param res_dir result directory, default is directory where BAM file is located, to specify other/create new directory, enter path as character string.
#' @param normalized normalized to total reads.
#' @param stranded Whether to make separate bigwigs for each strand, TRUE or FALSE (default).
#' @param verbose print messages, TRUE or FALSE (default).
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta, overwrite=TRUE)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite = TRUE)
#' FqFile_FF <- ctk_fastqFilter(FqFile,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile_FF,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose = TRUE)
#' FqFile_ColStrip <- ctk_stripBarcode(FqFile_Col,linkerlength=5, inputFormat="fastq") 
#' ##map reads to genome
#' mapped <- bowtie_align(FqFile_ColStrip,myIndex, mode="genome_map", inputFormat="fastq")
#' wig <- CLIP_bw2(mapped)
#' @docType methods
#' @import  Rsamtools
#' @return bigwig.
#' @export
CLIP_bw2 <- function(sort_bam,res_dir=dirname(sort_bam),normalized=TRUE, stranded=FALSE,verbose=FALSE){
  if(!dir.exists(res_dir)){
    dir.create(res_dir) 
    if(verbose){
    message("creating output directory ", res_dir)}
  }else if(dir.exists(res_dir) & verbose){
      message("bigwig will be written to existing directory: ",res_dir)}
  samID <- gsub("_sort.bam","",basename(sort_bam)) 
  bw_out <- file.path(res_dir,paste0(samID,".bigwig"))
  if(!stranded){
    if(normalized){
      allChromosomeStats <- Rsamtools::idxstatsBam(sort_bam)
      mapped <- sum(allChromosomeStats[,"mapped"])
      toRun <- coverage(Rsamtools::BamFile(sort_bam),weight = (10^6)/mapped)
      export.bw(toRun,con=bw_out)
      }
  #     else if(notStranded==FALSE){
  #     mapped <- scanBam(BamFile(sort_bam), flag = scanBamFlag(isUnmappedQuery = FALSE)))
  #     mapped_pos <- length(which(mapped$strand=="+"))
  #     toRun_pos <- Rsamtools::BamFile(sort_bam),weight = (10^6)/mapped)
  # scanBam(BamFile(sort_bam),param=ScanBamParam(what = c("qname","seq"),flag = scanBamFlag(isUnmappedQuery = FALSE)))
  #     mapped_neg <- length(which(mapped$strand=="-"))
  #     }
    else{
      toRun <- coverage(BamFile(sort_bam))
      export.bw(toRun,con=bw_out)
      }
  }else{
    if(normalized){
      allChromosomeStats <- Rsamtools::idxstatsBam(sort_bam)
      mapped <- sum(allChromosomeStats[,"mapped"])
      myparam <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))
      toRun <- coverage(Rsamtools::BamFile(sort_bam),param=myparam,weight = (10^6)/mapped)
      export.bw(toRun,con=gsub("\\.bigwig","_Pos.bigwig",bw_out))
      myparam <- ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE))
      toRun <- coverage(Rsamtools::BamFile(sort_bam),param=myparam,weight = (10^6)/mapped)
      export.bw(toRun,con=gsub("\\.bigwig","_Neg.bigwig",bw_out))
    }
    #     else if(notStranded==FALSE){
    #     mapped <- scanBam(BamFile(sort_bam), flag = scanBamFlag(isUnmappedQuery = FALSE)))
    #     mapped_pos <- length(which(mapped$strand=="+"))
    #     toRun_pos <- Rsamtools::BamFile(sort_bam),weight = (10^6)/mapped)
    # scanBam(BamFile(sort_bam),param=ScanBamParam(what = c("qname","seq"),flag = scanBamFlag(isUnmappedQuery = FALSE)))
    #     mapped_neg <- length(which(mapped$strand=="-"))
    #     }
    else{
      myparam <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))
      toRun <- coverage(Rsamtools::BamFile(sort_bam),param=myparam)
      export.bw(toRun,con=gsub("\\.bigwig","_Pos.bigwig",bw_out))
      myparam <- ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE))
      toRun <- coverage(Rsamtools::BamFile(sort_bam),param=myparam)
      export.bw(toRun,con=gsub("\\.bigwig","_Neg.bigwig",bw_out))
    }
  }
  
}

#' Build index for reference genome
#'
#' @rdname CLIP_buildIDX
#' @param ref_seq path to reference genome.
#' @param res_dir result directory.
#' @param preFIX file name prefix for indexed genome.
#' @param aligner Assign aligner bowtie2/subread/bwa/hisat2.
#' @param thread threads used in processing, default is 4.
#' @param bwa executable bwa file path.
#' @import Rbowtie2 Rsubread Rhisat2
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' @docType methods
#' @return indexed reference genome.
#' @export
CLIP_buildIDX <- function(ref_seq=NULL,res_dir=NULL,preFIX=NULL,aligner=NULL,thread=4,bwa=NULL){
  if(aligner=="bowtie2"){
    bowtie2_build(references = ref_seq,bt2Index=file.path(res_dir,paste0(preFIX,"_bowtie2")),
                  overwrite = TRUE,"--threads 4 --quiet")
  }else if(aligner=="subread"){
    buildindex(basename = file.path(res_dir,paste0(preFIX,"_subread")),reference = ref_seq)
  }else if(aligner=="bwa"){
    system(paste(bwa,"index","-p",paste0(preFIX,"_bwa"),ref_seq,sep=" "))
  }else if(aligner=="hisat2"){
    # system(paste(hisat_build,"-p",thread,ref_seq,file.path(res_dir,paste0(preFIX,"_hisat2")),sep=" "))
    hisat2_build(references=ref_seq,outdir=res_dir, prefix=paste0(preFIX,"_hisat2"),p=thread,
                 force=TRUE, strict=TRUE, execute=TRUE)
  }else{"Please assign an available aligner: bowtie2, subread, bwa and hisat2"}}

#' Align to Genome
#'
#' @rdname CLIP_align
#' @param samID sample id in the sample sheet.
#' @param samSheet sample sheet.
#' @param res_dir result directory.
#' @param genome_idx path to indexed genome.
#' @param aligner Assign aligner bowtie2/bwa_mem/bwa_alb/subread/subjunc/hisat2
#' @param thread thread used in the processing, default is 4
#' @param bwa executable bwa file path. 
#' @import Rbowtie2 Rsamtools Rsubread ShortRead 
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' @docType methods
#' @return path to sorted BAM file.
#' @export
CLIP_align <- function(samID,res_dir=NULL,genome_idx=NULL,aligner=NULL,samSheet=NULL,thread=4,bwa=NULL){
  # samID <- gsub("_trimmed.fq","",basename(fastq))
  fastq <- samSheet$FQLocation[samSheet$GenomicsID==samID] 
  # samID <- gsub(".fq","",basename(fastq)) # automatic get sample ID from fastq name
  if(aligner=="bowtie2"){
    ci <- unlist(strsplit(basename(genome_idx),split="_"))
    if(ci[length(ci)]=="bowtie2"){
      file.copy(fastq,res_dir,recursive = TRUE)
      fq <- file.path(res_dir,basename(fastq))
      system(paste("gunzip",fq,sep=" "))
      fq <- file.path(res_dir,gsub(".gz","",basename(fq)))
      aln_sam <- file.path(res_dir,paste0(samID,"_BW2_aln.sam"))
      aln_bam <- file.path(res_dir,paste0(samID,"_BW2_aln"))
      sort_bam <- file.path(res_dir,paste0(samID,"_BW2_sort"))
      bowtie2(bt2Index = genome_idx,samOutput=aln_sam,
              seq1=fq,seq2=NULL,overwrite=TRUE,paste0("--threads ",thread))
      asBam(aln_sam,aln_bam)
      sortBam(paste0(aln_bam,".bam"),sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      # unlink(aln_sam)
      unlink(fq)
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by bowtie2.")}
  }else if(aligner=="subread"){
    ci <- unlist(strsplit(genome_idx,split="_"))
    if(ci[length(ci)]=="subread"){
      aln_bam <- file.path(res_dir,paste0(samID,"_SR_aln.bam"))
      sort_bam <- file.path(res_dir,paste0(samID,"_SR_sort"))
      align(index=genome_idx,readfile1=fastq,
            output_file=aln_bam,unique=FALSE,type="dna")
      sortBam(aln_bam,sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by subread.")}
  }else if(aligner=="subjunc"){
    ci <- unlist(strsplit(genome_idx,split="_"))
    if(ci[length(ci)]=="subread"){
      aln_bam <- file.path(res_dir,paste0(samID,"_SJ_aln.bam"))
      sort_bam <- file.path(res_dir,paste0(samID,"_SJ_sort"))
      subjunc(genome_idx,fastq,output_file=aln_bam,
              nthreads=10,unique=FALSE)
      sortBam(aln_bam,sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by subread.")}
  }else if(aligner=="bwa_mem"){
    ci <- unlist(strsplit(basename(genome_idx),split="_"))
    if(ci[length(ci)]=="bwa"){
      aln_sam <- file.path(res_dir,paste0(samID,"_BM_aln.sam"))
      aln_bam <- file.path(res_dir,paste0(samID,"_BM_aln"))
      sort_bam <- file.path(res_dir,paste0(samID,"_BM_sort"))
      system(paste(bwa,"mem","-t",thread,genome_idx,fastq,">",aln_sam,sep=" "))
      asBam(aln_sam,aln_bam)
      sortBam(paste0(aln_bam,".bam"),sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      # unlink(aln_sam)
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by bwa.")}
  }else if(aligner=="bwa_aln"){
    ci <- unlist(strsplit(basename(genome_idx),split="_"))
    if(ci[length(ci)]=="bwa"){
      aln_sai <- file.path(res_dir,paste0(samID,"_BA_aln.sai"))
      aln_sam <- file.path(res_dir,paste0(samID,"_BA_aln.sam"))
      aln_bam <- file.path(res_dir,paste0(samID,"_BA_aln"))
      sort_bam <- file.path(res_dir,paste0(samID,"_BA_sort"))
      system(paste(bwa,"aln","-t",thread,genome_idx,fastq,">",aln_sai,sep=" "))
      system(paste(bwa,"samse","-f",aln_sam,genome_idx,aln_sai,fastq,sep=" "))
      asBam(aln_sam,aln_bam)
      sortBam(paste0(aln_bam,".bam"),sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      # unlink(aln_sam)
      # unlink(aln_sai)
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by bwa.")}
  }else if(aligner=="hisat2"){
    ci <- unlist(strsplit(basename(genome_idx),split="_"))
    if(ci[length(ci)]=="hisat2"){
      file.copy(fastq,res_dir,recursive = TRUE)
      fq <- file.path(res_dir,basename(fastq))
      system(paste("gunzip",fq,sep=" "))
      fq <- file.path(res_dir,gsub(".gz","",basename(fq)))
      aln_sam <- file.path(res_dir,paste0(samID,"_HS2_aln.sam"))
      aln_bam <- file.path(res_dir,paste0(samID,"_HS2_aln"))
      sort_bam <- file.path(res_dir,paste0(samID,"_HS2_sort"))
      # system(paste(hisat,"-p",thread,"-x",genome_idx,"-U",fastq,"-S",aln_sam,sep=" "))
      hisat2(sequences=fq, index=genome_idx,p=thread,
             type="single", outfile=aln_sam,
             force=TRUE, strict=TRUE, execute=TRUE)
      asBam(aln_sam,aln_bam)
      sortBam(paste0(aln_bam,".bam"),sort_bam)
      indexBam(paste0(sort_bam,".bam"))
      # unlink(aln_sam)
      unlink(fq)
      unlink(dir(res_dir,pattern=basename(aln_bam),full.names = TRUE))
    }else{print("The index is not built by hisat2.")}
  }else{print("Please assign an available aligner: bowtie2, subread, subjunc, bwa, and hisat2")}}

#' Revmap_process
#'
#' Process reads to remove linker artifacts or filter  by length
#'
#'
#' @docType methods
#' @name revmap_process
#' @rdname revmap_process
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas path to FASTA file to process.
#' @param length_min minimum length required during fasta processing, default is NULL, set to integer to specify.
#' @param length_max maximum length required during fasta processing, default is NULL, set to integer to specify. 
#' @param linkers contaminating sequences to remove, default is NULL, set to character string to specify.
#' @examples
#' testFastq <- system.file("extdata/SRR1742056.fastq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFastq,overwrite = TRUE)
#' FaFile <- fastx_qtoa(FqFile)
#' processed_Fa <-  revmap_process(FaFile, linkers="GGACGATGC", length_min=18, length_max=30)
#' @return path to processed FASTA.
#' @import Biostrings IRanges
#' @export
revmap_process <- function(fastas, linkers = NULL, length_max = NULL, length_min = NULL) { 
  fa <- readDNAStringSet(fastas, format = "fasta", nrec = -1L)
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
  outname = paste(fastas, sep = "")
  outname = gsub(".fa$", "_processed.fa", outname)
  writeXStringSet(fa, outname)
  return(outname)
}



#' revmap_count
#'
#'Reverse map small RNAs to processed or unprocessed reads
#'
#'
#' @docType methods
#' @name revmap_count
#' @rdname revmap_count
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas path to FASTA files to reverse map.
#' @param knownMiRNAs path to FASTA file containing annotated miRNA or other small RNA sequences.
#' @param linkers contaminating sequences to remove, default is NULL, set to character string to specify. 
#' @param length_min minimum length required during fasta processing, default is NULL, set to integer to specify. 
#' @param length_max maximum length required during fasta processing, default is NULL, set to integer to specify. 
#' @param removedups remove multiple miRNAs or small RNAs mapping to the same read, TRUE (default) or FALSE. If TRUE, known miRNAs will be prioritized and known miRNA names must be in the format "miR-", "let-", "bantam-", "iab-".
#' @param verbose print messages, TRUE or FALSE (default).
#' @param bpparam TRUE or FALSE (default).
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).  
#' @return path to BEDs where small RNAs mapped to reads. 
#' @examples 
#' testFastq <- system.file("extdata/SRR1742056.fastq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFastq,overwrite = TRUE)
#' FaFile <- fastx_qtoa(FqFile)
#' FaFile_clip <- fastx_clipper(FaFile, writelog = FALSE)
#' testMiRNA <- system.file("extdata/hsa_mature.fa",package="CLIPflexR")
#' revmap <- revmap_count(FaFile_clip, testMiRNA, 
#' linkers="GGACGATGC", length_min=18, length_max=30, overwrite=TRUE)
#' @import GenomicAlignments BiocParallel GenomicRanges Biostrings
#' @importFrom magrittr "%>%"
#' @importFrom rtracklayer import
#' @export
revmap_count <- function(fastas, knownMiRNAs, bpparam=NULL,verbose=FALSE, linkers = NULL, length_max = NULL, length_min = NULL, removedups =TRUE, overwrite = FALSE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  pro_fastas <- vector("list",length = length(fastas)) 
  if(!is.null(linkers) | !is.null(length_max) | !is.null(length_min)) { 
    if (verbose) message("Processing fastas..",appendLF = FALSE)
    if (!is.null(linkers) & verbose) message("removing linkers... ", linkers)
    if (!is.null(length_max) & verbose) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min) & verbose) message ("getting sequences longer than ", as.character(length_min), " nt")
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }} else {
      pro_fastas <-  fastas
    }
  if(verbose) message("Creating indices from FASTA files..",appendLF = FALSE)
  indicies <- bplapply(pro_fastas,bowtie2_index,overwrite=overwrite, BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Mapping miRNAs to processed reads..",appendLF = FALSE)
  revBams <- bplapply(indicies,
                      function(x,knownMiRNAs){
                        bowtie_align(knownMiRNAs,index=x, mode = "reverse_map", overwrite = overwrite)
                      },knownMiRNAs=knownMiRNAs, BPPARAM=bpparam)
  if(verbose) message("done")
    if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
  bedpath <- bplapply(revBams, bamtobed,BPPARAM=bpparam)
  dedup<- lapply(bedpath ,read.delim,  header = F, sep = "\t")
  names(dedup) <- bedpath
  if (removedups &  verbose) { message("deduplicating...") 
   for (i in 1:length(dedup)) {
      dedup[[i]]$miRNAnum  <- ifelse(!grepl("miR|let|iab|bantam", dedup[[i]]$V4), stringr::str_sub(dedup[[i]]$V4, 2, 7), NA)
      dedup[[i]]$miRNAnum <- ifelse(grepl("miR|let|iab|bantam", dedup[[i]]$V4), as.numeric(apply(dedup[[i]],1,function(x) regmatches(x["name"],regexpr("[0-9]+",x["name"])))), paste(dedup[[i]]$miRNAnum))
      dedup[[i]] <- dedup[[i]][order(dedup[[i]]$V1, dedup[[i]]$V2, dedup[[i]]$miRNAnum,dedup[[i]]$V4),]
      dedup[[i]] <- dedup[[i]][!duplicated(dedup[[i]]$V1),] 
      dedup[[i]]$miRNAnum <- NULL
    } 
    names(dedup) <- gsub(".bed", "_deduplicated.bed", bedpath)
    for (i in seq_along(dedup)) {
      write.table(dedup[i],names(dedup)[i], col.names =  FALSE, row.names= FALSE, sep = "\t", quote = F)
    } }
  count <- lapply(dedup, function(x) x %>% dplyr::group_by(V4) %>% dplyr::summarize(count=dplyr::n()))
  count_stat <- Reduce(function(x, y) merge(x, y, by = "V4", all = TRUE), count )
  count_stat[is.na(count_stat)] <- 0 
  col.names<- basename(unlist(bedpath)) 
  col.names <- c("miRNA", col.names)
  count_stat <- setNames(count_stat, col.names)
  write.table(count_stat,file.path(dirname(fastas)[1], "revmap_counts.txt"), col.names = TRUE, row.names= FALSE, sep = "\t", quote = F)
  if(verbose) message("done")  
   p <- c(names(dedup), file.path(dirname(fastas)[1], "revmap_counts.txt"))
    return(p)
  }

#' Ranges_count
#'
#' Map unprocessed or processed reads to genome and count small RNAs by location
#'
#'
#' @docType methods
#' @name Ranges_count
#' @rdname Ranges_count
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas path to FASTA files to process.
#' @param miRNA_ranges BED file containing annotated miRNA or other small RNA genomic ranges.
#' @param linkers contaminating sequences to remove, default is NULL, set to character string to specify. 
#' @param length_min minimum length required during fasta processing, default is NULL, set to integer to specify. 
#' @param length_max maximum length required during fasta processing, default is NULL, set to integer to specify. 
#' @param genomeIndex path to genome index. 
#' @param mode mapping mode, default is NULL to allow custom mapping with default arguments set below (0 mismatches in seed alignment, seed substring length 18, disallow alignments with mismatches); see bowtie_align for other mapping modes.
#' @param maxMismatches max mismatches in seed alignment (default is 0).
#' @param seedSubString length of seed substrings (default is 18).
#' @param threads number of threads to use (default is 1).
#' @param report_k number of alignments to report, default is NULL (the best alignment is reported), set to integer to specify.
#' @param keep_all set to TRUE to report all alignments, default is NULL (the best alignment is reported).
#' @param soft_clip set to TRUE to allow soft clipping, default is no soft clipping.
#' @param additional_Args any additional mapping arguments, default is "--score-min C,0,0" to disallow alignments with mismatches; run "Rbowtie2::bowtie2_usage()" to see all options; please note that due to a Rbowtie2 bug, the "--no-1mm-upfront" is not available.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param overwrite overwrite existing output files, TRUE or FALSE (default). 
#' @param bpparam TRUE or FALSE (default)
#' @examples
#' testFastq <- system.file("extdata/SRR1742056.fastq.gz",package="CLIPflexR")
#' @return path to count matrix.
#' @import GenomicAlignments BiocParallel GenomicRanges Biostrings
#' @importFrom rtracklayer import
#' @export
#' 
Ranges_count <- function(fastas,miRNA_ranges,genomeIndex,linkers = NULL, length_max = NULL, length_min = NULL, 
                         mode = NULL, maxMismatches=0, threads=1, report_k = NULL, keep_all  = NULL, soft_clip = NULL,
                         additional_Args = "--score-min C,0,0", seedSubString = 18, overwrite=FALSE, 
                         bpparam=NULL,verbose=FALSE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(!is.null(linkers) | !is.null(length_max) | !is.null(length_min)) { 
    if (verbose) message("Processing fastas..",appendLF = FALSE)
    if (!is.null(linkers) & verbose) message("removing linkers... ", linkers)
    if (!is.null(length_max) & verbose) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min) & verbose) message ("getting sequences longer than ", as.character(length_min), " nt")
    pro_fastas <- vector("list",length = length(fastas)) 
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }} else {
      pro_fastas <-  fastas
    }
  if(verbose) message("Mapping reads to genome..",appendLF = FALSE)
  mappedBams <- bplapply(pro_fastas, bowtie_align, index = genomeIndex, mode=mode, report_k=report_k, threads=threads, maxMismatches=maxMismatches, keep_all=keep_all, soft_clip=soft_clip, seedSubString=seedSubString, additional_Args=additional_Args, overwrite=overwrite, BPPARAM=bpparam)
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
  outname <- paste0(outname, "/range_count_matrix.txt")
  write.table(kksMat, outname,  col.names = TRUE,  row.names = FALSE,  sep =  "\t", quote = F)
  return(outname)
  if (verbose) message("done")
}

#' extract_unmapped
#'
#' Extracts unmapped reads from BAM and writes them to a FASTA file
#'
#'
#' @docType methods
#' @name extract_unmapped
#' @rdname extract_unmapped
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param bam path to BAM file.
#' @param outfa path to output FASTA (default is same directory as input BAM).
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta, overwrite=TRUE)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20")
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5, inputFormat="fastq")
#' bam <- bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq")
#' extract_unmapped(bam)
#' @return path to FASTA file of unmapped reads.
#' @import GenomicAlignments Biostrings
#' @export
extract_unmapped <- function(bam,outfa=NULL){
  inBAM <- scanBam(bam,param=ScanBamParam(what = c("qname","seq"),flag = scanBamFlag(isUnmappedQuery = TRUE)))
  if(is.null(outfa)) outfa <- gsub(".bam","_unmapped.fa",bam)
  toWrite <- inBAM[[1]]$seq
  names(toWrite) <- inBAM[[1]]$qname
  writeXStringSet(toWrite,filepath = outfa,append = FALSE)
  return(outfa)
}





#' chimera_Process
#'
#'Chimera process CLIP data
#'
#'
#' @docType methods
#' @name chimera_Process
#' @rdname chimera_Process
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param bams path to BAM files mapped to the genome, unmapped reads will be extracted for chimera processing.
#' @param knownMiRNAs path to FASTA file containing annotated miRNA or other small RNA sequence. Known miRNAs will be prioritized and known miRNA names must be in the format "miR-", "let-", "bantam-", "iab-".  
#' @param genomeIndex path to genome index.
#' @param removedups remove multiple small RNAs mapping to the same read, TRUE (default) or FALSE. If TRUE, known miRNAs will be prioritized and known miRNA names must be in the format "miR-", "let-", "bantam-", "iab-".
#' @param exclude names of small RNAs to remove, default is NULL, can be set to character vector to specify; must match names in knownMiRNAs file.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param overwrite overwrite existing output files, TRUE or FALSE (default).   
#' @param bpparam TRUE or FALSE (default).
#' @return path to chimera output table.
#' @examples 
#' testFastq <- system.file("extdata/SRR1742056.fastq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFastq,overwrite=TRUE)
#' FaFile <- fastx_qtoa(FqFile)
#' FaFile_clip <- fastx_clipper(FaFile, writelog=FALSE)
#' myGenome <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(myGenome, overwrite = TRUE)
#' myBam <- bowtie_align(FaFile_clip,myIndex, overwrite=TRUE)
#' miRNAs <- system.file("extdata/hsa_mature.fa",package="CLIPflexR")
#' chimera_bed <- chimera_Process(myBam, miRNAs, myIndex, exclude="hsa-miR-19a-3p", overwrite=TRUE)
#' @import GenomicAlignments BiocParallel stringr
#' @importFrom tibble rownames_to_column
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @importFrom purrr map2
#' @export
chimera_Process <- function(bams,knownMiRNAs,genomeIndex,exclude=NULL, bpparam=NULL,verbose=FALSE, removedups = TRUE, overwrite=FALSE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(verbose) message("Extracting unmapped reads to FASTA...",appendLF = FALSE)
  fastas <- bplapply(bams,extract_unmapped,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Creating indices from FASTA files...",appendLF = FALSE)
  indicies <- bplapply(fastas,bowtie2_index,overwrite=overwrite, BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Mapping miRNAs to unmapped reads...",appendLF = FALSE)
  newBams <- bplapply(indicies,
                      function(x,knownMiRNAs){
                        bowtie_align(knownMiRNAs,index=x,mode = "reverse_map", overwrite=overwrite)
                      },knownMiRNAs=knownMiRNAs,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Converting BAMs to BEDs...",appendLF = FALSE)
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
  if (verbose) message("Checking read names...", appendLF = FALSE)
  checkreads <- lapply(fasta, function(x) duplicated(x$rowname))
  checkreads <- unlist(lapply(checkreads, function(x) length(which(x==TRUE))))
  idx <- which(checkreads > 0)
  if (verbose & isTRUE(length(idx)==0)) message("read names ok") else message(paste0("Warning, the following sample(s) have duplicate read names! \n", names(idx)))
  #merge together read sequence and bed by rowname (read name)
  if (verbose) message("Merging reads with mapped small RNAs to read sequences...", appendLF = FALSE)
  chimera <- purrr::map2(fasta,chimera, ~merge(.x,.y, by = "rowname"))
  if(verbose) message("done") 
  #write files to chimera inputs: bed with read name, start, stop, strand, read name, and read sequence
  # setwd(Dir)
  if (verbose) message("Writing merged output files...", appendLF = FALSE)
  fastaTxts <- vector("list",length = length(chimera))
  for (x in 1:length(chimera)) {
    write.table(chimera[[x]], file=paste0(names(chimera)[x],".txt"), sep="\t", quote = FALSE)
    fastaTxts[[x]] <- paste0(names(chimera)[x],".txt")
  }
  if(verbose) message("done")
  names(fastaTxts) <- names(chimera)
  chimerafiles <- vector("list",length = length(fastaTxts))
  if (verbose) message("Extracting read sequences up- and downstream of small RNA sequences & writing output files...", appendLF = FALSE)
  for (i in 1:length(fastaTxts)) {
    chimerafiles[[i]] <- chimeraProcess(fastaTxts[[i]], exclude=exclude, removedups=removedups)
  }
  if(verbose) message("done")
  rffiles <- vector("list",length = length(chimerafiles)) 
  if(verbose) message("Writing downstream sequences as fastas...", appendLF = FALSE)
  for (i in 1:length(chimerafiles)) {
    rffiles[[i]] <- reformat(chimerafiles[[i]])
  }
  if(verbose) message("done")
  
  if(verbose) message("Remapping downstream sequences to genome...", appendLF = FALSE)
  remappedBams <- bplapply(rffiles, bowtie_align, index =genomeIndex,mode =  "genome_map", overwrite=overwrite, BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Converting remapped BAMs to BEDs...", appendLF = FALSE)
  myBeds <- bplapply(remappedBams, bamtobed,BPPARAM=bpparam)
  if(verbose) message("done")
  
  return(myBeds)
  
}

chimeraProcess <- function(input, exclude, removedups) {
  BR1 <- read.delim(input, header=T)
  BR1 <- BR1[!BR1$name %in% exclude,] 
  BR1<- BR1[BR1$strand=="+",]
  if(removedups){
  BR1$miRNAnum <- ifelse(!grepl("miR|let|iab|bantam", BR1$name), stringr::str_sub(BR1$name, 2, 7), NA)
  BR1$miRNAnum <- ifelse(grepl("miR|let|iab|bantam", BR1$name), as.numeric(apply(BR1,1,function(x) regmatches(x["name"],regexpr("[0-9]+",x["name"])))), paste(BR1$miRNAnum))
  BR1 <- BR1[order(BR1$rowname, BR1$start, BR1$miRNAnum, BR1$name),]
  BR1 <- BR1[!duplicated(BR1$rowname),] 
  BR1$miRNAnum <- NULL}
  BR1$length <- sapply(as.character(BR1$x), nchar)
  BR1$ups.seq <- mapply(substr, x=BR1$x, start=0, stop=BR1$start)
  BR1$dns.seq <- mapply(substr, x=BR1$x, start=BR1$stop+1, stop=BR1$length)
  outname = paste(input, sep = "")
  outname = gsub(".fa.txt", "_chimera.txt", outname)
  write.table(BR1, outname, quote=F, sep="\t", col.names = T, row.names = F)
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

