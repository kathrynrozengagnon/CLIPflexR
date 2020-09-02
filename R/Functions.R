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
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
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
#' @param fileToDecompress File to bzip
#' @param outDir Output directory
#' @param keep keep (don't delete) input files
#' @param overwrite overwrite existing output files
#' @return Path to unzipped file
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
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
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
#' @name fetchSequencesForCLIP
#' @rdname fetchSequencesForCLIP
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param peaks CLIP peaks to process.
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
#' @param peaks CLIP peaks to process
#' @param fasta Genome fasta file
#' @param patterns Patterns to scan in CLIP peaks
#' @param resize Size of window in bp for resizing around peak center
#' @param add5 Bp to add to 5' of resized peak
#' @param add3 Bp to add to 3' of resized peak
#' @param verbose Print messages, TRUE (default) or FALSE
#' @param checkReverse Check the reverse strand for pattern.
#' @param bedHeader TRUE if peak file contains column headers.
#' @param bedSep Seperator in BED file.
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
    sitesFA <- c(fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta),reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta)))
    sitesFAExtended <- c(fetchSequencesForCLIP(sitesGRpos,resize=NULL,fasta,add5=add5,add3=add3),
                         reverseComplement(fetchSequencesForCLIP(sitesGRneg,resize=NULL,fasta,add5=add5,add3=add3))
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
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
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
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
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

#' make count matrix from bed
#'
#' make count matrix from bed
#'
#'
#' @docType methods
#' @name countFromBed
#' @rdname countFromBed
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param Bed BED bed files containing mapped reads
#' @param GR BED file containing genomic ranges over which to count 
#' @param notStranded Stranded or not stranded TRUE (default) or FALSE
#' @param interFeature Count reads mapping to multiple features TRUE or FALSE (default)
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' bowtie2_index(testFasta)
#' @import rtracklayer GenomicRanges GenomicAlignments
#' @return counting matrix
#' @export
countFromBed <- function(Bed,GR,notStranded=TRUE,interFeature=FALSE){
  reads <- import.bed(Bed)
  fk <- summarizeOverlaps(GR,reads,ignore.strand = notStranded,inter.feature=interFeature)
  assay(fk)
}


#' Build bigwigs
#'
#' @rdname CLIP_bw2
#' @param sort_bam sorted bam file
#' @param res_dir result directory
#' @param normalized Normalised to total reads
#' 
#' @docType methods
#' @import  Rsamtools
#' @return bigwig
#' @export
CLIP_bw2 <- function(sort_bam,res_dir=NULL,normalized=TRUE){
  if(dir.exists(res_dir)){dir.create(res_dir)}else{print(paste0("folder is exist: ",res_dir))}
  samID <- gsub("_sort.bam","",basename(sort_bam))
  bw_out <- file.path(res_dir,paste0(samID,".bigwig"))
  if(isTRUE(normalized)){
    allChromosomeStats <- idxstatsBam(sort_bam)
    mapped <- sum(allChromosomeStats[,"mapped"])
    toRun <- coverage(BamFile(sort_bam),weight = (10^6)/mapped)
    export.bw(toRun,con=bw_out)
  }else{
    toRun <- coverage(BamFile(sort_bam))
    export.bw(toRun,con=bw_out)}
}

#' Build index for reference genome
#'
#' @rdname CLIP_buildIDX
#' @param ref_seq reference genome
#' @param res_dir result directory
#' @param preFIX file name prefix for indexed genome
#' @param aligner Assign aligner bowtie2/subread/bwa/hisat2
#' @param thread threads used in processing, default is 4
#' @param bwa exercutable bwa file path
#' @import Rbowtie2 Rsubread Rhisat2
#' 
#' @docType methods
#' @return indexed reference genome
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
#' @param samID sample id in the sample sheet
#' @param samSheet sample sheet
#' @param res_dir result directory
#' @param genome_idx file path of indexed genome
#' @param aligner Assign aligner bowtie2/bwa_mem/bwa_alb/subread/subjunc/hisat2
#' @param thread thread used in the processing, default is 4
#' @param bwa executable file path of bwa
#' @import Rbowtie2 Rsamtools Rsubread ShortRead 
#' 
#' @docType methods
#' @return sorted bam file
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

#' Install ctk pipeline
#'
#' Install ctk pipeline
#'
#'
#' @docType methods
#' @name install_ctk
#' @rdname install_ctk
#'
#' @author Kathryn Rozen-Gagnon Thomas Carroll Ji-Dung Luo
#'
#' @param path Path to where to install ctk and czplib
#' @import utils
#' @examples 
#' install_ctk()
#' getOption("CLIPflexR.condaEnv")
#' getOption("CLIPflexR.ctk")
#' getOption("CLIPflexR.czplib")
#' @export
install_ctk <- function(path=NULL){
  tempdir <- tempdir()
  miniCondaPath <- miniconda_path()
  miniCondaPathExists <- miniconda_exists(miniCondaPath)
  clipr <- file.path(miniCondaPath,"envs",paste0("CLIPflexR","_",packageVersion("CLIPflexR")))
  if(dir.exists(clipr)) path <- clipr
  if(is.null(path) & !is.null(getOption("CLIPflexR.condaEnv"))) path <- clipr
  if(is.null(path) & is.null(getOption("CLIPflexR.condaEnv"))) path <- getwd()
  download.file(url = "https://github.com/chaolinzhanglab/czplib/archive/master.zip",
                destfile = file.path(tempdir,"czplib-master.zip"))
  download.file(url = "https://github.com/chaolinzhanglab/ctk/archive/master.zip",
                destfile = file.path(tempdir,"ctk-master.zip"))
  utils::unzip(file.path(tempdir,"czplib-master.zip"),exdir = file.path(tempdir))
  utils::unzip(file.path(tempdir,"ctk-master.zip"),exdir = file.path(tempdir))
  czplipCopy <- list.files(file.path(tempdir,"czplib-master"),recursive=TRUE)
  ctkCopy <- list.files(file.path(tempdir,"ctk-master"),recursive=TRUE)
  dir.create(file.path(path,"lib","czplib"),recursive = TRUE,showWarnings = FALSE)
  dir.create(file.path(path,"bin","ctk"),recursive = TRUE,showWarnings = FALSE)
  file.copy(file.path(tempdir,"czplib-master",czplipCopy),file.path(path,"lib","czplib"),recursive = TRUE,copy.mode = TRUE)
  file.copy(file.path(tempdir,"ctk-master",ctkCopy),file.path(path,"bin","ctk"),recursive = TRUE,copy.mode = TRUE)
  Sys.chmod(dir(file.path(path,"lib","czplib"),include.dirs = TRUE,recursive = TRUE,full.names = TRUE),mode = "0755")
  Sys.chmod(dir(file.path(path,"bin","ctk"),include.dirs = TRUE,recursive = TRUE,full.names = TRUE),mode = "0755")
}

#' Process reads for small RNA counting
#'
#' Process reads for small RNA counting
#'
#'
#' @docType methods
#' @name revmap_process
#' @rdname revmap_process
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas FASTA file to process
#' @param length_min minimum length required during fasta processing, integer (default is NULL)
#' @param length_max maximum length required during fasta processing, integer (default is NULL) 
#' @param linkers contaminating seqences to remove (default is NULL)
#' @examples
#' bam <- system.file("extdata/xxbam",package="CLIPflexR")
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped,paste0(FqFile_QF,".gz"))
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5)
#' bam <- bowtie_align(FqFile_QFColStripped,myIndex)
#' unbam(bam)
#' @return Processed FASTAS
#' @import reticulate Biostrings IRanges
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



#' Reverse map small RNAs to processed reads
#'
#'Reverse map small RNAs to processed reads
#'
#'
#' @docType methods
#' @name revmap_count
#' @rdname revmap_count
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fastas FASTA files to reverse map
#' @param knownMiRNAs FASTA file containing annotated miRNA or other small RNA sequences
#' @param linkers contaminating sequences to remove 
#' @param length_min minimum length required during fasta processing, integer (default is NULL)
#' @param length_max maximum length required during fasta processing, integer (default is NULL)
#' @param verbose print messages, TRUE (default) or FALSE 
#' @param bpparam TRUE or FALSE (default)
#' @param removedups remove multiple miRNAs or small RNAs mapping to the same read, TRUE or FALSE (FALSE)
#' @return BEDs containing counts of  
#' @import GenomicAlignments BiocParallel stringr Biostrings
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @export
revmap_count <- function(fastas, knownMiRNAs, bpparam=NULL,verbose=TRUE, linkers = NULL, length_max = NULL, length_min = NULL, removedups =FALSE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(!is.null(linkers) | !is.null(length_max) | !is.null(length_min)) { 
    if (verbose) message("Processing fastas..",appendLF = FALSE)
    if (!is.null(linkers)) message("removing linkers... ", linkers)
    if (!is.null(length_max)) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min)) message ("getting sequences longer than ", as.character(length_min), " nt")
    pro_fastas <- vector("list",length = length(fastas)) 
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }} else {
      pro_fastas <-  fastas
    }
  if(verbose) message("Creating indices from FASTA files..",appendLF = FALSE)
  indicies <- bplapply(pro_fastas,bowtie2_index,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Mapping miRNAs to processed reads..",appendLF = FALSE)
  revBams <- bplapply(indicies,
                      function(x,knownMiRNAs){
                        bowtie_align(knownMiRNAs,index=x, bam = paste0(x, ".bam"), maxMismatches = 0, report_k = 1000000)
                      },knownMiRNAs=knownMiRNAs,BPPARAM=bpparam)
  if(verbose) message("done")
  if(verbose) message("Converting BAMs to BEDs..",appendLF = FALSE)
  beds <- bplapply(revBams, bamtobed,BPPARAM=bpparam)
  if (removedups) { message("deduplicating...") 
    dedup <- lapply(beds,read.delim,  header = F)
    for (i in 1:length(dedup)) {
      dedup[[i]]$miRNAnum  <- ifelse(!grepl("miR|let|iab", dedup[[i]]$V4), stringr::str_sub(dedup[[i]]$V4, 2, 7), NA)
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

#' Map unprocessed or processed reads to genome and count small RNAs by location
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
#' @param fastas FASTA files to process.
#' @param miRNA_ranges BED file containing annotated miRNA or other small RNA genomic ranges
#' @param linkers contaminating sequences to remove (default is NULL) 
#' @param length_min minimum length required during fasta processing (default is NULL)
#' @param length_max maximum length required during fasta processing (default is NULL)
#' @param genomeIndex full genome index 
#' @param verbose print messages, TRUE (default) or FALSE 
#' @param bpparam TRUE or FALSE (default)
#' @return Path 
#' @import GenomicAlignments BiocParallel stringr rtracklayer GenomicRanges
#' @importMethodsFrom rtracklayer export.bed export.bw mcols
#' @export
#' 
Ranges_count <- function(fastas,miRNA_ranges,genomeIndex,linkers = NULL, length_max = NULL, length_min = NULL, bpparam=NULL,verbose=TRUE){
  if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
  if(!is.null(linkers) | !is.null(length_max) | !is.null(length_min)) { 
    if (verbose) message("Processing fastas..",appendLF = FALSE)
    if (!is.null(linkers)) message("removing linkers... ", linkers)
    if (!is.null(length_max)) message ("getting sequences shorter than ", as.character(length_max), " nt")
    if (!is.null(length_min)) message ("getting sequences longer than ", as.character(length_min), " nt")
    pro_fastas <- vector("list",length = length(fastas)) 
    for (i in 1:length(fastas)) {
      pro_fastas[[i]] <- revmap_process(fastas[[i]], length_max = length_max, linkers =  linkers, length_min = length_min)
    }} else {
      pro_fastas <-  fastas
    }
  if(verbose) message("Mapping reads to genome..",appendLF = FALSE)
  mappedBams <- bplapply(pro_fastas, bowtie_align, index = genomeIndex, maxMismatches=0,BPPARAM=bpparam)
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
  outname <- paste0(outname, "/miRNA_count_matrix.txt")
  write.table(kksMat, outname,  col.names = TRUE,  row.names = FALSE,  sep =  "\t", quote = F)
  return(outname)
  if (verbose) message("done")
}

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
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <- bowtie2_index(testFasta)
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
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
  inBAM <- scanBam(bam,param=ScanBamParam(what = c("qname","seq"),flag = scanBamFlag(isUnmappedQuery = TRUE)))
  if(is.null(outfa)) outfa <- gsub(".bam","_unmapped.fa",bam)
  toWrite <- inBAM[[1]]$seq
  names(toWrite) <- inBAM[[1]]$qname
  writeXStringSet(toWrite,filepath = outfa,append = FALSE)
  return(outfa)
}





#' Chimera process CLIP data
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
#' @param bams BAM files mapped to the genome to process
#' @param knownMiRNAs FASTA file containing annotated miRNA or other small RNA sequences
#' @param genomeIndex Full genome index
#' @param exclude Names of miRNAs or small RNA sequences to remove, character vector 
#' @param verbose print messages, TRUE (default) or FALSE 
#' @param bpparam TRUE or FALSE (default)
#' @return Path 
#' @import GenomicAlignments BiocParallel purrr stringr tibble
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
  newBams <- bplapply(indicies,
                      function(x,knownMiRNAs){
                        bowtie_align(knownMiRNAs,index=x, bam = paste0(x, ".bam"),maxMismatches=0,seedSubString=18,report_k=1000000)
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

