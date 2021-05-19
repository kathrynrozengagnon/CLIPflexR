#' Wrapper function for ctk's stripBarcode
#'
#' Wrapper function for ctk's stripBarcode
#'
#'
#' @docType methods
#' @name ctk_stripBarcode
#' @rdname ctk_stripBarcode
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (fastq or fasta).
#' @param outFile path to output file (fastq or fasta).
#' @param sb path to stripBarcode.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param linkerlength length of barcode/linker sequences (default is 27).
#' @param inputFormat input file format, "fasta" (default) or "fastq"
#' @param barcodeStartWith filter sequences based on the starting nucleotides in the barcode.
#' @param barcodeEndWith filter sequences based on the ending nucleotides in the barcode.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' @return path to file with index stripped, in the same format as the input file.
#' @importFrom reticulate miniconda_path
#' @export
ctk_stripBarcode <- function(filesToRun,
                             outFile=paste(file_path_sans_ext(fileToRun),"_rm5.",file_ext(fileToRun),sep=""),
                             sb="stripBarcode.pl",
                             perl="perl",
                             PATHTOPERLLIB=NULL,
                             linkerlength=27,
                             inputFormat="fasta",
                             barcodeStartWith=NULL,
                             barcodeEndWith=NULL,
                             stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_stripBarcode_stderr.txt")),
                             stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_stripBarcode_stdout.txt")),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  if (grepl("fa|fasta|fastq|fq",file_ext(outFile))) {
    outFile <- outFile
  } else if(inputFormat=="fasta") { outFile <- paste0(outFile,  "fa")}
  else if(inputFormat=="fastq"){ outFile <- paste0(outFile,  "fastq")}
  args <- c(cmd,
            paste0("-len ",linkerlength),
            paste0("-format ",inputFormat),
            ifelse(!is.null(barcodeStartWith),paste0("-barcode-start-with ",barcodeStartWith),""),
            ifelse(!is.null(barcodeEndWith),paste0("-barcode-end-with ",barcodeEndWith),""),
            fileToRun,
            gsub("\\.gz$","",outFile))
  args <- args[args != ""]
  if(verbose){
    
    message("stripBarcode.pl command is ",cmd)
    message("stripBarcode.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}


#' Wrapper function for ctk's CIMS
#'
#' Wrapper function for ctk's CIMS
#'
#'
#' @docType methods
#' @name ctk_cims
#' @rdname ctk_cims
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param mutationBedFile mutation BED file.
#' @param outFile path to output file (BED & text).
#' @param sb path to CIMS.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param bigFile big file, TRUE or FALSE (default).
#' @param mutationSize mutation size (default is 1).
#' @param permutations number of iterations for permutation (default is 5).
#' @param trackMutationPos track mutation position relative to read start, TRUE or FALSE (default).
#' @param noSparseCorrect no sparcity correction, TRUE or FALSE (default).
#' @param FDR threshold of FDR (default is 1).
#' @param mfr threshold of m-over-k-ratio (default is 0).
#' @param cacheDir name for cache directory.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' \dontrun{
#' mutations <- system.file("extdata/BrdU.Fox.pool.tag.uniq.mutation.small.txt",package="CLIPflexR")
#' delBed <- ctk_getMutationType(mutations)
#' ctk_cims("~/Downloads/uniq_tags_mutations/Fox.pool.tag.uniq.rgb.bed",delBed,verbose=TRUE)
#' }
#' @return path to CIMS text & bed files.
#' @export
ctk_cims <- function(filesToRun,
                     mutationBedFile,
                     outFile=paste(file_path_sans_ext(fileToRun),"CIMS","txt",sep="."),
                     sb="CIMS.pl",
                     perl="perl",
                     PATHTOPERLLIB=NULL,
                     bigFile=FALSE,
                     mutationSize=1,
                     permutations=5,
                     trackMutationPos=FALSE,
                     noSparseCorrect=FALSE,
                     FDR=1,
                     mfr=0,
                     cacheDir=paste0(filesToRun,
                                     "_cache",paste(sample(letters,10),collapse="")),
                     stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_cims_stderr.txt")),
                     stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_cims_stdout.txt")),
                     useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                     additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            paste0("-w ",mutationSize),
            paste0("-n ",permutations),
            paste0("-FDR ",FDR),
            paste0("-mkr ",mfr),
            paste0("-c ",cacheDir),
            ifelse(bigFile,"-big",""),
            ifelse(trackMutationPos,"-p",""),
            ifelse(noSparseCorrect,"---no-sparse-correct",""),
            fileToRun,
            mutationBedFile,
            outFile
  )
  args <- args[args!=""]
  if(verbose){
    
    message("CIMS.pl command is ",cmd)
    message("CIMS.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
    
  }
  return(outFile)
}


#' Wrapper function for ctk's CITS
#'
#' Wrapper function for ctk's CITS
#'
#'
#' @docType methods
#' @name ctk_cits
#' @rdname ctk_cits
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile output file (BED).
#' @param sb path to CITS.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param bigFile big file, TRUE or FALSE (default).
#' @param pCutOff p-value cut off (default is 0.01).
#' @param multiTest perform multiple testing correction, TRUE (default) or FALSE.
#' @param gap gap size (default is 25).
#' @param cacheDir name for cache directory.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_Col<- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' @return path to CITS bed file.
#' @export
ctk_cits <- function(filesToRun,
                     outFile=paste(file_path_sans_ext(fileToRun),"CITS","bed",sep="."),
                     sb="CITS.pl",
                     perl="perl",
                     PATHTOPERLLIB=NULL,
                     bigFile=FALSE,
                     pCutOff=0.01,
                     multiTest=TRUE,
                     gap=25,
                     cacheDir=paste0(filesToRun,
                                     "_cache",paste(sample(letters,10),collapse="")),
                     stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_cits_stderr.txt")),
                     stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_cits_stdout.txt")),
                     useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                     additional_Args=NULL,verbose=FALSE, writelog=T){
  
  
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  if(useClipRConda) perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  if(!is.null(PATHTOPERLLIB)){
    on.exit((function(perl5libPathOld) {
      if(!is.na(perl5libPathOld)){
        Sys.setenv("PERL5LIB"=perl5libPathOld)
      }else{
        Sys.unsetenv("PERL5LIB")
      }
    })(perl5libPathOld))
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            paste0("-p ",pCutOff),
            paste0("--gap ",gap),
            paste0("-c ",cacheDir),
            ifelse(bigFile,"-big ",""),
            ifelse(multiTest,"--multi-test",""),
            fileToRun,
            outFile
  )
  args <- args[args!=""]
  if(verbose){
    
    message("CITS.pl command is ",perl)
    message("CITS.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
    # unlink(cacheDir,recursive = TRUE)
  }
  return(outFile)
}

#' Wrapper function for ctk's getMutationType
#'
#' Wrapper function for ctk's getMutationType
#'
#'
#' @docType methods
#' @name ctk_getMutationType
#' @rdname ctk_getMutationType
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile output file (BED).
#' @param sb path to ctk_getMutationType.pl from CTK toolkit.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param mutationType mutation type, "del" (deletions; default),"ins" (insertions), "sub" (substitutions).
#' @param summaryStat create summary file, TRUE or FALSE (default).
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' mutations <- system.file("extdata/BrdU.Fox.pool.tag.uniq.mutation.small.txt",package="CLIPflexR")
#' ctk_getMutationType(mutations)
#' @return path to getMutationType bed file.
#' @export
ctk_getMutationType <- function(filesToRun,
                                outFile=paste(file_path_sans_ext(fileToRun),mutationType,"bed",sep="."),
                                sb="getMutationType.pl",
                                perl="perl",
                                PATHTOPERLLIB=NULL,
                                mutationType="del",
                                summaryStat=FALSE,
                                stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_getMutationType_stderr.txt")),
                                stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_getMutationType_stdout.txt")),
                                useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  if(useClipRConda) perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  if(!is.null(PATHTOPERLLIB)){
    on.exit((function(perl5libPathOld) {
      if(!is.na(perl5libPathOld)){
        Sys.setenv("PERL5LIB"=perl5libPathOld)
      }else{
        Sys.unsetenv("PERL5LIB")
      }
    })(perl5libPathOld))
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            paste0("-t ",mutationType),
            ifelse(summaryStat,paste0("--summary ",paste(file_path_sans_ext(fileToRun),mutationType,"summary",sep=".")),""),
            fileToRun,
            outFile
  )
  args <- args[args!=""]
  
  if(verbose){
    
    message("getMutationType.pl command is ",cmd)
    message("getMutationType.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
  }
  return(outFile)
}



#' Wrapper function for ctk's fastq_filter
#'
#' Wrapper function for ctk's fastq_filter
#'
#'
#' @docType methods
#' @name ctk_fastqFilter
#' @rdname ctk_fastqFilter
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process  (fastq).
#' @param outFile output file (fastq).
#' @param sb path to fastq_filter.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param fastqFormat fastq format used, can be "sanger" (default) or "solexa".
#' @param indexPosition position and sequence of index in read, default is NULL; set "position:seqeunce" to specify (e.g. "1:CATCGC").
#' @param qsFilter set quality score filter, default is NULL; set method:start-end:score to specify (e.g. "mean:0-29:20"; starts/ends are 0-based).
#' @param maxN maximum number of unknown nucleotides (N) allowed, default is NULL; set to integer to specify. 
#' @param outputFormat output format, "fastq" (default) or "fasta".
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' @return Path to filtered file in specified format.
#' @export
ctk_fastqFilter <- function(filesToRun,
                            outFile=file.path(dirname(fileToRun),paste("FF_",basename(fileToRun),sep="")),
                            sb="fastq_filter.pl",
                            perl="perl",
                            PATHTOPERLLIB=NULL,
                            fastqFormat="sanger",
                            indexPosition=NULL,
                            qsFilter=NULL,
                            maxN=NULL,
                            outputFormat="fastq",
                            stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_fastqFilter_stderr.txt")),
                            stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_fastqFilter_stdout.txt")),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }

  args <- c(cmd,
            paste0("-if ",fastqFormat),
            paste0("-of ",outputFormat),
            ifelse(!is.null(indexPosition),paste0("-index ",indexPosition),""),
            ifelse(!is.null(qsFilter),paste0("-f ",qsFilter),""),
            ifelse(!is.null(maxN),paste0("-maxN ",maxN),""),
            fileToRun,
            outFile
  )
  args <- args[args!=""]
  if(verbose){
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
  }
  return(outFile)
}


#' Wrapper function for ctk's fastq2collapse
#'
#' Wrapper function for ctk's fastq2collapse
#'
#'
#' @docType methods
#' @name ctk_fastq2collapse
#' @rdname ctk_fastq2collapse
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (fastq).
#' @param outFile output file (fastq).
#' @param sb path to fastq2collapse.pl from CTK.
#' @param perl path to PERL
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' @return Path to collapsed fastq file.
#' @export
ctk_fastq2collapse <- function(filesToRun,
                               outFile=file.path(dirname(fileToRun),paste("Collapsed_",basename(fileToRun),sep="")),
                               sb="fastq2collapse.pl",
                               perl="perl",
                               PATHTOPERLLIB=NULL,
                               stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_fastq2collapse_stderr.txt")),
                               stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_fastq2collapse_stdout.txt")),
                               useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                               additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1]
  
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            fileToRun,
            gsub("\\.gz$","",outFile))
  
  args <- args[args!=""]
  if(verbose){
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")}
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}


#' Wrapper function for ctk's parseAlignment
#'
#' Wrapper function for ctk's parseAlignment
#'
#'
#' @docType methods
#' @name ctk_parseAlignment
#' @rdname ctk_parseAlignment
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (SAM).
#' @param outFile path to output file (BED).
#' @param sb path to parseAlignment.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param mutationFile mutation file path, default is NULL.
#' @param mapQual minimum map quality, default is NULL.
#' @param minLen minimum length of read, default is NULL.
#' @param indelToEnd minimum distance from indel to end of read, default is 5.
#' @param splitDel whether to split reads with deletions, TRUE or FALSE (default).
#' @param indelInScore include indels in mutation score count, TRUE or FALSE (default).
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose = TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' parsedAlignment <- ctk_parseAlignment(bam)
#' @return path to BED file.
#' @export
ctk_parseAlignment <- function(filesToRun,
                               outFile=paste0(file_path_sans_ext(filesToRun),".bed"),
                               sb="parseAlignment.pl",
                               perl="perl",
                               PATHTOPERLLIB=NULL,
                               mutationFile=NULL,
                               mapQual=NULL,
                               minLen=NULL,
                               indelToEnd=5,
                               splitDel=FALSE,
                               indelInScore=FALSE,
                               stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_parseAlignment_stderr.txt")),
                               stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_parseAlignment_stdout.txt")),
                               useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                               additional_Args=NULL,verbose=FALSE, writelog  =T){
  
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1] 
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  fileToRun <- Rsamtools::asSam(fileToRun,overwrite=TRUE)
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            ifelse(!is.null(mutationFile),paste0("--mutation-file ",mutationFile),""),
            ifelse(!is.null(mapQual),paste0("--map-qual ",mapQual),""),
            ifelse(!is.null(minLen),paste0("--min-len ",minLen),""),
            paste0("--indel-to-end ",indelToEnd),
            ifelse(!is.null(splitDel),paste0("--split-del "),""),
            ifelse(!is.null(indelInScore),paste0("--indel-in-score "),""),
            fileToRun,
            outFile)
  
  args <- args[args!=""]
  if(verbose){
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")
      }
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}


#' Wrapper function for ctk's tag2collapse
#'
#' Wrapper function for ctk's tag2collapse
#'
#'
#' @docType methods
#' @name ctk_tag2collapse
#' @rdname ctk_tag2collapse
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile path to output file
#' @param sb path to tag2collapse.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param keepMaxScore keep the tag with the most weight (instead of the longest one) as representative.
#' @param keepTagName do not change tag name (no extra information).
#' @param weight consider the weight of each tag.
#' @param bigFile Set to TRUE when files are big.
#' @param weightInName find weight in name.
#' @param randomBarcode random barcode exists, no collapse for different barcodes.
#' @param seqErrorModel sequencing error model to use (alignment or em-local or em-global or fix=0.01).
#' @param outputSeqError output sequencing errors estimated by the EM algorithm.
#' @param em EM threshold to infer reliability of each collapsed read (when have random linker, -1=no EM).
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args Additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' parsedAlignment <- ctk_parseAlignment(bam)
#' ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName = FALSE,verbose = TRUE)
#' @return path to collapsed BED file.
#' @export
ctk_tag2collapse <- function(filesToRun,
                             outFile=file.path(dirname(filesToRun),paste0("TC_",basename(filesToRun))),
                             sb="tag2collapse.pl",
                             perl="perl",
                             PATHTOPERLLIB=NULL,
                             keepMaxScore=TRUE,
                             keepTagName=TRUE,
                             weight=TRUE,
                             bigFile=FALSE,
                             weightInName=TRUE,
                             randomBarcode=TRUE,
                             seqErrorModel=TRUE,
                             outputSeqError=NULL,
                             em=NULL,
                             stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2collapse_stderr.txt")),
                             stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2collapse_stdout.txt")),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additional_Args=NULL,verbose=FALSE, writelog = T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1] 
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  args <- c(cmd,
            ifelse(keepTagName,paste0("--keep-tag-name "),""),
            ifelse(keepMaxScore,paste0("--keep-max-score "),""),
            ifelse(!is.null(outputSeqError),paste0("--output-seq-error ",outputSeqError),""),
            ifelse(!is.null(em),paste0("-EM ",em),paste0("-EM -1")),
            ifelse(!is.null(seqErrorModel),paste0("--seq-error-model ",seqErrorModel),""),
            ifelse(bigFile,paste0("-big  "),""),
            ifelse(weight,paste0("-weight  "),""),
            ifelse(randomBarcode,paste0("--random-barcode  "),""),
            ifelse(weightInName,paste0("--weight-in-name  "),""),
            fileToRun,
            outFile)        
  
  
  args <- args[args!=""]
  if(verbose){
    
    message("tag2collapse.pl command is ",cmd)
    message("tag2collapse.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")
      }
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}

#' Wrapper function for galaxy's joinWrapper
#'
#' Wrapper function for galaxy's joinWrapper
#'
#'
#' @docType methods
#' @name ctk_joinWrapper
#' @rdname ctk_joinWrapper
#'
#' @author Kathryn Rozen-Gagnon
#' 
#' @param file1 path to file 1 to join (BED or tab delimited text).
#' @param file2 path to file 2 to join (BED or tab delimited text).
#' @param field1 field/column for file 1 join.
#' @param field2 field/column for file 2 join.
#' @param mode join mode.
#' @param outFile path to output file.
#' @param sb path to joinWrapper.py from Galaxy
#' @param python path to python.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' \dontrun{
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' mutationFile <- system.file("extdata/Fox3_Std_small_mutation.txt",package="CLIPflexR")
#' parsedAlignment <- ctk_parseAlignment(bam,mutationFile=mutationFile)
#' uniqueTags <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' ctk_joinWrapper(mutationFile,uniqueTags,4,4,"N",verbose=TRUE)
#' }
#' @return path to joined file.
#' @export
ctk_joinWrapper <- function(file1,
                            file2,
                            field1,
                            field2,
                            mode,
                            outFile=file.path(dirname(file1),paste0("Unique_",basename(file1))),
                            sb="joinWrapper.py",
                            python="python",
                            stderr=file.path(dirname(file1),paste0(basename(file1),"_ctk_joinWrapper_stderr.txt")),
                            stdout=file.path(dirname(file1),paste0(basename(file1),"_ctk_joinWrapper_stdout.txt")),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additional_Args=NULL,verbose=FALSE, writelog = T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  # perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    python <- file.path(getOption("CLIPflexR.condaEnv"),"bin",python)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(system.file("extdata/",package="CLIPflexR"),cmd),sep = " ")
  
  fileToRun <- file1
  
  if(!file.exists(fileToRun) | !file.exists(file2)) stop("File/s do not exist")
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(pathOld) {
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(pathOld))
  
  
  
  args <- c(cmd,
            file1,
            file2,
            field1,
            field2,
            mode,
            outFile)
  
  args <- args[args!=""]
  if(verbose){
    
    message("joinWrapper.py command is ",cmd)
    message("joinWrapper.py arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(python,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(python,
                args,
                stdout="",
                stderr="")
      }
  }
  
  return(outFile)
}

#' Wrapper function for ctk's bed2rgb
#'
#' Wrapper function for ctk's bed2rgb
#'
#'
#' @docType methods
#' @name ctk_bed2rgb
#' @rdname ctk_bed2rgb
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile path to output file (BED).
#' @param sb path to bed2rgb.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param col color to include in BED rgb column (default is "blue").
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package ="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package ="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20")
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollapsed <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' ctk_bed2rgb(myCollapsed,col="128,0,0")
#' @return path to BED file with color specified.
#' @export
ctk_bed2rgb <- function(filesToRun,
                        outFile=paste0(file_path_sans_ext(filesToRun),".RGB.",file_ext(filesToRun)),
                        sb="bed2rgb.pl",
                        perl="perl",
                        PATHTOPERLLIB=NULL,
                        col="blue",
                        stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_bed2rgb_stderr.txt")),
                        stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_bed2rgb_stdout.txt")),
                        useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                        additional_Args=NULL,verbose=FALSE, writelog = T){
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1] 
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  
  args <- c(cmd,
            paste0("-col ",col),
            fileToRun,
            outFile)
  
  args <- args[args!=""]
  if(verbose){
    
    message("bed2rgb.pl command is ",cmd)
    message("bed2rgb.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")
      }
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}


#' Wrapper function for ctk's ctk_tag2profile
#'
#' Wrapper function for ctk's ctk_tag2profile
#'
#'
#' @docType methods
#' @name ctk_tag2profile
#' @rdname ctk_tag2profile
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile path to output file (BED).
#' @param outFile2 path to output 2 file (BED; specify two output files to separate strands).
#' @param bigFile TRUE when working with a big file, TRUE or FALSE (default).
#' @param weight weight counts according to the score of each tag, TRUE or FALSE (default).
#' @param weightAvg weight average the score of each tag, TRUE or FALSE (default).
#' @param ss separate strand, TRUE (default) or FALSE.
#' @param exact exact count at each nucleotide, TRUE (default) or FALSE.
#' @param nz don't print zeroes (works for sgr and bed), TRUE or FALSE (default).
#' @param ext5 extension of tags at the 5' end, default is NULL (set to integer to specify).
#' @param ext3 extension of tags at the 3' end, default is NULL (set to integer to specify).
#' @param chromLen chrom length file, default is NULL (set file path to specify).
#' @param region a bed file with regions to count tag numbers; if not specified  (NULL), count in moving windows
#' @param minBlockSize minimum number of lines to read in each block for a big file, default is 2000000.
#' @param windowSize window size, default is 100.
#' @param stepSize step size, default is 20.
#' @param outputFormat output format, "bed" or "bedgraph" (default) or "sgr".
#' @param normalization normalization, "none" (default) or "rpkm" or multiply={1.3}).
#' @param sb path to tag2profile.pl from CTK.
#' @param perl path to PERL
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args Additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollaped <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' myrgbBed <- ctk_bed2rgb(myCollaped,col="128,0,0")
#' ctk_tag2profile(myrgbBed,verbose=TRUE)
#' @return path to BED file.
#' @export
ctk_tag2profile <- function(filesToRun,
                            outFile=paste0(file_path_sans_ext(filesToRun),".out"),
                            outFile2=NULL,
                            sb="tag2profile.pl",
                            perl="perl",
                            PATHTOPERLLIB=NULL,
                            bigFile=FALSE,
                            weight=FALSE,
                            weightAvg=FALSE,
                            ss=TRUE,
                            exact=TRUE,
                            nz=FALSE,
                            ext5=NULL,
                            ext3=NULL,
                            chromLen=NULL,
                            region=NULL,
                            minBlockSize=2000000,
                            windowSize=100,
                            stepSize=20,
                            outputFormat="bedgraph",
                            normalization="none",
                            stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2profile_stderr.txt")),
                            stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2profile_stdout.txt")),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additional_Args=NULL,verbose=FALSE, writelog = T){
  
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1] 
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  
  args <- c(cmd,
            ifelse(bigFile,"-big ",""),
            ifelse(weight,"-weight ",""),
            ifelse(weightAvg,"-weight-avg ",""),
            ifelse(ss,"-ss ",""),
            ifelse(exact,"-exact ",""),
            ifelse(nz,"-nz ",""),
            ifelse(!is.null(ext5),paste0("-ext5 ",ext5),""),
            ifelse(!is.null(ext3),paste0("-ext3 ",ext3),""),
            ifelse(!is.null(chromLen),paste0("--chrLen ",chromLen),""),
            ifelse(!is.null(region),paste0("--region ",region),""),
            paste0("-minBlockSize ",as.integer(minBlockSize)),
            paste0("-w ",windowSize),
            paste0("-s ",stepSize),
            paste0("-of ",outputFormat),
            paste0("-normalize ",normalization),
            fileToRun,
            outFile,
            ifelse(!is.null(outFile2),outFile2,""))
  
  args <- args[args!=""]
  if(verbose){
    
    message("tag2profile.pl command is ",cmd)
    message("tag2profile.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )}  else {
        system2(perl,
                args,
                stdout="",
                stderr="")
      }
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}



#' Wrapper function for ctk's bed2rgb
#'
#' Wrapper function for ctk's bed2rgb
#'
#'
#' @docType methods
#' @name ctk_tag2peak
#' @rdname ctk_tag2peak
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun path to file to process (BED).
#' @param outFile path to output file (BED).
#' @param sb path to tag2peak.pl from CTK.
#' @param perl path to PERL.
#' @param PATHTOPERLLIB path to PERL5LIB.
#' @param outBoundary output cluster boundaries.
#' @param outHalfPH output half peak height boundaries.
#' @param bigFile big input file, TRUE or FALSE (default).
#' @param ss separate the two strands, TRUE (default) or FALSE.
#' @param valleySeeking find candidate peaks by valley seeking.
#' @param valleyDepth depth of valley if valley seeking (between 0.5 and 1, default is 0.9).
#' @param genes custom gene bed file for scan statistics (will override --dbkey)
#' @param multiTest do Bonferroni multiple test correction.
#' @param useExpr use expression levels given in the score column in the gene bed file for normalization.
#' @param skipOutOfRangePeaks Remove out of bounds ranges.
#' @param pCutOff threshold of p-value to call peak (e.g. 0.01).
#' @param minPH min peak height.
#' @param maxPH max peak height.
#' @param gap merge cluster peaks closer than the gap (-1, no merge if < 0).
#' @param peakPrefix prefix of peak id (Peak) (so output file will look like Peak1, Peak2, etc).
#' @param stderr path to stdout file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @examples
#' testFasta <- system.file("extdata/hg19Small.fa",package="CLIPflexR")
#' myIndex <-suppressWarnings(bowtie2_index(testFasta, overwrite = TRUE))
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter="mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' bam <- suppressWarnings(bowtie_align(FqFile_QFColStripped,myIndex, overwrite=TRUE, inputFormat="fastq"))
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollaped <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' myrgbBed <- ctk_bed2rgb(myCollaped,col="128,0,0")
#' ctk_tag2profile(myrgbBed,verbose=TRUE)
#' @return path to BED file.
#' @export
ctk_tag2peak <- function(filesToRun,
                         outFile=paste0(file_path_sans_ext(filesToRun),".peak.bed"),
                         outBoundary=paste0(file_path_sans_ext(filesToRun),".boundary.bed"),
                         outHalfPH=paste0(file_path_sans_ext(filesToRun),".halfPF.bed"),
                         sb="tag2peak.pl",
                         perl="perl",
                         PATHTOPERLLIB=NULL,
                         bigFile=FALSE,
                         ss=TRUE,
                         valleySeeking=TRUE,
                         valleyDepth=0.9,
                         genes=NULL,
                         multiTest=FALSE,
                         useExpr=FALSE,
                         skipOutOfRangePeaks=FALSE,
                         pCutOff=0.01,
                         minPH=2,
                         maxPH=-1,
                         gap=-1,
                         peakPrefix="Peak",
                         stderr=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2peak_stderr.txt")),
                         stdout=file.path(dirname(fileToRun),paste0(basename(fileToRun),"_ctk_tag2peak_stdout.txt")),
                         useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                         additional_Args=NULL,verbose=FALSE, writelog = T){
  
  
  pathOld <- Sys.getenv("PATH",unset = NA)
  perl5libPathOld <- Sys.getenv("PERL5LIB",unset=NA)
  cmd <- sb
  
  if(useClipRConda){
    path <- paste(file.path(getOption("CLIPflexR.condaEnv"),"bin"),
                  ifelse(!is.na(pathOld),pathOld,""), 
                  collapse = .Platform$path.sep,sep=.Platform$path.sep)
    Sys.setenv("PATH"= path)
    perl <- file.path(getOption("CLIPflexR.condaEnv"),"bin",perl)
  }
  
  if(!is.null(getOption("CLIPflexR.ctk")) & useClipRConda) cmd <- paste(file.path(getOption("CLIPflexR.ctk"),cmd),sep = " ")
  if(is.null(PATHTOPERLLIB) & !is.null(getOption("CLIPflexR.czplib")) & useClipRConda) PATHTOPERLLIB <- getOption("CLIPflexR.czplib")
  
  fileToRun <- filesToRun[1] 
  
  if(!file.exists(fileToRun)) stop("File does not exist")
  
  
  # exportPATH <- ifelse(!is.null(PATHTOPERLLIB),paste0("export PERL5LIB=",PATHTOPERLLIB,";"),"")
  on.exit((function(perl5libPathOld,pathOld) {
    if(!is.na(perl5libPathOld)){
      Sys.setenv("PERL5LIB"=perl5libPathOld)
    }else{
      Sys.unsetenv("PERL5LIB")
    }
    if(!is.na(pathOld)){
      Sys.setenv("PATH"=pathOld)
    }else{
      Sys.unsetenv("PATH")
    }
  })(perl5libPathOld,pathOld))
  
  if(!is.null(PATHTOPERLLIB)){
    Sys.setenv("PERL5LIB"=PATHTOPERLLIB)
  }
  
  
  args <- c(cmd,
            ifelse(bigFile,"-big ",""),
            ifelse(ss,"-ss ",""),
            ifelse(valleySeeking,"--valley-seeking ",""),
            ifelse(!is.null(valleyDepth),paste0("--valley-depth ",valleyDepth),""),
            ifelse(!is.null(outBoundary),paste0("--out-boundary ",outBoundary),""),
            ifelse(!is.null(outHalfPH),paste0("--out-half-PH ",outHalfPH),""),
            ifelse(!is.null(genes),paste0("--gene ",genes),""),
            ifelse(multiTest,"--multi-test ",""),
            ifelse(useExpr,"--use-expr ",""),
            ifelse(skipOutOfRangePeaks,"--skip-out-of-range-peaks ",""),
            paste0("-p ",pCutOff),
            paste0("-minPH ",minPH),
            paste0("-maxPH ",maxPH),
            paste0("-gap ",gap),
            paste0("--prefix ",peakPrefix),
            fileToRun,
            outFile)
  
  
  
  args <- args[args!=""]
  if(verbose){
    
    message("tag2peak.pl command is ",cmd)
    message("tag2peak.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    if(writelog){
      system2(perl,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(perl,
                args,
                stdout="",
                stderr="")
      }
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}

