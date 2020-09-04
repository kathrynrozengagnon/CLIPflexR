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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to stripBarcode.pl from CTK.
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param linkerlength length of barcode/linker sequences.
#' @param inputFormat Input file format (fasta|fastq)
#' @param barcodeStartWith Filter sequences based on the starting nucleotides in the barcode
#' @param barcodeEndWith Filter sequences based on the ending nucleotides in the barcode
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped,paste0(FqFile_QF,".gz"))
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_Col,linkerlength=5,inputFormat="fastq")
#' @return Path to unzipped file
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
                             stderr=paste0(getwd(),"stripBarcode_stderr"),
                             stdout=paste0(getwd(),"stripBarcode_stdout"),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additionalArgumements=NULL,verbose=FALSE){
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)
  if (grepl("fa|fasta|fastq",file_ext(outFile))) {
    outFile <- outFile
  } else { outFile <- paste0(outFile,  "fa")}
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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param mutationBedFile Mutation BED file.
#' @param outFile Output file
#' @param sb Path to fastx_barcode_splitter.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param bigFile big file
#' @param mutationSize mutation size
#' @param permutations number of iterations for permutation
#' @param trackMutationPos track mutation position relative to read start
#' @param noSparseCorrect no sparcity correction
#' @param FDR threshold of FDR
#' @param mfr threshold of m-over-k-ratio
#' @param cacheDir Name for cache directory.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' \dontrun{
#' mutations <- system.file("extdata/BrdU.Fox.pool.tag.uniq.mutation.small.txt",package="CLIPflexR")
#' delBed <- ctk_getMutationType(mutations)
#' ctk_cims("~/Downloads/uniq_tags_mutations/Fox.pool.tag.uniq.rgb.bed",delBed,verbose=TRUE)
#' }
#' @return Path to unzipped file
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
                             stderr=paste0(getwd(),"stripBarcode_stderr"),
                             stdout=paste0(getwd(),"stripBarcode_stdout"),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additionalArgumements=NULL,verbose=FALSE){
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)


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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
    # unlink(cacheDir,recursive = TRUE)
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to fastx_barcode_splitter.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param bigFile big file
#' @param pCutOff p-value cut off
#' @param multiTest perform multiple testing correction.
#' @param gap gap size
#' @param cacheDir Name for cache directory.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' FqFile_QFColStripped <- ctk_stripBarcode(FqFile_QFCollapsed)
#' @return Path to unzipped file
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
                     stderr=paste0(getwd(),"stripBarcode_stderr"),
                     stdout=paste0(getwd(),"stripBarcode_stdout"),
                     useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                     additionalArgumements=NULL,verbose=FALSE){
  
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)
  
  # PATHTOPERLLIB=NULL,
  # bigFile=FALSE,
  # pCutOff=0.01,
  # multiTest=TRUE,
  
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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to getMutationType.pl from CTK toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param mutationType mutation size
#' @param summaryStat Create summary file
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' mutations <- system.file("extdata/BrdU.Fox.pool.tag.uniq.mutation.small.txt",package="CLIPflexR")
#' ctk_getMutationType(mutations)
#' @return Path to unzipped file
#' @export
ctk_getMutationType <- function(filesToRun,
                     outFile=paste(file_path_sans_ext(fileToRun),mutationType,"bed",sep="."),
                     sb="getMutationType.pl",
                     perl="perl",
                     PATHTOPERLLIB=NULL,
                     mutationType="del",
                     summaryStat=FALSE,
                     stderr=paste0(getwd(),"stripBarcode_stderr"),
                     stdout=paste0(getwd(),"stripBarcode_stdout"),
                     useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                     additionalArgumements=NULL,verbose=FALSE){
  
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)
  
  
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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to fastx_barcode_splitter.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param fastqFormat The fastq format used.
#' @param indexPosition Position of index in read
#' @param qsFilter Quality score filter
#' @param maxN Maximum number of bases below qsFilter 
#' @param outputFormat Output format (fastq | fasta)
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' @return Path to unzipped file
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
                             stderr=paste0(getwd(),"stripBarcode_stderr"),
                             stdout=paste0(getwd(),"stripBarcode_stdout"),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additionalArgumements=NULL,verbose=FALSE){
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)
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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to fastx_barcode_splitter.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped,paste0(FqFile_QF,".gz"))
#' FqFile_Col <- ctk_fastq2collapse(FqFile_QF,verbose=TRUE)
#' @return Path to unzipped file
#' @export
ctk_fastq2collapse <- function(filesToRun,
                            outFile=file.path(dirname(fileToRun),paste("Collapsed_",basename(fileToRun),sep="")),
                            sb="fastq2collapse.pl",
                            perl="perl",
                            PATHTOPERLLIB=NULL,
                            stderr=paste0(getwd(),"stripBarcode_stderr"),
                            stdout=paste0(getwd(),"stripBarcode_stdout"),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additionalArgumements=NULL,verbose=FALSE){
  
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
  # cmd2 <- paste0(exportPATH," ",
  #                cmd," ",
  #                " -len ",linkerlength," -v ",
  #                " ",fileToRun," ",
  #                paste(fileToRun,"_rm5",sep=""))
  # #"temp.rm")
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)
  args <- c(cmd,
            fileToRun,
            gsub("\\.gz$","",outFile))
  
  args <- args[args!=""]
  if(verbose){
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to parseAlignment.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param mutationFile Mutation file path.
#' @param mapQual Minimum map quality
#' @param minLen Minimum length of read
#' @param indelToEnd Minimum distance from indel to end of read.
#' @param splitDel Whether to split reads with deletions.
#' @param indelInScore Include indels in mutation score count.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
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
#' parsedAlignment <- ctk_parseAlignment(bam)
#' @return Path to unzipped file
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
                               stderr=paste0(getwd(),"stripBarcode_stderr"),
                               stdout=paste0(getwd(),"stripBarcode_stdout"),
                               useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                               additionalArgumements=NULL,verbose=FALSE){
  # args <- c(cmd,
  #           ifelse(!is.null(mutationFile),paste0("--mutation-file ",mutationFile),""),
  #           ifelse(!is.null(mapQual),paste0("--map-qual ",mapQual),""),
  #           ifelse(!is.null(minLen),paste0("--min-len ",minLen),""),
  #           paste0("--indel-to-end ",indelToEnd),
  #           ifelse(!splitDel,paste0("--split-del "),""),
  #           ifelse(!indelInScore,paste0("--indel-in-score "),""),
  #           fileToRun,
  #           gsub("\\.gz$","",outFile))
  
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
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to tag2collapse.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param keepMaxScore keep the tag with the most weight (instead of the longest one) as representative
#' @param keepTagName do not change tag name (no extra information)
#' @param weight consider the weight of each tag
#' @param bigFile Set to TRUE when files are big.
#' @param weightInName find weight in name
#' @param randomBarcode random barcode exists, no collapse for different barcodes
#' @param seqErrorModel sequencing error model to use (alignment or em-local or em-global or fix=0.01)
#' @param outputSeqError output sequencing errors estimated by the EM algorithm
#' @param em EM threshold to infer reliability of each collapsed read (when have random linker, -1=no EM)
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
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
#' parsedAlignment <- ctk_parseAlignment(bam)
#' ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,weightInName=FALSE,verbose=TRUE)
#' @return Path to unzipped file
#' @export
ctk_tag2collapse <- function(filesToRun,
                               outFile=file.path(dirname(filesToRun),paste0("TC_",basename(filesToRun))),
                               sb="tag2collapse.pl",
                               perl="perl",
                               PATHTOPERLLIB=NULL,
                               keepMaxScore=TRUE,
                               keepTagName=TRUE,
                               weight=TRUE,
                               bigFile=TRUE,
                               weightInName=TRUE,
                               randomBarcode=TRUE,
                               seqErrorModel=NULL,
                               outputSeqError=NULL,
                               em=NULL,
                               stderr=paste0(getwd(),"stripBarcode_stderr"),
                               stdout=paste0(getwd(),"stripBarcode_stdout"),
                               useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                               additionalArgumements=NULL,verbose=FALSE){

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
            ifelse(keepMaxScore,paste0("--keep-max-score "),""),
            ifelse(keepTagName,paste0("--keep-tag-name "),""),
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
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param file1 File 1 to join.
#' @param file2 File 2 to join.
#' @param field1 Field/column for file 1 join.
#' @param field2 Field/column for file 2 join.
#' @param mode Join mode.
#' @param outFile Output file
#' @param sb Path to tag2collapse.pl from FastX toolkit
#' @param python Path to python
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @examples
#' \dontrun{
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
#' mutationFile <- system.file("extdata/Fox3_Std_small_mutation.txt",package="CLIPflexR")
#' parsedAlignment <- ctk_parseAlignment(bam,mutationFile=mutationFile)
#' uniqueTags <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' ctk_joinWrapper(mutationFile,uniqueTags,4,4,"N",verbose=TRUE)
#' }
#' @return Path to unzipped file
#' @export
ctk_joinWrapper <- function(file1,
                            file2,
                            field1,
                            field2,
                            mode,
                             outFile=file.path(dirname(file1),paste0("Unique_",basename(file1))),
                             sb="joinWrapper.py",
                             python="python",
                             stderr=paste0(getwd(),"stripBarcode_stderr"),
                             stdout=paste0(getwd(),"stripBarcode_stdout"),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additionalArgumements=NULL,verbose=FALSE){
  # args <- c(cmd,
  #           ifelse(keepMaxScore,paste0("--keep-max-score ",keepMaxScore),""),
  #           ifelse(keepTagName,paste0("--keep-tag-name ",keepTagName),""),
  #           ifelse(!is.null(outputSeqError),paste0("--output-seq-error ",outputSeqError),""),
  #           ifelse(!is.null(em),paste0("-EM ",em),paste0("-EM -1")),
  #           ifelse(!is.null(bigFile),paste0("-big  "),""),
  #           ifelse(!is.null(weight),paste0("-weight  "),""),
  #           ifelse(!is.null(randomBarcode),paste0("--random-barcode  "),""),
  #           ifelse(!is.null(weightInName),paste0("--weight-in-name  "),""),
  #           fileToRun,
  #           outFile)
  # 
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
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(python,
            args,
            stdout=stdout,
            stderr=stderr
    )
    # if(file_ext(outFile) == "gz"){
    #   R.utils::gzip(gsub("\\.gz$","",outFile),
    #                 destname=outFile)
    # }
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to tag2collapse.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param col Colour to include in BED rgb column
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
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
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollaped <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' ctk_bed2rgb(myCollaped,col="128,0,0")
#' @return Path to unzipped file
#' @export
ctk_bed2rgb <- function(filesToRun,
                             outFile=paste0(file_path_sans_ext(filesToRun),".RGB.",file_ext(filesToRun)),
                             sb="bed2rgb.pl",
                             perl="perl",
                             PATHTOPERLLIB=NULL,
                             col="blue",
                             stderr=paste0(getwd(),"stripBarcode_stderr"),
                             stdout=paste0(getwd(),"stripBarcode_stdout"),
                             useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                             additionalArgumements=NULL,verbose=FALSE){
  # args <- c(cmd,
  #           ifelse(keepMaxScore,paste0("--keep-max-score ",keepMaxScore),""),
  #           ifelse(keepTagName,paste0("--keep-tag-name ",keepTagName),""),
  #           ifelse(!is.null(outputSeqError),paste0("--output-seq-error ",outputSeqError),""),
  #           ifelse(!is.null(em),paste0("-EM ",em),paste0("-EM -1")),
  #           ifelse(!is.null(bigFile),paste0("-big  "),""),
  #           ifelse(!is.null(weight),paste0("-weight  "),""),
  #           ifelse(!is.null(randomBarcode),paste0("--random-barcode  "),""),
  #           ifelse(!is.null(weightInName),paste0("--weight-in-name  "),""),
  #           fileToRun,
  #           outFile)
  # 
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
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @name ctk_tag2profile
#' @rdname ctk_tag2profile
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param outFile2 Output 2 file
#' @param bigFile TRUE when working with a big file
#' @param weight weight counts according to the score of each tag
#' @param weightAvg weight average the score of each tag
#' @param ss separate strand
#' @param exact exact count at each nucleotide
#' @param nz don't print zeroes (works for sgr and bed)
#' @param ext5 extension of tags at the 5' end
#' @param ext3 extension of tags at the 3' end
#' @param chromLen chrom length file
#' @param region a bed file with regions to count tag numbers. If not specified, count in moving windows
#' @param minBlockSize minimum number of lines to read in each block for a big file
#' @param windowSize Window size
#' @param stepSize Step size
#' @param outputFormat output format ([bed] or bedgraph or sgr)
#' @param normalization normalization ([none] or rpkm or multiply={1.3})
#' @param sb Path to tag2collapse.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
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
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollaped <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' myrgbBed <- ctk_bed2rgb(myCollaped,col="128,0,0")
#' ctk_tag2profile(myrgbBed,verbose=TRUE)
#' @return Path to unzipped file
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
                        stderr=paste0(getwd(),"stripBarcode_stderr"),
                        stdout=paste0(getwd(),"stripBarcode_stdout"),
                        useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                        additionalArgumements=NULL,verbose=FALSE){

  
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
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param filesToRun File to process.
#' @param outFile Output file
#' @param sb Path to tag2collapse.pl from FastX toolkit
#' @param perl Path to PERL
#' @param PATHTOPERLLIB Path to PERL5LIB.
#' @param outBoundary output cluster boundaries
#' @param outHalfPH output half peak height boundaries
#' @param bigFile big input file
#' @param ss separate the two strands
#' @param valleySeeking find candidate peaks by valley seeking
#' @param valleyDepth depth of valley if valley seeking (between 0.5 and 1, default=0.9)
#' @param genes custom gene bed file for scan statistics (will override --dbkey)
#' @param multiTest do Bonferroni multiple test correction
#' @param useExpr use expression levels given in the score column in the gene bed file for normalization
#' @param skipOutOfRangePeaks Remove out of bounds ranges.
#' @param pCutOff threshold of p-value to call peak (e.g. 0.01)
#' @param minPH min peak height
#' @param maxPH max peak height
#' @param gap merge cluster peaks closer than the gap (-1, no merge if < 0)
#' @param peakPrefix prefix of peak id (Peak) (so output file will look like Peak1, Peak2, etc)
#' @param stderr Path to stdout file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
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
#' parsedAlignment <- ctk_parseAlignment(bam)
#' myCollaped <- ctk_tag2collapse(parsedAlignment,weight=FALSE,randomBarcode=FALSE,
#' weightInName=FALSE,verbose=TRUE)
#' myrgbBed <- ctk_bed2rgb(myCollaped,col="128,0,0")
#' ctk_tag2profile(myrgbBed,verbose=TRUE)
#' @return Path to unzipped file
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
                           stderr=paste0(getwd(),"stripBarcode_stderr"),
                          stdout=paste0(getwd(),"stripBarcode_stdout"),
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additionalArgumements=NULL,verbose=FALSE){

  
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
    
    message("fastq_filter.pl command is ",cmd)
    message("fastq_filter.pl arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  if(!file.exists(outFile)){
    system2(perl,
            args,
            stdout=stdout,
            stderr=stderr
    )
    if(file_ext(outFile) == "gz"){
      R.utils::gzip(gsub("\\.gz$","",outFile),
                    destname=outFile)
    }
  }
  
  return(outFile)
}

