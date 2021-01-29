#' Wrapper function for bzip2
#'
#' Wrapper function for bzip2
#'
#'
#' @docType methods
#' @name bzip2
#' @rdname bzip2
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToBzip2 File to bzip.
#' @param bzip2 Path to bzip2.
#' @param keep keep (don't delete) input files.
#' @param force overwrite existing output files.
#' @param small use less memory (at most 2500k).
#' @param blockSize  set block size to 100k .. 900k.
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper.
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @return Path to unzipped file.
#' @import Rsamtools GenomicAlignments
#' @export
bzip2 <- function(fileToBzip2,bzip2="Bzip2",
                  keep=TRUE,
                  force=FALSE,
                  small=FALSE,
                  blockSize=1,
                  stderr=file.path(dirname(fileToBzip2),paste0(basename(fileToBzip2),"_Bzip2_stderr.txt")),
                  stdout=file.path(dirname(fileToBzip2),paste0(basename(fileToBzip2),"_Bzip2_stdout.txt")),
                  useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                  additionalArgumements=NULL,
                  verbose=FALSE){
  cmd <- bzip2
  
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileToBzip2))stop("File does not exist")
  
  fileWithoutExtension <- file_path_sans_ext(fileToBzip2)
  if(file.exists(fileToBzip2) & !file.exists(fileWithoutExtension)){
    args <- c(fileWithoutExtension,
              ifelse(keep,"-k",""),
              ifelse(force,"-f",""),
              ifelse(small,"-s",""),
              paste0("-",blockSize),
              ifelse(!is.null(additionalArgumements),additionalArgumements,""))
    
    if(verbose){      
      message("Bzip2 command is ",cmd)
      message("Bzip2 arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    
    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
  }
  return(fileWithoutExtension)
}


#' Wrapper function for fastq_quality_filter
#'
#' Wrapper function for fastq_quality_filter
#'
#'
#' @docType methods
#' @name fastq_quality_filter
#' @rdname fastq_quality_filter
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqf File to process.
#' @param outFile Output file name.
#' @param fqf Path to fastq_quality_filter from FastX toolkit
#' @param qEncoding Quality encoding
#' @param minimumQuality Minimum quality score to keep.
#' @param minimumPercentOfRead Minimum percent of bases that must have [-q] quality.
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' @export
fastq_quality_filter <- function(fileTofqf,
                                 outFile=file.path(dirname(fileTofqf),paste0("QF_",basename(fileTofqf))),
                                 fqf="fastq_quality_filter",qEncoding=33,
                                 minimumQuality=20,
                                 minimumPercentOfRead=80,
                                 stderr=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_filter_stderr.txt")),
                                 stdout=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_filter_stdout.txt")),
                                 useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                 additionalArgumements=NULL,
                                 verbose=FALSE, writelog= T){
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")
  
  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){
    if(verbose){
      
      args <- c(
        paste0("-Q ",qEncoding),
        paste0("-q ",minimumQuality),
        paste0("-p ",minimumPercentOfRead),
        paste0("-v "),
        paste0("-i  ",fileTofqf),
        paste0("-o ",file_fqf)) } else {
          args <- c(
            paste0("-Q ",qEncoding),
            paste0("-q ",minimumQuality),
            paste0("-p ",minimumPercentOfRead),
            paste0("-i  ",fileTofqf),
            paste0("-o ",file_fqf) )
        }
    if(verbose){      
      message("fastq_quality_filter command is ",cmd)
      message("fastq_quality_filter arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    if(writelog){
      system2(cmd,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(cmd,
                args,
                stdout="",
                stderr="") 
      }
  }
  return(file_fqf)
}


#' Wrapper function for fastq_quality_trimmer
#'
#' Wrapper function for fastq_quality_trimmer
#'
#'
#' @docType methods
#' @name fastq_quality_trimmer
#' @rdname fastq_quality_trimmer
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqf File to process.
#' @param outFile Output file name.
#' @param fqf Path to fastq_quality_trimmer from FastX toolkit
#' @param qualityThreshold Minimum quality score to keep.
#' @param minimumLength Minimum percent of bases that must have [-q] quality.
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' @export
fastq_quality_trimmer <- function(fileTofqf,
                                  outFile=file.path(dirname(fileTofqf),paste0("QF_",basename(fileTofqf))),
                                  fqf="fastq_quality_trimmer",
                                  qualityThreshold=5,
                                  minimumLength=20,
                                  stderr=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_trimmer_stderr.txt")),
                                  stdout=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_trimmer_stdout.txt")),
                                  useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                  additionalArgumements=NULL,
                                  verbose=FALSE, writelog = T){
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")
  
  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){
    if (verbose){
      args <- c(
        paste0("-t ",qualityThreshold),
        paste0("-l ",minimumLength),
        paste0("-v "),
        paste0("-i  ",fileTofqf),
        paste0("-o ",gsub("\\.gz$","",outFile))
      )
      
    }else {
      args <- c(
        paste0("-t ",qualityThreshold),
        paste0("-l ",minimumLength),
        paste0("-i  ",fileTofqf),
        paste0("-o ",gsub("\\.gz$","",outFile))
      ) 
    }
    
    
    if(verbose){      
      message("fastq_quality_trimmer command is ",cmd)
      message("fastq_quality_trimmer arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    if(writelog){
      system2(cmd,
              args,
              stdout=stdout,
              stderr=stderr
      )} else {
        system2(cmd,
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

#' Wrapper function for fastx_quality_stats
#'
#' Wrapper function for fastx_quality_stats
#'
#'
#' @docType methods
#' @name fastx_quality_stats
#' @rdname fastx_quality_stats
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqs File to process.
#' @param outFile Output file name.
#' @param fqs Path to fastx_quality_stats from FastX toolkit
#' @param qEncoding Quality encoding
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' Fq_Stats <- fastx_quality_stats(FqFile_QF)
#' @export
fastx_quality_stats <- function(fileTofqs,
                                outFile=file.path(dirname(fileTofqs),gsub("\\.fastq|fq|fq\\.gz|fastq\\.gz",".txt",basename(fileTofqs))),
                                fqs="fastx_quality_stats",qEncoding=33,
                                stderr=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_quality_stats_stderr.txt")),
                                stdout=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_quality_stats_stdout.txt")),
                                useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                additionalArgumements=NULL,verbose=FALSE){
  
  
  cmd <- fqs
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqs))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofqs) & !file.exists(file_fqs)){
    
    
    args <- c(
      paste0("-Q ",qEncoding),
      paste0("-i  ",fileTofqs),
      paste0("-o ",file_fqs)
    )
    if(verbose){      
      message("fastx_quality_stats command is ",cmd)
      message("fastx_quality_stats arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    
    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
  }
  return(file_fqs)
}

#' Wrapper function for fastx_collapser
#'
#' Wrapper function for fastx_collapser
#'
#'
#' @docType methods
#' @name fastx_collapser
#' @rdname fastx_collapser
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofxc File to process.
#' @param outFile Output file name.
#' @param fxc Path to fastx_collapser from FastX toolkit
#' @param qEncoding Quality encoding
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE
#' @return Path to unzipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' @export
fastx_collapser <- function(fileTofxc,
                            outFile=file.path(dirname(fileTofxc),gsub("\\.fastq|\\.fq","_collapse.fasta",basename(fileTofxc))),
                            fxc="fastx_collapser",qEncoding=33,
                            stderr=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_collapse_stderr.txt")),
                            stdout=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_collapse_stdout.txt")),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additionalArgumements=NULL,verbose=FALSE, writelog = T){
  
  
  cmd <- fxc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofxc))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofxc) & !file.exists(file_fqs)){
    
    if(verbose) {
      args <- c(
        paste0("-Q ",qEncoding),
        paste0("-v "),
        paste0("-i  ",fileTofxc),
        paste0("-o ",file_fqs)
      )} else {
        args <- c(
          paste0("-Q ",qEncoding),
          paste0("-i  ",fileTofxc),
          paste0("-o ",file_fqs))   
      }
    if(verbose){      
      message("fastx_collapser command is ",cmd)
      message("fastx_collapser arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    if(writelog) {
      system2(cmd,
              args,
              stdout=stdout,
              stderr=stderr)} else {
                system2(cmd,
                        args,
                        stdout="",
                        stderr="")            }
  }
  return(file_fqs)
}



#' Wrapper function for fastx_barcode_splitter
#'
#' Wrapper function for fastx_barcode_splitter
#'
#'
#' @docType methods
#' @name fastx_barcode_splitter
#' @rdname fastx_barcode_splitter
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofxc File to process.
#' @param bcFile Barcode file
#' @param mismatches Number of mismatches allowed.
#' @param fbs Path to fastx_barcode_splitter.pl from FASTX toolkit
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @return output3, path to split files
#' @export
fastx_barcode_splitter <- function(fileTofxc,bcFile,mismatches=0,
                                   fbs="fastx_barcode_splitter.pl",
                                   stderr=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_barcode_splitter_stderr.txt")),
                                   stdout=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_barcode_splitter_stdout.txt")),
                                   useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                   additionalArgumements=NULL,verbose=FALSE){
  
  cmd <- fbs
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofxc))stop("File does not exist")
  
  prefix <- file.path(dirname(fileTofxc), gsub("QF_|\\.fa|\\.fasta|\\.fastq","",basename(fileTofxc)))
  
  cmd2 <- c(paste0("cat ", fileTofxc," | ", cmd, " --bcfile ",bcFile,
                   " --bol --mismatches ",mismatches, " --prefix '",prefix,"_'"))
  if(verbose){      
    message("barcode_splitter command is ",cmd2)
    message("barcode_splitter arguments are ",paste0(" --bcfile ",bcFile),
            "--bol",paste0("--mismatches ",mismatches),
            paste0(" --prefix '",prefix,"_'"))
  }
  
  output <- system(cmd2, wait = TRUE, intern = TRUE)
  samplestats <- as.data.frame(str_split_fixed(output, "\t",  3))
  write.table(samplestats, file = file.path(dirname(fileTofxc),gsub(".txt",  "_stats.txt", basename(bcFile))), col.names = F, row.names = F, sep = "\t",  quote = F)
  BCFILE  <- read.delim(bcFile, header = F, sep =  "\t")
  output2 <- samplestats[samplestats$V1 %in% BCFILE$V1,]
  output3 <- output2$V3 
  return(output3)
  if(verbose){
    print(samplestats)   
  }
}

# cmd2 <- paste0(cmd," ",
#                " --bcfile ",bcFile," ",
#                "--bol --mismatches ",mismatches," ",
#                "--prefix '",prefix,"_' ")

# temp <- system(cmd2,wait = TRUE,intern = TRUE)




#' Wrapper function for fastx_clipper
#'
#' Wrapper function for fastx_clipper
#'
#'
#' @docType methods
#' @name fastx_clipper
#' @rdname fastx_clipper
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqs File to process.
#' @param outFile Output file path
#' @param fqc Path to fastx_clipper from FastX toolkit
#' @param length Length for fastx_clipper
#' @param adaptor Adapter to remove
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print messages to screen
#' @return Path to clipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' @export
fastx_clipper <- function(fileTofqs,
                          outFile=paste0(file_path_sans_ext(fileTofqs),"_clip.",file_ext(fileTofqs)),
                          fqc="fastx_clipper",length=18,
                          adaptor="GTGTCAGTCACTTCCAGCGG", writelog = T,
                          stderr=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_clipper_stderr.txt")),
                          stdout=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_clipper_stdout.txt")),
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additionalArgumements=NULL,verbose=FALSE){
  
  
  cmd <- fqc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqs))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofqs) & !file.exists(file_fqs)){
    if(verbose) {
      args <- c(
        paste0("-l ",length),
        paste0("-a  ",adaptor),
        paste0("-v" ),
        paste0("-o ",file_fqs),
        paste0("-i ",fileTofqs))
    } else { 
      args <- c(
        paste0("-l ",length),
        paste0("-a  ",adaptor),
        paste0("-o ",file_fqs),
        paste0("-i ",fileTofqs))
    }
    if(verbose){      
      message("fastx_clipper command is ",cmd)
      message("fastx_clipper arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    if(writelog) {
      system2(cmd,
              args,
              stdout=stdout,
              stderr=stderr)}  else {
                system2(cmd,
                        args,
                        stdout="",
                        stderr="")
              }
  }
  return(file_fqs)
}

#' Wrapper function for fastx_trimmer
#'
#' Wrapper function for fastx_trimmer
#'
#'
#' @docType methods
#' @name fastx_trimmer
#' @rdname fastx_trimmer
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofqt File to process
#' @param outFile Output file path
#' @param fqt Path to fastx_trimmer from FastX toolkit
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE 
#' @param read_start read starting base
#' @param read_end read ending base
#' @param stderr Path to stderr file
#' @param stdout Path to stdout file
#' @param useClipRConda Boolean on whether to use conda environment install by Herper
#' @param additionalArgumements Additional arguments to be passed to system call
#' @param verbose Print more message to screen
#' @return Path to trimmed file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' @export
fastx_trimmer <- function(fileTofqt,fqt="fastx_trimmer",read_start = 10, read_end = NULL,
                          outFile=paste0(file_path_sans_ext(fileTofqt),"_trim.",file_ext(fileTofqt)), writelog = T,
                          stderr=file.path(dirname(fileTofqt),paste0(basename(fileTofqt),"_trimmer_stats_stderr.txt")),
                          stdout=file.path(dirname(fileTofqt),paste0(basename(fileTofqt),"_trimmer_stats_stdout.txt")),  
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additionalArgumements=NULL,verbose=FALSE){
  cmd <- fqt
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqt))stop("File does not exist")
  if (grepl("fa|fasta|fastq",file_ext(outFile))) {
    outFile <- outFile
  } else { outFile <- paste0(outFile,  "fa")}
  
  file_fqt <- outFile
  if(file.exists(fileTofqt) & !file.exists(file_fqt)){
    if(verbose) {
      args <- c(
        paste0("-f ",read_start),
        paste0("-v "),
        paste0("-o ",file_fqt),
        paste0("-i ",fileTofqt)) } else {
          args <- c(
            paste0("-f ",read_start),
            paste0("-o ",file_fqt),
            paste0("-i ",fileTofqt))
        }
    if(verbose){      
      message("fastx_trimmer command is ",cmd)
      message("fastx_trimmer arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    if(writelog) {
      system2(cmd,
              args,
              stdout=stdout,
              stderr=stderr)} else {
                system2(cmd,
                        args,
                        stdout="",
                        stderr=""  )        
              }
    
  }
  return(file_fqt)
}

