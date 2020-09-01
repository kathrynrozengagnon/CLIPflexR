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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs.
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file.
#' @import Rsamtools GenomicAlignments
#' @export
bzip2 <- function(fileToBzip2,bzip2="Bzip2",
                  keep=TRUE,
                  force=FALSE,
                  small=FALSE,
                  blockSize=1,
                  stderr=paste0(getwd(),"gunzip_stderr"),
                  stdout=paste0(getwd(),"gunzip_stdout"),
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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
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
                                 stderr=paste0(getwd(),"fastq_quality_filter_stderr"),
                                 stdout=paste0(getwd(),"fastq_quality_filter_stdout"),
                                 useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                 additionalArgumements=NULL,
                                 verbose=FALSE){

  # [-h]         = This helpful help screen.
  # [-q N]       = Minimum quality score to keep.
  # [-p N]       = Minimum percent of bases that must have [-q] quality.
  # [-z]         = Compress output with GZIP.
  # [-i INFILE]  = FASTA/Q input file. default is STDIN.
  # [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
  # [-v]         = Verbose - report number of sequences.
  # If [-o] is specified,  report will be printed to STDOUT.
  # If [-o] is not specified (and outp
                            
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")

  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){


    args <- c(
      paste0("-Q ",qEncoding),
      paste0("-q ",minimumQuality),
      paste0("-p ",minimumPercentOfRead),
      paste0("-i  ",fileTofqf),
      paste0("-o ",file_fqf)
    )
    if(verbose){      
      message("fastq_quality_filter command is ",cmd)
      message("fastq_quality_filter arguments are ",paste0(args,sep=" ",collapse=" "))
    }

    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
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
                                 stderr=paste0(getwd(),"fastq_quality_filter_stderr"),
                                 stdout=paste0(getwd(),"fastq_quality_filter_stdout"),
                                 useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                 additionalArgumements=NULL,
                                 verbose=FALSE){
  
  # [-h]         = This helpful help screen.
  # [-q N]       = Minimum quality score to keep.
  # [-p N]       = Minimum percent of bases that must have [-q] quality.
  # [-z]         = Compress output with GZIP.
  # [-i INFILE]  = FASTA/Q input file. default is STDIN.
  # [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
  # [-v]         = Verbose - report number of sequences.
  # If [-o] is specified,  report will be printed to STDOUT.
  # If [-o] is not specified (and outp
  
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")
  
  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){
    
    
    args <- c(
      paste0("-t ",qualityThreshold),
      paste0("-l ",minimumLength),
      paste0("-i  ",fileTofqf),
      paste0("-o ",gsub("\\.gz$","",outFile))
    )
    if(verbose){      
      message("fastq_quality_trimmer command is ",cmd)
      message("fastq_quality_trimmer arguments are ",paste0(args,sep=" ",collapse=" "))
    }
    
    system2(cmd,
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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' Fq_Stats <- fastx_quality_stats(FqFile_QF)
#' @export
fastx_quality_stats <- function(fileTofqs,
                                outFile=file.path(dirname(fileTofqs),gsub("\\.fastq|fq|fq\\.gz|fastq\\.gz",".txt",basename(fileTofqs))),
                                fqs="fastx_quality_stats",qEncoding=33,
                                stderr=paste0(getwd(),"fastq_quality_stats_stderr"),
                                stdout=paste0(getwd(),"fastq_quality_stats_stdout"),
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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' @export
fastx_collapser <- function(fileTofxc,
                            outFile=file.path(dirname(fileTofxc),gsub("\\.fastq|\\.fq","_collapse.fastq",basename(fileTofxc))),
                            fxc="fastx_collapser",qEncoding=33,
                            stderr=paste0(getwd(),"fastq_collapse_stderr"),
                            stdout=paste0(getwd(),"fastq_collapse_stdout"),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additionalArgumements=NULL,verbose=FALSE){
  
  
  cmd <- fxc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofxc))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofxc) & !file.exists(file_fqs)){
    
    
    args <- c(
      paste0("-Q ",qEncoding),
      paste0("-i  ",fileTofxc),
      paste0("-o ",file_fqs)
    )
    if(verbose){      
      message("fastx_collapser command is ",cmd)
      message("fastx_collapser arguments are ",paste0(args,sep=" ",collapse=" "))
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
#' @name fastx_barcode_splitter
#' @rdname fastx_barcode_splitter
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofxc File to process.
#' @param bcFile Barcode file
#' @param mismatches Number of mismatches allowed.
#' @param fbs Path to fastx_barcode_splitter.pl from FastX toolkit
#' @param stderr Path to stderr file.
#' @param stdout Path to stdout file.
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @export
fastx_barcode_splitter <- function(fileTofxc,bcFile,mismatches=0,
                                   fbs="fastx_barcode_splitter.pl",
                                   stderr=paste0(getwd(),"fastx_barcode_splitter_stderr"),
                                   stdout=paste0(getwd(),"fastx_barcode_splitter_stdout"),
                                   useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                   additionalArgumements=NULL,verbose=FALSE){
                                   


  cmd <- fbs
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)

  if(!file.exists(fileTofxc))stop("File does not exist")

  prefix <- gsub("QF_|\\.fasta|\\.fastq","",basename(fileTofxc))

  args <- c(
    paste0("--bcfile ",bcFile),
    "--bol",
    paste0("--mismatches ",mismatches),
    paste0("--prefix '",prefix,"_' ")
  )
  if(verbose){      
    message("fastx_collapser command is ",cmd)
    message("fastx_collapser arguments are ",paste0(args,sep=" ",collapse=" "))
  }
  
  system2(cmd,
          args,
          stdout=stdout,
          stderr=stderr
  )
  
  # cmd2 <- paste0(cmd," ",
  #                " --bcfile ",bcFile," ",
  #                "--bol --mismatches ",mismatches," ",
  #                "--prefix '",prefix,"_' ")
  
  # temp <- system(cmd2,wait = TRUE,intern = TRUE)

  return(NULL)
}

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
#' @param useClipRConda Boolean on whether to use conda environment install by CondaSysReqs
#' @param additionalArgumements Additional arguments to be passed to system call.
#' @param verbose Print more message to screen.
#' @return Path to unzipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="clipR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' @export
fastx_clipper <- function(fileTofqs,
                          outFile=paste0(file_path_sans_ext(fileTofqs),"_clip.",file_ext(fileTofqs)),
                          fqc="fastx_clipper",length=18,
                          adaptor="GTGTCAGTCACTTCCAGCGG",
                          stderr=paste0(getwd(),"clipper_stats_stderr"),
                          stdout=paste0(getwd(),"clipper_stats_stdout"),
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additionalArgumements=NULL,verbose=FALSE){


  cmd <- fqc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)

  if(!file.exists(fileTofqs))stop("File does not exist")

  file_fqs <- outFile
  if(file.exists(fileTofqs) & !file.exists(file_fqs)){

    args <- c(
      paste0("-l ",length),
      paste0("-a  ",adaptor),
      paste0("-o ",file_fqs),
      paste0("-i ",fileTofqs)
    )
    if(verbose){      
      message("fastx_clipper command is ",cmd)
      message("fastx_clipper arguments are ",paste0(args,sep=" ",collapse=" "))
    }

    system2(cmd,
            args,
            stdout=stdout,
            stderr=stderr
    )
  }
  return(file_fqs)
}


