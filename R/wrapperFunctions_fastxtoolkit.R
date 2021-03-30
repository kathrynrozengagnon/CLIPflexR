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
#' @param fileToBzip2 path to file to bzip.
#' @param bzip2 path to bzip2.
#' @param keep keep (don't delete) input files, TRUE (default) or FALSE.
#' @param force overwrite existing output files.
#' @param small use less memory (at most 2500k).
#' @param blockSize set block size to 100k .. 900k.
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to bzipped file.
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
                  additional_Args=NULL,
                  verbose=FALSE){
  cmd <- bzip2
  
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileToBzip2))stop("File does not exist")
  
  fileWithoutExtension <- file_path_sans_ext(fileToBzip2)
  if(file.exists(fileToBzip2) & !file.exists(fileWithoutExtension)){
    args <- c(fileToBzip2,
              ifelse(keep,"-k",""),
              ifelse(force,"-f",""),
              ifelse(small,"-s",""),
              paste0("-",blockSize),
              ifelse(!is.null(additional_Args),additional_Args,""))
    
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
  return(paste0(fileWithoutExtension, ".bz2"))
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
#' @param fileTofqf path to file to process (fastq).
#' @param outFile output file name.
#' @param fqf path to fastq_quality_filter from FastX toolkit.
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq; default) or NULL (Phred+64 fastq).
#' @param minimumQuality minimum quality score to keep (default is 20).
#' @param minimumPercentOfRead minimum percent of bases that must have [-q] quality (default is 80).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
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
                                 additional_Args=NULL,
                                 verbose=FALSE, writelog= T){
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")
  
  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){
      
      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        paste0("-q ",minimumQuality),
        paste0("-p ",minimumPercentOfRead),
        ifelse(verbose, "-v ", ""),
        paste0("-i  ",fileTofqf),
        paste0("-o ",file_fqf)) 
      args <- args[!args %in% ""]
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
#' @param fileTofqf path to file to process (fastq).
#' @param outFile output file name.
#' @param fqf path to fastq_quality_trimmer from FastX toolkit.
#' @param qualityThreshold minimum quality score to keep.
#' @param minimumLength minimum percent of bases that must have [-q] quality.
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq; default) or NULL (Phred+64 fastq).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE.
#' @return path to unzipped file
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples  
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20)
#' FqFile_QF <- fastq_quality_trimmer(FqFile_clipped)
#' @export
fastq_quality_trimmer <- function(fileTofqf,
                                  outFile=file.path(dirname(fileTofqf),paste0("QT_",basename(fileTofqf))),
                                  fqf="fastq_quality_trimmer",
                                  qualityThreshold=5,
                                  minimumLength=20, qEncoding = 33,
                                  stderr=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_trimmer_stderr.txt")),
                                  stdout=file.path(dirname(fileTofqf),paste0(basename(fileTofqf),"_fastq_quality_trimmer_stdout.txt")),
                                  useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                  additional_Args=NULL,
                                  verbose=FALSE, writelog = T){
  cmd <- fqf
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  if(!file.exists(fileTofqf))stop("File does not exist")
  
  file_fqf <- outFile
  if(file.exists(fileTofqf) & !file.exists(file_fqf)){
      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        paste0("-t ",qualityThreshold),
        paste0("-l ",minimumLength),
        ifelse(verbose, "-v ", ""),
        paste0("-i  ",fileTofqf),
        paste0("-o ",gsub("\\.gz$","",outFile))
      )
      args <- args[!args %in% ""]
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
#' @param fileTofqs path to file to process (fastq or fasta - fasta will only return nucleotide distribution).
#' @param outFile output file name.
#' @param fqs path to fastx_quality_stats from FastX toolkit.
#' @param qEncoding quality encoding,  set to either 33 (Sanger Phred+33 encoding/Illumina fastq; default) or NULL (Phred+64 fastq or fasta).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to unzipped file.
#' @import Rsamtools Rbowtie2  GenomicAlignments
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' Fq_Stats <- fastx_quality_stats(FqFile_QF)
#' @export
fastx_quality_stats <- function(fileTofqs,
                                outFile=file.path(dirname(fileTofqs),gsub("\\.fastq|\\.fq|fq\\.gz|fastq\\.gz|\\.fasta|\\.fa|fa\\.gz|fasta\\.gz",".txt",basename(fileTofqs))),
                                fqs="fastx_quality_stats",qEncoding=33,
                                stderr=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_quality_stats_stderr.txt")),
                                stdout=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_quality_stats_stdout.txt")),
                                useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                additional_Args=NULL,verbose=FALSE){
  
  
  cmd <- fqs
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqs))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofqs) & !file.exists(file_fqs)){
    
    
    args <- c(
      ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
      paste0("-i  ",fileTofqs),
      paste0("-o ",file_fqs)
    )
    args <- args[!args %in% ""]
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
#' @param fileTofxc path to file to process (fastq or fasta).
#' @param outFile output file name (will be formatted as a fasta).
#' @param fxc path to fastx_collapser from FastX toolkit.
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq; default) or NULL (Phred+64 fastq or fasta)
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE.
#' @return path to unzipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile <- decompress(testFQ,overwrite=TRUE)
#' FqFile_QF <- fastq_quality_filter(FqFile)
#' FqFile_QFCollapsed <- fastx_collapser(FqFile_QF)
#' @export
fastx_collapser <- function(fileTofxc,
                            outFile=file.path(dirname(fileTofxc),gsub("\\.fastq|\\.fq|\\.fa|\\.fasta","_collapse.fasta",basename(fileTofxc))),
                            fxc="fastx_collapser", qEncoding = NULL,
                            stderr=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_collapse_stderr.txt")),
                            stdout=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_collapse_stdout.txt")),
                            useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                            additional_Args=NULL,verbose=FALSE, writelog = T){
  
  
  cmd <- fxc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofxc))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofxc) & !file.exists(file_fqs)){

      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        ifelse(verbose, "-v ", ""),
        paste0("-i  ",fileTofxc),
        paste0("-o ",file_fqs))
        args <- args[!args %in% ""]
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
  return(file_fqs)
}



#' Wrapper function for fastx_barcode_splitter
#'
#' Split multiplexed samples by user-defined indics using fastx_barcode_splitter
#'
#'
#' @docType methods
#' @name fastx_barcode_splitter
#' @rdname fastx_barcode_splitter
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileTofxc path to file to process (fastq or fasta).
#' @param bcFile tab-delimited barcode file.
#' @param mismatches number of mismatches allowed.
#' @param fbs path to fastx_barcode_splitter.pl from FASTX toolkit.
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to split files; files will be written without any extension, in the same format as the input file. 
#' @export
fastx_barcode_splitter <- function(fileTofxc,bcFile,mismatches=0,
                                   fbs="fastx_barcode_splitter.pl",
                                   stderr=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_barcode_splitter_stderr.txt")),
                                   stdout=file.path(dirname(fileTofxc),paste0(basename(fileTofxc),"_fastx_barcode_splitter_stdout.txt")),
                                   useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                                   additional_Args=NULL,verbose=FALSE){
  
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
  output3 <- as.character(output2$V3) 
  return(output3)
  if(verbose){
    print(samplestats)   
  }
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
#' @param fileTofqs path to file to process (fastq or fasta).
#' @param outFile output file path.
#' @param fqc path to fastx_clipper from FastX toolkit.
#' @param length miminum sequence length, default is 18 (set to .
#' @param adapter adapter to remove (default is "GTGTCAGTCACTTCCAGCGG", specify a string to change adapter sequence).
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq) or NULL (Phred+64 fastq or fasta; default).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to clipped file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_clipped <- fastx_clipper(FqFile,length=20, adapter = "GTGTCAG")
#' @export
fastx_clipper <- function(fileTofqs,
                          outFile=paste0(file_path_sans_ext(fileTofqs),"_clip.",file_ext(fileTofqs)),
                          fqc="fastx_clipper",length=18,
                          adapter="GTGTCAGTCACTTCCAGCGG", qEncoding = NULL, writelog = T,
                          stderr=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_clipper_stderr.txt")),
                          stdout=file.path(dirname(fileTofqs),paste0(basename(fileTofqs),"_fastx_clipper_stdout.txt")),
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additional_Args=NULL,verbose=FALSE){
  
  
  cmd <- fqc
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqs))stop("File does not exist")
  
  file_fqs <- outFile
  if(file.exists(fileTofqs) & !file.exists(file_fqs)){
      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        paste0("-l ",length),
        paste0("-a  ",adapter),
        ifelse(verbose, "-v ", ""),
        paste0("-o ",file_fqs),
        paste0("-i ",fileTofqs))
      args <- args[!args %in% ""]
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
#' @param fileTofqt path to file to process (fastq or fasta).
#' @param outFile output file path.
#' @param fqt path to fastx_trimmer from FastX toolkit.
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @param read_start read starting base (default is 10).
#' @param read_end read ending base (default is NULL, set to integer to trim).
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq) or NULL (Phred+64 fastq or fasta; default).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to trimmed file.
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_trimmed <- fastx_trimmer(FqFile,read_start = 5, read_end = 50)
#' @export
fastx_trimmer <- function(fileTofqt,fqt="fastx_trimmer",read_start = 10, read_end = NULL, qEncoding = NULL, 
                          outFile=paste0(file_path_sans_ext(fileTofqt),"_trim.",file_ext(fileTofqt)), writelog = T,
                          stderr=file.path(dirname(fileTofqt),paste0(basename(fileTofqt),"_trimmer_stats_stderr.txt")),
                          stdout=file.path(dirname(fileTofqt),paste0(basename(fileTofqt),"_trimmer_stats_stdout.txt")),  
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additional_Args=NULL,verbose=FALSE){
  cmd <- fqt
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileTofqt))stop("File does not exist")
  if (grepl("fa|fasta|fastq|fq",file_ext(outFile))) {
    outFile <- outFile
  } else { outFile <- paste0(outFile,  "fa")}
  
  file_fqt <- outFile
  if(file.exists(fileTofqt) & !file.exists(file_fqt)){
      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        paste0("-f ",read_start),
        ifelse(verbose, "-v ", ""),
        ifelse(!is.null(read_end), paste0("-l ",read_end),""),
        paste0("-o ",file_fqt),
        paste0("-i ",fileTofqt)) 
      args <- args[!args %in% ""]
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
  return(file_fqt)
}


#' Wrapper function for fastq_to_fasta
#'
#' Wrapper function for fastq_to_fasta
#'
#'
#' @docType methods
#' @name fastq_to_fasta
#' @rdname fastq_to_fasta
#'
#' @author Kathryn Rozen-Gagnon
#'
#' @param fileToc path to file to process (fastq).
#' @param outFile output file path.
#' @param fqa path to fastx_trimmer from FastX toolkit.
#' @param writelog write stderr/stdout logs, TRUE (default) or FALSE. 
#' @param discard_N discard reads with unknown nucleotides (N), TRUE (default)  or FALSE.
#' @param rename replace sequence names with number identifiers, TRUE or FALSE (default).
#' @param qEncoding quality encoding, set to either 33 (Sanger Phred+33 encoding/Illumina fastq; default) or NULL (Phred+64 fastq).
#' @param stderr path to stderr file.
#' @param stdout path to stdout file.
#' @param useClipRConda use conda environment installed by Herper, TRUE (default) or FALSE.
#' @param additional_Args additional arguments to be passed to system call.
#' @param verbose print messages to screen, TRUE or FALSE (default).
#' @return path to trimmed file
#' @examples
#' testFQ <- system.file("extdata/Fox3_Std_small.fq.gz",package="CLIPflexR")
#' FqFile_FF <- ctk_fastqFilter(testFQ,qsFilter = "mean:0-29:20",verbose=TRUE)
#' FqFile <- decompress(FqFile_FF,overwrite=TRUE)
#' FqFile_fa <- fastx_qtoa(FqFile)
#' @export
fastx_qtoa <- function(fileToc,fqa="fastq_to_fasta",discard_N = T, rename = F,
                          outFile=paste0(file_path_sans_ext(fileToc),".fa"), writelog = T, qEncoding = 33,
                          stderr=file.path(dirname(fileToc),paste0(basename(fileToc),"_qtoa_stderr.txt")),
                          stdout=file.path(dirname(fileToc),paste0(basename(fileToc),"_qtoa_stdout.txt")),  
                          useClipRConda=ifelse(is.null(getOption("CLIPflexR.condaEnv")),FALSE,TRUE),
                          additional_Args=NULL,verbose=FALSE){
  cmd <- fqa
  if(useClipRConda) cmd <- file.path(getOption("CLIPflexR.condaEnv"),"bin",cmd)
  
  if(!file.exists(fileToc))stop("File does not exist")
  if (!grepl("fastq|fq",file_ext(fileToc))) {
    message("Warning: make sure file is in fastq format")
  } 
   file_fqc <- outFile
  if(file.exists(fileToc) & !file.exists(file_fqc)){
      args <- c(
        ifelse(!is.null(qEncoding), paste0("-Q ",qEncoding),""),
        ifelse(verbose, paste0("-v "), ""),
        ifelse(!discard_N, paste0("-n "), ""),
        ifelse(rename, paste0("-r "), ""),
        paste0("-i ", fileToc),paste0("-o ", file_fqc))
      args <- args[!args %in% ""]
      }

    if(verbose){      
      message("fastx_qtoa command is ",cmd)
      message("fastx_qtoa arguments are ",paste0(args,sep=" ",collapse=" "))
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
    
  return(file_fqc)
}
