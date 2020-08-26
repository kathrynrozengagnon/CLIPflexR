#' Align to Genome
#'
#' @rdname ClIP_align
#' @param samID sample id in the sample sheet
#' @param samSheet sample sheet
#' @param res_dir result directory
#' @param genome_idx file path of indexed genome
#' @param aligner Assign aligner bowtie2/bwa_mem/bwa_alb/subread/subjunc/hisat2
#' @param thread thread used in the processing, default is 4
#' @param bwa exercutable file path of bwa
#' @import Rbowtie2 Rsamtools Rsubread ShortRead 
#' 
#' @docType methods
#' @return sorted bam file
#' @export
ClIP_align <- function(samID,res_dir=NULL,genome_idx=NULL,aligner=NULL,samSheet=NULL,thread=4,bwa=NULL){
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