#' Build index for referene genome
#'
#' @rdname ClIP_biuldIDX
#' @param ref_seq reference genome
#' @param res_dir result direcotry
#' @param preFIX file name prefix for indexed genome
#' @param aligner Assign aligner bowtie2/subread/bwa/hisat2
#' @param thread threads used in processing, default is 4
#' @param bwa exercutable bwa file path
#' @import Rbowtie2 Rsubread Rhisat2
#' 
#' @docType methods
#' @return indexed reference genome
#' @export
ClIP_buildIDX <- function(ref_seq=NULL,res_dir=NULL,preFIX=NULL,aligner=NULL,thread=4,bwa=NULL){
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