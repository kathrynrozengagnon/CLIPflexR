
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
#' @importFrom utils download.file packageVersion read.delim read.table unzip write.table
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