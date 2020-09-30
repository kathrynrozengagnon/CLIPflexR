# Taken from reticulate as not exported
miniconda_path <- function () {
  Sys.getenv("RETICULATE_MINICONDA_PATH", unset = miniconda_path_default())
}

# Taken from reticulate as not exported
miniconda_conda <- function (path = miniconda_path()) 
{
  exe <- if (is_windows()) 
    "condabin/conda.bat"
  else "bin/conda"
  file.path(path, exe)
}

# Taken from reticulate as not exported
is_osx <- function (){
  Sys.info()["sysname"] == "Darwin"
}

# Taken from reticulate as not exported
is_windows <- function () 
{
  identical(.Platform$OS.type, "windows")
}

# Taken from reticulate as not exported
#' @import rappdirs 
miniconda_path_default <- function () {
  if (is_osx()) 
    return(path.expand("~/Library/r-miniconda"))
  root <- normalizePath(rappdirs::user_data_dir(), winslash = "/", 
                        mustWork = FALSE)
  file.path(root, "r-miniconda")
}

# Taken from reticulate as not exported
miniconda_exists <- function (path = miniconda_path()) {
  conda <- miniconda_conda(path)
  file.exists(conda)
}

.onAttach <- function(libname, pkgname) {
  miniCondaPath <- miniconda_path()
  miniCondaPathExists <- miniconda_exists(miniCondaPath)
  clipr <- file.path(miniCondaPath,"envs",paste0(pkgname,"_",packageVersion(pkgname)))
  op <- options()
  op.CLIPflexR <- list(
    CLIPflexR.condaEnv = NULL,
    CLIPflexR.ctk = NULL,
    CLIPflexR.czplib = NULL
  )

  

  if(dir.exists(clipr)){
    packageStartupMessage(basename(clipr)," conda env found at ",clipr)
    op.CLIPflexR$CLIPflexR.condaEnv <- clipr
  }else{
    packageStartupMessage(basename(clipr)," conda env was not found found at ",clipr)
    packageStartupMessage("For more information on setting up conda environments for CLIPflexR, see the Herper package")
  }
  ctk <- file.path(clipr,"bin","ctk")
  czplib <- file.path(clipr,"lib","czplib")
  ctkInstalled <- dir.exists(ctk)
  czplibInstalled <- dir.exists(czplib)
  if(ctkInstalled & czplibInstalled){
    packageStartupMessage("ctk found  at ",ctk)
    packageStartupMessage("czplib found  at ",czplib)
    op.CLIPflexR$CLIPflexR.ctk <- ctk
    op.CLIPflexR$CLIPflexR.czplib <- czplib
  }else{
    packageStartupMessage("ctk not found at ",ctk)
    packageStartupMessage("czplib not found at ",czplib)
  }

  toset <- !(names(op.CLIPflexR) %in% names(op))
  if(any(toset)) options(op.CLIPflexR[toset])
  invisible()
}
