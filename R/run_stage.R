args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Expected one request.json")
request <- jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
root <- request$config$project$project_root
if(!is.null(request$config$paths$r_library)) .libPaths(c(request$config$paths$r_library,.libPaths()))
source(file.path(root,"R","common.R"))
source(file.path(root,"R","normative.R"))
source(file.path(root,"R","pca.R"))
source(file.path(root,"R","clinical.R"))
dir.create(request$output_dir,recursive=TRUE,showWarnings=FALSE)
if(startsWith(request$action,"norm_")) {
  run_normative(request)
} else if(startsWith(request$action,"clinical_")) {
  run_clinical(request)
} else {
  run_pca(request)
}
