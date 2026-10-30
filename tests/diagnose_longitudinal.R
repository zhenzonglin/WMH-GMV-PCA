# Read a saved synthetic request and preserve optimizer evidence without relaxing checks.
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
req<-jsonlite::read_json(args[1],simplifyVector=FALSE)
.libPaths(c(req$config$paths$r_library,.libPaths()))
source(file.path(req$config$project$project_root,"R","clinical.R"))
req$output_dir<-args[2];dir.create(req$output_dir,recursive=TRUE,showWarnings=FALSE)
original_check<-clin_ready
check_index<-0L
clin_ready<-function(fit,gradient_tolerance) {
  check_index<<-check_index+1L
  saveRDS(fit,file.path(req$output_dir,paste0("diagnostic_fit_",check_index,".rds")))
  eig<-eigen(fit$Hessian,symmetric=TRUE,only.values=TRUE)$values
  cat("CHECK",check_index,"randomSD",as.numeric(unlist(fit$ST))[1],"eigenrange",range(eig),"maxgradient",max(abs(fit$gradient)),"\n")
  print(fit$optRes)
  original_check(fit,gradient_tolerance)
}
run_clinical(req)
