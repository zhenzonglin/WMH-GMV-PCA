args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1L) stop("Usage: Rscript scripts/install_r.R RESOLVED_CONFIG.json")
cfg <- jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
Sys.setenv(PATH=paste(dirname(cfg$paths$python),Sys.getenv("PATH"),sep=.Platform$path.sep))
lib <- cfg$paths$r_library
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN=cfg$runtime$environment$cran_repository), timeout=600)
pkgs <- unlist(cfg$runtime$environment$r_packages)
missing <- setdiff(pkgs,rownames(installed.packages(lib.loc=lib)))
lockfile <- file.path(cfg$project$project_root,"renv.lock")
if(file.exists(lockfile)) {
  renv::restore(project=cfg$project$project_root,lockfile=lockfile,library=lib,prompt=FALSE)
} else if(length(missing)) install.packages(missing,lib=lib,Ncpus=cfg$runtime$environment$install_jobs)
bad <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(bad)) stop(paste("R dependency installation incomplete:",paste(bad,collapse=", ")))
versions <- setNames(lapply(pkgs,function(p) as.character(packageVersion(p))),pkgs)
jsonlite::write_json(list(R=R.version.string,packages=versions),file.path(lib,"installed-versions.json"),pretty=TRUE,auto_unbox=TRUE)
