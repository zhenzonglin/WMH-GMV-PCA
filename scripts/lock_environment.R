args <- commandArgs(trailingOnly=TRUE)
cfg <- jsonlite::fromJSON(args[[1]],simplifyVector=FALSE)
.libPaths(c(cfg$paths$r_library,.libPaths()))
options(repos=c(CRAN=cfg$runtime$environment$cran_repository))
packages <- unlist(cfg$runtime$environment$r_packages)
deps <- tools::package_dependencies(packages,db=installed.packages(),which=c("Depends","Imports","LinkingTo"),recursive=TRUE)
packages <- sort(setdiff(unique(c(packages,unlist(deps))),"R"))
renv::snapshot(project=cfg$project$project_root,lockfile=file.path(cfg$project$project_root,"renv.lock"),
               library=.libPaths(),packages=packages,prompt=FALSE)
versions <- setNames(lapply(packages,function(p) as.character(packageVersion(p))),packages)
jsonlite::write_json(list(R=R.version.string,packages=versions),
    file.path(cfg$project$project_root,"config","verified_r_versions.json"),pretty=TRUE,auto_unbox=TRUE)
