# 公共I/O；所有标识符先按字符串读入，科学参数只来自配置。
wg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))
}

wg_cfg <- function(request, ...) {
  keys <- c(...)
  value <- request$config$analysis
  for (key in keys) {
    if (is.null(value[[key]])) stop("Missing analysis configuration: ", paste(keys, collapse = "."))
    value <- value[[key]]
  }
  value
}

wg_read <- function(path) {
  if (is.null(path) || !file.exists(path)) stop("Missing table: ", path)
  out <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = "character", na.strings = c("", "NA"), quote = "")
  if (!"subject_key" %in% names(out)) stop("Table lacks subject_key: ", path)
  if (anyNA(out$subject_key) || anyDuplicated(out$subject_key)) stop("Missing or duplicate subject_key")
  out
}

wg_num <- function(x) suppressWarnings(as.numeric(x))

wg_ids <- function(data) data[, intersect(c("subject_key", "participant_id", "cohort"), names(data)), drop = FALSE]

wg_write <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

wg_json <- function(value, path) {
  wg_require("jsonlite")
  # 非有限数不能进入机器可读JSON。
  clean <- function(x) {
    if (is.list(x)) return(lapply(x, clean))
    if (is.numeric(x)) x[!is.finite(x)] <- NA_real_
    x
  }
  jsonlite::write_json(clean(value), path, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = 16)
}

wg_result <- function(request, status, reason = "", artifacts = character(), n_subjects = NULL, metrics = list()) {
  dir.create(request$output_dir, recursive = TRUE, showWarnings = FALSE)
  result <- list(status = status, reason = reason, artifacts = as.list(artifacts),
                 n_subjects = n_subjects, metrics = metrics)
  wg_json(result, file.path(request$output_dir, "result.json"))
  invisible(result)
}

wg_seed <- function(request, offset = 0L) {
  seed <- as.integer(wg_cfg(request, "seed"))
  set.seed(as.integer((as.double(seed) + as.double(offset)) %% .Machine$integer.max))
}

wg_assert_cohort <- function(data, cohort) {
  if ("cohort" %in% names(data) && any(is.na(data$cohort) | data$cohort != cohort)) {
    stop("Training table must contain only cohort ", cohort)
  }
}

wg_versions <- function(packages) setNames(lapply(packages, function(p) as.character(utils::packageVersion(p))), packages)

wg_features <- function(request) {
  features <- as.character(unlist(request$features, use.names = FALSE))
  if (anyNA(features) || anyDuplicated(features)) stop("Invalid or repeated feature names")
  features
}

wg_matrix <- function(data, features) {
  if (!all(features %in% names(data))) stop("Missing features: ", paste(setdiff(features, names(data)), collapse = ", "))
  if (!length(features)) return(matrix(numeric(), nrow(data), 0L, dimnames = list(data$subject_key, character())))
  X <- vapply(data[features], wg_num, numeric(nrow(data)))
  X <- matrix(X, nrow = nrow(data), ncol = length(features), dimnames = list(data$subject_key, features))
  X
}

wg_complete_matrix <- function(data, features) {
  X <- wg_matrix(data, features)
  good <- if (ncol(X)) rowSums(!is.finite(X)) == 0 else rep(TRUE, nrow(X))
  list(data = data[good, , drop = FALSE], X = X[good, , drop = FALSE], valid = good)
}

wg_prefix <- function(modality) {
  if (!modality %in% c("wmh", "gmv", "acute")) stop("Unknown modality: ", modality)
  toupper(modality)
}

wg_pc_names <- function(modality, k) if (k > 0L) sprintf("%s_PC%02d", wg_prefix(modality), seq_len(k)) else character()

wg_png <- function(request, filename, draw) {
  cfg <- wg_cfg(request, "figures")
  grDevices::png(file.path(request$output_dir, filename), width = as.integer(cfg$width_px),
                 height = as.integer(cfg$height_px), res = as.integer(cfg$resolution))
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  invisible(filename)
}
