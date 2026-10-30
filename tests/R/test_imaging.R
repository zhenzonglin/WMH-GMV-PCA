# 在项目根目录运行：Rscript tests/R/test_imaging.R
# 仅合成数据；降低重复次数不等于正式参数的统计验证。
source("R/common.R")
source("R/normative.R")
source("R/pca.R")
wg_require(c("jsonlite", "gamlss", "gamlss.dist", "FactoMineR", "missMDA", "clue"))
scratch <- Sys.getenv("WMH_GMV_TEST_OUTPUT_DIR", unset = tempfile("wmh_gmv_r_tests_"))
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
config <- list(analysis = list(
  seed = 431L, age_range = c(50, 75),
  normative = list(folds = 5L, age_spline_df = 3L, min_positive = 50L, cdf_clip = 1e-6, prediction_interval = .90, max_cycles = 100L),
  pca = list(max_components = 2L, cv_repeats = 2L, mask_fraction = .05, bootstrap_repeats = 2L,
             convergence_tolerance = 1e-6, near_degenerate_ratio = 1.10,
             angle_median_max = 15, angle_p95_max = 30, loading_similarity_min = .90, score_icc_min = .80),
  acute = list(max_components = 2L, volume_spline_df = 3L),
  figures = list(width_px = 1000L, height_px = 700L, resolution = 100L, max_heatmap_features = 40L)
))
request <- function(action, tag, inputs, ...) {
  out <- file.path(scratch, tag)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  c(list(action = action, output_dir = out, inputs = inputs, config = config), list(...))
}
read_result <- function(req) jsonlite::read_json(file.path(req$output_dir, "result.json"), simplifyVector = TRUE)
read_scores <- function(req) wg_read(file.path(req$output_dir, "scores.tsv"))

set.seed(11)
n <- 240L
ls <- data.frame(subject_key = sprintf("LS:%05d", seq_len(n)), participant_id = sprintf("%05d", seq_len(n)),
                 cohort = "LS", age = runif(n, 50, 75), sex = rep(c(0, 1), n/2), icv_ml = rnorm(n, 1500, 80))
ls$wmh_region_ml <- rgamma(n, shape = 2, scale = exp(-.5 + .025 * (ls$age - 60)))
ls$wmh_region_ml[seq(1, n, by = 4)] <- 0
ls$wmh_allzero_ml <- 0
ls$wmh_positive_ml <- rgamma(n, shape = 2, scale = .5)
ls$muse_001_ml <- exp(3 - .006 * (ls$age - 60) + .3 * log(ls$icv_ml/1500) + rnorm(n, 0, .04))
ls_path <- file.path(scratch, "ls.tsv")
wg_write(ls, ls_path)

wmh_fit <- request("norm_fit", "wmh_fit", list(train = ls_path), modality = "wmh", feature = "wmh_region_ml")
run_normative(wmh_fit)
stopifnot(read_result(wmh_fit)$status == "completed", read_result(wmh_fit)$metrics$validation_n == n)
zero_fit <- request("norm_fit", "zero_fit", list(train = ls_path), modality = "wmh", feature = "wmh_allzero_ml")
run_normative(zero_fit)
stopifnot(read_result(zero_fit)$status == "not_estimable", !file.exists(file.path(zero_fit$output_dir, "model.rds")))
positive_fit <- request("norm_fit", "positive_fit", list(train = ls_path), modality = "wmh", feature = "wmh_positive_ml")
run_normative(positive_fit)
stopifnot(read_result(positive_fit)$status == "completed", read_result(positive_fit)$metrics$family == "GA")
gm_fit <- request("norm_fit", "gm_fit", list(train = ls_path), modality = "gmv", feature = "muse_001_ml")
run_normative(gm_fit)
stopifnot(read_result(gm_fit)$status == "completed", read_result(gm_fit)$metrics$validation_n == n)

query <- ls[1:8, ]
query$cohort <- "GB"
query$subject_key <- sub("LS", "GB", query$subject_key)
query$age[8] <- 85
query_path <- file.path(scratch, "query.tsv")
wg_write(query, query_path)
norm_apply <- request("norm_apply", "wmh_apply", list(data = query_path, model = file.path(wmh_fit$output_dir, "model.rds")))
run_normative(norm_apply)
norm_scores <- read_scores(norm_apply)
stopifnot(norm_scores$participant_id[1] == "00001", wg_num(norm_scores$wmh_region_ml_dev[1]) < 0,
          is.na(wg_num(norm_scores$wmh_region_ml_dev[8])), norm_scores$wmh_region_ml_valid[8] == "FALSE")

n <- 120L
latent <- rnorm(n)
train <- data.frame(subject_key = sprintf("C3:%05d", seq_len(n)), participant_id = sprintf("%05d", seq_len(n)), cohort = "C3")
for (j in 1:4) train[[paste0("f", j)]] <- j * latent + rnorm(n, 0, .10)
train_path <- file.path(scratch, "train.tsv")
wg_write(train, train_path)
features <- paste0("f", 1:4)
reference <- request("pca_fit", "reference", list(train = train_path), modality = "wmh", features = features, k = "max")
run_pca(reference)
reference_model <- file.path(reference$output_dir, "model.rds")
cv_files <- character()
bs_files <- character()
for (replicate in 1:2) {
  cv <- request("pca_cv", paste0("cv", replicate), list(train = train_path), modality = "wmh", features = features, replicate = replicate)
  run_pca(cv)
  cv_files <- c(cv_files, file.path(cv$output_dir, "criterion.tsv"))
  for (k in 1:2) {
    bs <- request("pca_bootstrap", paste0("bs", k, "_", replicate), list(train = train_path, reference = reference_model), modality = "wmh", k = k, replicate = replicate)
    run_pca(bs)
    bs_files <- c(bs_files, file.path(bs$output_dir, "bootstrap.tsv"))
  }
}
select <- request("pca_select", "select", list(cv = as.list(cv_files), bootstrap = as.list(bs_files)))
run_pca(select)
stopifnot(read_result(select)$status == "completed")
choice <- jsonlite::read_json(file.path(select$output_dir, "selection.json"), simplifyVector = TRUE)$selected_k
stopifnot(choice %in% 0:2)
incomplete_select <- request("pca_select", "incomplete_select", list(cv = as.list(cv_files[1]), bootstrap = as.list(bs_files)))
run_pca(incomplete_select)
stopifnot(read_result(incomplete_select)$status == "waiting_input", !file.exists(file.path(incomplete_select$output_dir, "selection.json")))
one_fit <- request("pca_fit", "one_fit", list(train = train_path), modality = "wmh", features = features, k = 1)
run_pca(one_fit)
one_model <- file.path(one_fit$output_dir, "model.rds")
hash <- unname(tools::md5sum(one_model))
new <- train[1:8, ]
new$subject_key <- sub("C3", "GB", new$subject_key)
new$cohort <- "GB"
new$f1[8] <- NA
new_path <- file.path(scratch, "new.tsv")
wg_write(new, new_path)
projection <- request("pca_project", "projection", list(data = new_path, model = one_model))
run_pca(projection)
projected <- read_scores(projection)
original <- read_scores(one_fit)
stopifnot(isTRUE(all.equal(wg_num(projected$WMH_PC01[1:7]), wg_num(original$WMH_PC01[1:7]), tolerance = 1e-10)),
          projected$WMH_projection_valid[8] == "FALSE", identical(hash, unname(tools::md5sum(one_model))))
zero_pca <- request("pca_fit", "zero_pca", list(train = train_path), modality = "wmh", features = features, k = 0)
run_pca(zero_pca)
stopifnot(read_result(zero_pca)$metrics$k == 0, !any(grepl("^WMH_PC", names(read_scores(zero_pca)))))

gm_pca <- request("pca_fit", "gm_pca", list(train = train_path), modality = "gmv", features = features, k = 1)
run_pca(gm_pca)
joint <- request("joint", "joint", list(wmh_scores = file.path(one_fit$output_dir, "scores.tsv"), gmv_scores = file.path(gm_pca$output_dir, "scores.tsv")))
run_pca(joint)
joined <- read_scores(joint)
stopifnot(abs(stats::var(wg_num(joined$JOINT_WMH_PC01)) - 1) < 1e-10,
          abs(stats::var(wg_num(joined$JOINT_GMV_PC01)) - 1) < 1e-10)

train$acute_volume_ml <- exp(rnorm(n, 2, .4))
train$age <- seq(50, 75, length.out = n)
train$r1 <- .03 * log1p(train$acute_volume_ml) + rnorm(n, 0, .005)
train$r2 <- .06 * log1p(train$acute_volume_ml) + rnorm(n, 0, .005)
acute_path <- file.path(scratch, "acute.tsv")
wg_write(train, acute_path)
acute <- request("acute_residualize", "acute", list(train = acute_path), features = c("r1", "r2"))
run_pca(acute)
residuals <- wg_read(file.path(acute$output_dir, "residuals.tsv"))
stopifnot(abs(stats::cor(wg_num(residuals$r1_resid), log1p(train$acute_volume_ml))) < 1e-8)
# 必须在全新R进程中验证冻结基函数；同一R会话已加载splines会掩盖S3注册缺陷。
fresh <- request("acute_residualize", "acute_fresh_process", list(model = file.path(acute$output_dir, "model.rds"), data = acute_path))
fresh_request <- file.path(scratch, "fresh_request.rds")
saveRDS(fresh, fresh_request)
fresh_script <- file.path(scratch, "fresh_process.R")
writeLines(c(
  sprintf("setwd(%s)", deparse(normalizePath("."))),
  ".libPaths(c('.Rlib', .libPaths()))",
  "source('R/common.R'); source('R/pca.R')",
  sprintf("run_pca(readRDS(%s))", deparse(fresh_request))
), fresh_script)
fresh_output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(fresh_script)), stdout = TRUE, stderr = TRUE)
fresh_exit <- attr(fresh_output, "status")
if (!is.null(fresh_exit) && fresh_exit != 0L) stop("Fresh-process acute projection failed: ", paste(fresh_output, collapse = "\n"))
fresh_residuals <- wg_read(file.path(fresh$output_dir, "residuals.tsv"))
stopifnot(isTRUE(all.equal(wg_num(fresh_residuals$r1_resid), wg_num(residuals$r1_resid), tolerance = 1e-10)))
age_query <- train
age_query$age[1:3] <- c(49, 76, NA)
age_query_path <- file.path(scratch, "acute_age_query.tsv")
wg_write(age_query, age_query_path)
age_projection <- request("acute_residualize", "acute_age_projection", list(model = file.path(acute$output_dir, "model.rds"), data = age_query_path))
run_pca(age_projection)
age_result <- wg_read(file.path(age_projection$output_dir, "residuals.tsv"))
stopifnot(all(is.na(wg_num(age_result$r1_resid[1:3]))), all(age_result$acute_residual_valid[1:3] == "FALSE"),
          read_result(age_projection)$metrics$age_ineligible_n == 3)
age_training <- request("acute_residualize", "acute_age_training", list(train = age_query_path), features = c("r1", "r2"))
run_pca(age_training)
stopifnot(readRDS(file.path(age_training$output_dir, "model.rds"))$training_n == n - 3L)
age_query$age <- NULL
no_age_path <- file.path(scratch, "acute_no_age.tsv")
wg_write(age_query, no_age_path)
no_age <- request("acute_residualize", "acute_no_age", list(model = file.path(acute$output_dir, "model.rds"), data = no_age_path))
run_pca(no_age)
stopifnot(read_result(no_age)$status == "waiting_input", !file.exists(file.path(no_age$output_dir, "residuals.tsv")))
cat("PASS: normative folds, zero handling, frozen projection, PCA CV/bootstrap/selection, K0, block scaling, acute residualization\n")
cat("Synthetic artifacts:", scratch, "\n")
