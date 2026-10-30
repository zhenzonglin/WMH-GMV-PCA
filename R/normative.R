# GAMLSS健康参考。原始训练数据和年龄基函数一同保存，GB只预测不重拟合。
wg_norm_data <- function(data, feature, modality, age_range) {
  need <- c("age", "sex", "icv_ml", feature)
  if (!all(need %in% names(data))) stop("Missing normative columns: ", paste(setdiff(need, names(data)), collapse = ", "))
  d <- data.frame(age = wg_num(data$age), sex = wg_num(data$sex),
                  icv_ml = wg_num(data$icv_ml), observed = wg_num(data[[feature]]))
  valid <- is.finite(d$age) & d$age >= age_range[1] & d$age <= age_range[2] &
    is.finite(d$sex) & d$sex %in% c(0, 1) & is.finite(d$icv_ml) & d$icv_ml > 0 &
    is.finite(d$observed) & if (modality == "gmv") d$observed > 0 else d$observed >= 0
  d$log_icv <- log(d$icv_ml)
  d$response <- if (modality == "gmv") log(d$observed) else d$observed
  list(data = d, valid = valid)
}

wg_norm_fit_one <- function(data, modality, cfg) {
  if (length(unique(data$sex)) < 2L) stop("Reference training subset does not contain both sex categories")
  basis <- splines::ns(data$age, df = as.integer(cfg$age_spline_df))
  bnames <- paste0("age_ns", seq_len(ncol(basis)))
  data[bnames] <- as.data.frame(basis)
  mu_formula <- stats::reformulate(c(bnames, "sex", "log_icv"), response = "response")
  family_name <- if (modality == "gmv") "NO" else if (any(data$response == 0)) "ZAGA" else "GA"
  fam <- get(family_name, envir = asNamespace("gamlss.dist"))()
  args <- list(formula = mu_formula, sigma.formula = ~1, family = fam, data = data,
               control = gamlss::gamlss.control(n.cyc = as.integer(cfg$max_cycles), trace = FALSE))
  if (family_name == "ZAGA") args$nu.formula <- ~age + sex
  model <- do.call(gamlss::gamlss, args)
  if (!isTRUE(model$converged)) stop("GAMLSS did not converge")
  if (any(!is.finite(model$sigma.fv) | model$sigma.fv <= sqrt(.Machine$double.eps))) {
    stop("Degenerate reference residual scale; no artificial small-variance floor is applied")
  }
  list(model = model, data = data, age_basis = basis, basis_names = bnames, family = family_name)
}

wg_norm_predict <- function(bundle, data, modality, clip, interval_level) {
  # 全新R进程读取RDS时，需先注册splines的S3预测方法。
  wg_require("splines")
  if (!nrow(data)) return(data.frame())
  age_basis <- stats::predict(bundle$age_basis, data$age)
  data[bundle$basis_names] <- as.data.frame(age_basis)
  # 显式data防止RDS恢复后找不到原函数局部变量。
  parameters <- gamlss::predictAll(bundle$model, newdata = data, data = bundle$data,
                                    type = "response", output = "list")
  mu <- as.numeric(parameters$mu)
  sigma <- as.numeric(parameters$sigma)
  nu <- if (bundle$family == "ZAGA") as.numeric(parameters$nu) else rep(0, nrow(data))
  lower_probability <- (1 - interval_level) / 2
  upper_probability <- (1 + interval_level) / 2
  if (any(!is.finite(mu) | !is.finite(sigma) | sigma <= 0 | !is.finite(nu))) stop("Nonfinite reference prediction")
  if (modality == "gmv") {
    u <- stats::pnorm(data$response, mu, sigma)
    dev <- (mu - data$response) / sigma
    lower <- exp(stats::qnorm(lower_probability, mu, sigma))
    upper <- exp(stats::qnorm(upper_probability, mu, sigma))
    positive_u <- u
  } else {
    positive_u <- gamlss.dist::pGA(data$response, mu = mu, sigma = sigma)
    u <- nu + (1 - nu) * positive_u
    u[data$observed == 0] <- nu[data$observed == 0] / 2
    dev <- stats::qnorm(pmin(1 - clip, pmax(clip, u)))
    if (bundle$family == "ZAGA") {
      lower <- gamlss.dist::qZAGA(lower_probability, mu, sigma, nu)
      upper <- gamlss.dist::qZAGA(upper_probability, mu, sigma, nu)
    } else {
      lower <- gamlss.dist::qGA(lower_probability, mu, sigma)
      upper <- gamlss.dist::qGA(upper_probability, mu, sigma)
    }
  }
  data.frame(deviation = dev, percentile = u, cdf_clipped = u < clip | u > 1 - clip,
             predicted_zero = nu, predicted_mu = mu, predicted_sigma = sigma,
             lower_interval = lower, upper_interval = upper, positive_percentile = positive_u)
}

wg_norm_fit <- function(request) {
  wg_require(c("gamlss", "gamlss.dist"))
  cfg <- wg_cfg(request, "normative")
  age_range <- as.numeric(unlist(wg_cfg(request, "age_range")))
  data <- wg_read(request$inputs$train)
  wg_assert_cohort(data, "LS")
  feature <- request$feature
  modality <- request$modality
  nd <- wg_norm_data(data, feature, modality, age_range)
  d <- nd$data[nd$valid, , drop = FALSE]
  ids <- wg_ids(data[nd$valid, , drop = FALSE])
  if (nrow(d) < as.integer(cfg$folds) * 2L ||
      (modality == "wmh" && sum(d$observed > 0) < as.integer(cfg$min_positive))) {
    wg_json(list(feature = feature, n = nrow(d), positive_n = sum(d$observed > 0)), file.path(request$output_dir, "metrics.json"))
    return(wg_result(request, "not_estimable", "Insufficient healthy reference observations/positive WMH", "metrics.json", nrow(d)))
  }
  bundle <- wg_norm_fit_one(d, modality, cfg)
  bundle$feature <- feature
  bundle$modality <- modality
  bundle$age_range <- age_range
  bundle$cdf_clip <- as.numeric(cfg$cdf_clip)
  bundle$interval_level <- as.numeric(cfg$prediction_interval)
  bundle$versions <- wg_versions(c("gamlss", "gamlss.dist"))
  bundle$config <- cfg
  saveRDS(bundle, file.path(request$output_dir, "model.rds"))
  folds <- as.integer(cfg$folds)
  wg_seed(request)
  fold_id <- integer(nrow(d))
  strata <- interaction(d$sex, if (modality == "wmh") d$observed > 0 else rep(TRUE, nrow(d)), drop = TRUE)
  for (index in split(seq_len(nrow(d)), strata)) fold_id[index] <- sample(rep(seq_len(folds), length.out = length(index)))
  calibration <- cbind(ids, observed = d$observed, age = d$age, sex = d$sex, fold = fold_id)
  cols <- c("deviation", "percentile", "cdf_clipped", "predicted_zero", "predicted_mu", "predicted_sigma", "lower_interval", "upper_interval", "positive_percentile")
  for (name in cols) calibration[[name]] <- NA_real_
  calibration$fold_status <- "not_run"
  calibration$fold_reason <- ""
  for (fold in seq_len(folds)) {
    train <- which(fold_id != fold)
    test <- which(fold_id == fold)
    if (!length(test)) next
    value <- tryCatch({
      if (modality == "wmh" && sum(d$observed[train] > 0) < as.integer(cfg$min_positive)) stop("Training fold has insufficient positive WMH observations")
      fit <- wg_norm_fit_one(d[train, , drop = FALSE], modality, cfg)
      wg_norm_predict(fit, d[test, , drop = FALSE], modality, as.numeric(cfg$cdf_clip), as.numeric(cfg$prediction_interval))
    }, error = function(e) e)
    if (inherits(value, "error")) {
      calibration$fold_status[test] <- "failed"
      calibration$fold_reason[test] <- conditionMessage(value)
    } else {
      calibration[test, cols] <- value[, cols, drop = FALSE]
      calibration$fold_status[test] <- "completed"
    }
  }
  wg_write(calibration, file.path(request$output_dir, "calibration.tsv"))
  good <- calibration$fold_status == "completed"
  covered <- calibration$observed >= calibration$lower_interval & calibration$observed <= calibration$upper_interval
  metrics <- list(feature = feature, family = bundle$family, n = nrow(d),
                  positive_n = sum(d$observed > 0), n_excluded = sum(!nd$valid),
                  validation_n = sum(good), validation_failed_n = sum(!good),
                  interval_level = as.numeric(cfg$prediction_interval),
                  coverage_interval = if (any(good)) mean(covered[good]) else NA_real_,
                  deviation_mean = if (any(good)) mean(calibration$deviation[good]) else NA_real_,
                  deviation_sd = if (sum(good) > 1L) stats::sd(calibration$deviation[good]) else NA_real_,
                  calibration_status = if (all(good)) "evaluated_not_automatically_certified" else "incomplete",
                  zero_brier = if (modality == "wmh" && any(good)) mean(((calibration$observed[good] == 0) - calibration$predicted_zero[good])^2) else NA_real_)
  wg_json(metrics, file.path(request$output_dir, "metrics.json"))
  wg_png(request, "calibration.png", function() {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0), cex = .85)
    graphics::plot(d$age, d$observed, pch = 16, cex = .55, col = grDevices::adjustcolor("#32658c", .5),
                   xlab = "Age", ylab = "Observed volume (mL)", main = "Healthy reference observations")
    if (any(good)) {
      pred <- if (modality == "wmh") calibration$predicted_mu * (1 - calibration$predicted_zero) else exp(calibration$predicted_mu + calibration$predicted_sigma^2 / 2)
      graphics::plot(pred[good], calibration$observed[good], pch = 16, cex = .55,
                     col = grDevices::adjustcolor("#32658c", .5), xlab = "Out-of-fold expected volume (mL)",
                     ylab = "Observed volume (mL)", main = "Held-out conditional mean")
      graphics::abline(0, 1, col = "#bd4538", lty = 2)
      pit <- calibration$positive_percentile[good & (modality == "gmv" | calibration$observed > 0)]
      graphics::hist(pit, breaks = seq(0, 1, length.out = 11), col = "#80abc3", border = "white",
                     main = if (modality == "wmh") "Positive-part held-out CDF" else "Held-out conditional CDF",
                     xlab = "CDF percentile", xlim = c(0, 1))
      fold_coverage <- tapply(covered[good], calibration$fold[good], mean)
      graphics::barplot(fold_coverage, ylim = c(0, 1), col = "#80abc3", xlab = "Validation fold",
                        ylab = "Prediction-interval coverage", main = "Held-out coverage")
      graphics::abline(h = cfg$prediction_interval, col = "#bd4538", lty = 2)
    } else {
      for (j in seq_len(3L)) { graphics::plot.new(); graphics::text(.5, .5, "Validation unavailable: inspect calibration.tsv") }
    }
    graphics::mtext(paste(feature, "|", modality, "|", metrics$calibration_status), outer = TRUE, font = 2)
  })
  wg_result(request, "completed", artifacts = c("model.rds", "calibration.tsv", "metrics.json", "calibration.png"), n_subjects = nrow(d), metrics = metrics)
}

wg_norm_apply <- function(request) {
  wg_require(c("gamlss", "gamlss.dist"))
  bundle <- readRDS(request$inputs$model)
  data <- wg_read(request$inputs$data)
  feature <- bundle$feature
  nd <- wg_norm_data(data, feature, bundle$modality, bundle$age_range)
  result <- wg_ids(data)
  result[[paste0(feature, "_dev")]] <- NA_real_
  result[[paste0(feature, "_percentile")]] <- NA_real_
  result[[paste0(feature, "_valid")]] <- nd$valid
  result[[paste0(feature, "_cdf_clipped")]] <- NA
  if (any(nd$valid)) {
    prediction <- wg_norm_predict(bundle, nd$data[nd$valid, , drop = FALSE], bundle$modality, bundle$cdf_clip, bundle$interval_level)
    result[[paste0(feature, "_dev")]][nd$valid] <- prediction$deviation
    result[[paste0(feature, "_percentile")]][nd$valid] <- prediction$percentile
    result[[paste0(feature, "_cdf_clipped")]][nd$valid] <- prediction$cdf_clipped
  }
  wg_write(result, file.path(request$output_dir, "scores.tsv"))
  wg_result(request, "completed", artifacts = "scores.tsv", n_subjects = nrow(data),
            metrics = list(valid_n = sum(nd$valid), invalid_n = sum(!nd$valid)))
}

run_normative <- function(request) {
  dir.create(request$output_dir, recursive = TRUE, showWarnings = FALSE)
  switch(request$action, norm_fit = wg_norm_fit(request), norm_apply = wg_norm_apply(request),
         stop("Unknown normative action: ", request$action))
}
