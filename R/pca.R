# 分模态PCA、矩阵遮蔽验证、患者bootstrap和冻结投射。
wg_pca_limit <- function(request, X) {
  limit <- if (request$modality == "acute") wg_cfg(request, "acute", "max_components") else wg_cfg(request, "pca", "max_components")
  max(0L, min(as.integer(limit), ncol(X) - 1L, nrow(X) - 2L,
              if (ncol(X) && nrow(X)) qr(scale(X, center = TRUE, scale = FALSE))$rank else 0L))
}

wg_pca_bundle <- function(X, modality, k, ratio) {
  wg_require("FactoMineR")
  if (any(!is.finite(X))) stop("PCA requires explicit complete-case handling; implicit mean imputation is forbidden")
  if (nrow(X) < 2L) stop("PCA requires at least two complete subjects")
  if (k < 0L || k > min(ncol(X), nrow(X) - 1L)) stop("Invalid PCA dimension")
  center <- colMeans(X)
  if (k == 0L) {
    return(list(model = NULL, features = colnames(X), center = center, k = 0L, modality = modality,
                loadings = matrix(numeric(), ncol(X), 0L, dimnames = list(colnames(X), character())),
                signs = numeric(), score_sd = numeric(), eigenvalues = numeric(), near_degenerate = logical(),
                versions = wg_versions("FactoMineR")))
  }
  fit <- FactoMineR::PCA(as.data.frame(X), scale.unit = FALSE, ncp = k, graph = FALSE)
  loadings <- as.matrix(fit$svd$V[, seq_len(k), drop = FALSE])
  rownames(loadings) <- colnames(X)
  signs <- vapply(seq_len(k), function(j) if (loadings[which.max(abs(loadings[, j])), j] < 0) -1 else 1, numeric(1))
  loadings <- sweep(loadings, 2, signs, "*")
  coords <- sweep(as.matrix(fit$ind$coord[, seq_len(k), drop = FALSE]), 2, signs, "*")
  eigenvalues <- as.numeric(fit$eig[, 1])
  near <- rep(FALSE, k)
  if (length(eigenvalues) > 1L) {
    for (j in seq_len(min(k, length(eigenvalues) - 1L))) {
      if (is.finite(eigenvalues[j + 1L]) && eigenvalues[j + 1L] > 0 && eigenvalues[j] / eigenvalues[j + 1L] < ratio) {
        near[j] <- TRUE
        if (j + 1L <= k) near[j + 1L] <- TRUE
      }
    }
  }
  list(model = fit, features = colnames(X), center = center, k = k, modality = modality,
       loadings = loadings, signs = signs, score_sd = apply(coords, 2, stats::sd),
       eigenvalues = eigenvalues, near_degenerate = near, versions = wg_versions("FactoMineR"))
}

wg_pca_coords <- function(bundle, X) {
  if (!identical(colnames(X), bundle$features)) stop("PCA feature order does not match frozen model")
  if (!bundle$k) return(matrix(numeric(), nrow(X), 0L))
  if (any(!is.finite(X))) stop("PCA projection cannot silently impute")
  projected <- stats::predict(bundle$model, newdata = as.data.frame(X))$coord[, seq_len(bundle$k), drop = FALSE]
  sweep(as.matrix(projected), 2, bundle$signs, "*")
}

wg_score_table <- function(data, X, valid, bundle) {
  out <- wg_ids(data)
  names_pc <- wg_pc_names(bundle$modality, bundle$k)
  sd_names <- if (length(names_pc)) paste0(names_pc, "_C3SD") else character()
  for (name in c(names_pc, sd_names)) out[[name]] <- NA_real_
  if (bundle$k && any(valid)) {
    score <- wg_pca_coords(bundle, X[valid, , drop = FALSE])
    colnames(score) <- names_pc
    out[valid, names_pc] <- score
    if (any(!is.finite(bundle$score_sd) | bundle$score_sd <= 0)) stop("Invalid frozen component scale")
    out[valid, paste0(names_pc, "_C3SD")] <- sweep(score, 2, bundle$score_sd, "/")
  }
  out[[paste0(wg_prefix(bundle$modality), "_projection_valid")]] <- valid
  out
}

wg_pca_fit <- function(request) {
  data <- wg_read(request$inputs$train)
  wg_assert_cohort(data, "C3")
  features <- wg_features(request)
  complete <- wg_complete_matrix(data, features)
  X <- complete$X
  maximum <- wg_pca_limit(request, X)
  k <- request$k
  if (!is.null(request$inputs$selection)) k <- jsonlite::read_json(request$inputs$selection, simplifyVector = TRUE)$selected_k
  if (is.null(k)) stop("PCA fit requires k or inputs.selection")
  if (identical(k, "max")) k <- maximum
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 0L || k > maximum) {
    return(wg_result(request, "not_estimable", "Requested k exceeds permissible rank", n_subjects = nrow(X)))
  }
  if (nrow(X) < 2L) return(wg_result(request, "not_estimable", "Fewer than two complete C3 subjects", n_subjects = nrow(X)))
  bundle <- wg_pca_bundle(X, request$modality, k, as.numeric(wg_cfg(request, "pca", "near_degenerate_ratio")))
  bundle$training_subject_keys <- complete$data$subject_key
  saveRDS(bundle, file.path(request$output_dir, "model.rds"))
  all_X <- wg_matrix(data, features)
  scores <- wg_score_table(data, all_X, complete$valid, bundle)
  wg_write(scores, file.path(request$output_dir, "scores.tsv"))
  loadings <- data.frame(feature = features, bundle$loadings, check.names = FALSE)
  if (k) names(loadings)[-1] <- wg_pc_names(request$modality, k)
  wg_write(loadings, file.path(request$output_dir, "loadings.tsv"))
  eig <- data.frame(component = seq_along(bundle$eigenvalues), eigenvalue = bundle$eigenvalues)
  eig$explained_fraction <- if (sum(bundle$eigenvalues) > 0) eig$eigenvalue / sum(bundle$eigenvalues) else numeric(nrow(eig))
  wg_write(eig, file.path(request$output_dir, "eigenvalues.tsv"))
  groups <- integer(k)
  if (k) {
    groups[1] <- 1L
    if (k > 1L) for (j in 2:k) {
      adjacent <- bundle$eigenvalues[j] > 0 && bundle$eigenvalues[j - 1L] / bundle$eigenvalues[j] < as.numeric(wg_cfg(request, "pca", "near_degenerate_ratio"))
      groups[j] <- groups[j - 1L] + as.integer(!adjacent)
    }
  }
  group_size <- if (k) as.integer(table(groups)[as.character(groups)]) else integer()
  component_metadata <- data.frame(component = wg_pc_names(request$modality, k), index = seq_len(k),
                                    group_id = groups, group_size = group_size,
                                    individually_identifiable = !bundle$near_degenerate,
                                    score_sd = bundle$score_sd)
  wg_write(component_metadata, file.path(request$output_dir, "component_metadata.tsv"))
  wg_png(request, "scree.png", function() {
    graphics::par(mar = c(4, 4, 3, 1))
    if (nrow(eig)) {
      graphics::plot(eig$component, eig$explained_fraction, type = "b", pch = 16, col = "#32658c",
                     xlab = "Component", ylab = "Fraction of matrix variance", main = paste(wg_prefix(request$modality), "PCA | retained k =", k))
      graphics::abline(v = k + .5, lty = 2, col = "#bd4538")
    } else {
      graphics::plot.new(); graphics::text(.5, .5, "K = 0: mean-only representation; no retained spatial component")
    }
  })
  wg_png(request, "loadings.png", function() {
    if (!k) { graphics::plot.new(); graphics::text(.5, .5, "K = 0: no component loading map"); return(invisible(NULL)) }
    maximum <- as.integer(wg_cfg(request, "figures", "max_heatmap_features"))
    selected <- head(order(rowSums(bundle$loadings^2), decreasing = TRUE), maximum)
    matrix <- bundle$loadings[selected, , drop = FALSE]
    bound <- max(abs(matrix))
    graphics::par(cex = .75, mar = c(5, 13, 6, 2))
    # 一列成分也画为矩形栅格，避免image单列坐标范围退化。
    graphics::plot(NA, xlim = c(.5, k + .5), ylim = c(.5, nrow(matrix) + .5), xaxt = "n", yaxt = "n",
                   xlab = "Frozen component", ylab = "", main = "Feature loadings (not a voxelwise brain map)", xaxs = "i", yaxs = "i")
    palette <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    for (i in seq_len(nrow(matrix))) for (j in seq_len(k)) {
      color <- palette[1L + floor(100 * (matrix[i, j] + bound) / (2 * bound))]
      graphics::rect(j - .5, i - .5, j + .5, i + .5, col = color, border = "white")
    }
    graphics::axis(1, seq_len(k), wg_pc_names(request$modality, k), las = 1, cex.axis = .8)
    graphics::axis(2, seq_len(nrow(matrix)), rownames(matrix), las = 2, cex.axis = .7)
    graphics::mtext(sprintf("Blue = %.3f | white = 0 | red = %.3f | top %d loading-energy features; full table: loadings.tsv", -bound, bound, nrow(matrix)), side = 3, line = .2, cex = .7)
  })
  wg_result(request, "completed", artifacts = c("model.rds", "scores.tsv", "loadings.tsv", "eigenvalues.tsv", "component_metadata.tsv", "scree.png", "loadings.png"),
            n_subjects = nrow(X), metrics = list(k = k, permissible_max = maximum, excluded_n = sum(!complete$valid), near_degenerate = as.list(bundle$near_degenerate)))
}

wg_pca_project <- function(request) {
  wg_require("FactoMineR")
  bundle <- readRDS(request$inputs$model)
  data <- wg_read(request$inputs$data)
  complete <- wg_complete_matrix(data, bundle$features)
  all_X <- wg_matrix(data, bundle$features)
  output <- wg_score_table(data, all_X, complete$valid, bundle)
  wg_write(output, file.path(request$output_dir, "scores.tsv"))
  wg_result(request, "completed", artifacts = "scores.tsv", n_subjects = nrow(data),
            metrics = list(k = bundle$k, valid_n = sum(complete$valid), invalid_n = sum(!complete$valid)))
}

wg_pca_cv <- function(request) {
  wg_require(c("missMDA", "FactoMineR"))
  data <- wg_read(request$inputs$train)
  wg_assert_cohort(data, "C3")
  complete <- wg_complete_matrix(data, wg_features(request))
  X <- complete$X
  if (nrow(X) < 3L || ncol(X) < 1L) return(wg_result(request, "not_estimable", "Insufficient matrix for masked validation", n_subjects = nrow(X)))
  cfg <- wg_cfg(request, "pca")
  maximum <- wg_pca_limit(request, X)
  replicate <- if (is.null(request$replicate)) 0L else as.integer(request$replicate)
  repeats <- if (replicate > 0L) 1L else as.integer(cfg$cv_repeats)
  wg_seed(request, replicate)
  masked_n <- round(as.numeric(cfg$mask_fraction) * length(X))
  if (masked_n < 1L) stop("Configured masking fraction produces no masked cells")
  # 按K分别调用同一成熟接口，重置同一随机种子以使用相同遮蔽。
  # 上游将EM迭代上限作为warning返回；必须记录为该候选不收敛，不能当有效误差。
  criterion <- do.call(rbind, lapply(0:maximum, function(k) {
    wg_seed(request, replicate)
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      missMDA::estim_ncpPCA(X, ncp.min = k, ncp.max = k, method = "EM", scale = FALSE,
                            method.cv = "Kfold", nbsim = repeats, pNA = as.numeric(cfg$mask_fraction),
                            threshold = as.numeric(cfg$convergence_tolerance), verbose = FALSE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    bad <- inherits(fit, "error") || any(grepl("Stopped after|not converg|non.converg", warnings, ignore.case = TRUE))
    value <- if (inherits(fit, "error")) NA_real_ else as.numeric(fit$criterion)[1] / masked_n
    bad <- bad || !is.finite(value)
    reason <- paste(unique(c(warnings, if (inherits(fit, "error")) conditionMessage(fit))), collapse = " | ")
    data.frame(k = k, msep = if (bad) NA_real_ else value, replicate = replicate,
               internal_repeats = repeats, masked_cells = masked_n, converged = !bad, diagnostic = reason)
  }))
  wg_write(criterion, file.path(request$output_dir, "criterion.tsv"))
  valid <- criterion[criterion$converged, , drop = FALSE]
  wg_result(request, "completed", artifacts = "criterion.tsv", n_subjects = nrow(X),
            metrics = list(permissible_max = maximum,
                           minimum_msep_k = if (nrow(valid)) valid$k[which.min(valid$msep)] else NA_integer_,
                           nonconverged_candidates = as.list(criterion$k[!criterion$converged]), validation_unit = "masked_matrix_cells"))
}

wg_icc_absolute <- function(x, y) {
  good <- is.finite(x) & is.finite(y)
  ratings <- cbind(x[good], y[good])
  n <- nrow(ratings)
  if (n < 3L) return(NA_real_)
  grand <- mean(ratings)
  row_mean <- rowMeans(ratings)
  col_mean <- colMeans(ratings)
  ms_row <- 2 * sum((row_mean - grand)^2) / (n - 1)
  ms_col <- n * sum((col_mean - grand)^2)
  residual <- ratings - row_mean - rep(col_mean, each = n) + grand
  ms_error <- sum(residual^2) / (n - 1)
  denom <- ms_row + ms_error + 2 * (ms_col - ms_error) / n
  if (denom <= .Machine$double.eps) return(NA_real_)
  (ms_row - ms_error) / denom
}

wg_pca_bootstrap <- function(request) {
  wg_require(c("FactoMineR", "clue"))
  data <- wg_read(request$inputs$train)
  wg_assert_cohort(data, "C3")
  reference <- readRDS(request$inputs$reference)
  complete <- wg_complete_matrix(data, reference$features)
  X <- complete$X
  if (!setequal(complete$data$subject_key, reference$training_subject_keys)) {
    stop("Bootstrap population differs from the frozen C3 reference population")
  }
  k <- as.integer(request$k)
  replicate <- as.integer(request$replicate)
  if (k > reference$k || k < 0L) return(wg_result(request, "not_estimable", "Candidate k is outside reference rank", n_subjects = nrow(X)))
  if (k == 0L) {
    output <- data.frame(k = 0L, replicate = replicate, component = 0L, max_angle = 0,
                         loading_similarity = 1, score_icc = 1, near_degenerate = FALSE,
                         cuts_degenerate_group = FALSE, oob_n = nrow(X))
  } else {
    cfg <- wg_cfg(request, "pca")
    wg_seed(request, 100000L + replicate + 1000L * k)
    index <- sample(seq_len(nrow(X)), replace = TRUE)
    if (qr(scale(X[index, , drop = FALSE], center = TRUE, scale = FALSE))$rank < k) stop("Bootstrap sample has insufficient rank")
    boot <- wg_pca_bundle(X[index, , drop = FALSE], request$modality, k, as.numeric(cfg$near_degenerate_ratio))
    ref_load <- reference$loadings[, seq_len(k), drop = FALSE]
    cross <- crossprod(ref_load, boot$loadings)
    match <- as.integer(clue::solve_LSAP(abs(cross), maximum = TRUE))
    signs <- ifelse(diag(cross[, match, drop = FALSE]) < 0, -1, 1)
    boot_load <- sweep(boot$loadings[, match, drop = FALSE], 2, signs, "*")
    singular <- svd(crossprod(ref_load, boot_load), nu = 0, nv = 0)$d
    max_angle <- max(acos(pmin(1, pmax(0, singular)))) * 180 / pi
    oob <- setdiff(seq_len(nrow(X)), unique(index))
    eval_index <- oob
    ref_small <- reference
    ref_small$k <- k
    ref_small$signs <- reference$signs[seq_len(k)]
    ref_score <- if (length(eval_index)) wg_pca_coords(ref_small, X[eval_index, , drop = FALSE]) else matrix(numeric(), 0, k)
    boot_score <- if (length(eval_index)) sweep(wg_pca_coords(boot, X[eval_index, , drop = FALSE])[, match, drop = FALSE], 2, signs, "*") else matrix(numeric(), 0, k)
    icc <- vapply(seq_len(k), function(j) wg_icc_absolute(ref_score[, j], boot_score[, j]), numeric(1))
    near <- reference$near_degenerate[seq_len(k)]
    cuts <- length(reference$eigenvalues) > k && reference$eigenvalues[k + 1L] > 0 &&
      reference$eigenvalues[k] / reference$eigenvalues[k + 1L] < as.numeric(cfg$near_degenerate_ratio)
    output <- data.frame(k = k, replicate = replicate, component = seq_len(k), max_angle = max_angle,
                         loading_similarity = abs(diag(crossprod(ref_load, boot_load))), score_icc = icc,
                         near_degenerate = near, cuts_degenerate_group = cuts, oob_n = length(oob))
  }
  wg_write(output, file.path(request$output_dir, "bootstrap.tsv"))
  wg_result(request, "completed", artifacts = "bootstrap.tsv", n_subjects = nrow(X),
            metrics = list(k = k, replicate = replicate, max_angle = output$max_angle[1]))
}

wg_read_metrics <- function(paths, required) {
  paths <- as.character(unlist(paths, use.names = FALSE))
  if (!length(paths)) return(data.frame())
  pieces <- lapply(paths, function(path) {
    if (!file.exists(path)) stop("Missing completed validation output: ", path)
    d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(required %in% names(d))) stop("Invalid validation columns: ", path)
    d
  })
  do.call(rbind, pieces)
}

wg_pca_select <- function(request) {
  cfg <- wg_cfg(request, "pca")
  cv <- wg_read_metrics(request$inputs$cv, c("k", "msep", "replicate", "internal_repeats"))
  bs <- wg_read_metrics(request$inputs$bootstrap, c("k", "replicate", "component", "max_angle", "loading_similarity", "score_icc", "near_degenerate", "cuts_degenerate_group"))
  if (!nrow(cv)) return(wg_result(request, "waiting_input", "No CV results supplied"))
  candidates <- sort(unique(cv$k))
  rows <- lapply(candidates, function(k) {
    values <- cv[cv$k == k, , drop = FALSE]
    samples <- bs[bs$k == k, , drop = FALSE]
    if (anyDuplicated(values$replicate)) stop("Repeated CV replicate for same candidate")
    cv_ok <- sum(values$internal_repeats) >= as.integer(cfg$cv_repeats)
    cv_converged <- all(is.finite(values$msep)) && (!"converged" %in% names(values) || all(values$converged))
    if (!cv_converged) return(data.frame(k = k, msep = NA_real_, cv_complete = cv_ok,
                                        bootstrap_n = length(unique(samples$replicate)), angle_median = NA_real_, angle_p95 = NA_real_,
                                        loading_min_median = NA_real_, score_icc_min_median = NA_real_, stable = FALSE, reason = "cv_nonconvergence"))
    if (k == 0L) return(data.frame(k = k, msep = weighted.mean(values$msep, values$internal_repeats),
                                  cv_complete = cv_ok, bootstrap_n = 0L, angle_median = 0, angle_p95 = 0,
                                  loading_min_median = 1, score_icc_min_median = 1, stable = cv_ok, reason = "mean_only_candidate"))
    bs_n <- length(unique(samples$replicate))
    if (!nrow(samples)) return(data.frame(k = k, msep = mean(values$msep), cv_complete = cv_ok, bootstrap_n = 0L,
                                         angle_median = NA_real_, angle_p95 = NA_real_, loading_min_median = NA_real_,
                                         score_icc_min_median = NA_real_, stable = FALSE, reason = "missing_bootstrap"))
    if (anyDuplicated(paste(samples$replicate, samples$component))) stop("Duplicate bootstrap component records")
    angles <- samples$max_angle[!duplicated(samples$replicate)]
    separate <- samples[!samples$near_degenerate, , drop = FALSE]
    med_min <- function(column) {
      if (!nrow(separate)) return(1)
      medians <- tapply(separate[[column]], separate$component, stats::median, na.rm = FALSE)
      if (any(!is.finite(medians))) return(NA_real_)
      min(medians)
    }
    loading <- med_min("loading_similarity")
    icc <- med_min("score_icc")
    amed <- stats::median(angles)
    ap95 <- as.numeric(stats::quantile(angles, .95, names = FALSE))
    complete_bs <- bs_n >= as.integer(cfg$bootstrap_repeats) && all(table(samples$replicate) == k)
    stable <- cv_ok && complete_bs && !any(samples$cuts_degenerate_group) &&
      is.finite(amed) && is.finite(ap95) && is.finite(loading) && is.finite(icc) &&
      amed <= cfg$angle_median_max && ap95 <= cfg$angle_p95_max &&
      loading >= cfg$loading_similarity_min && icc >= cfg$score_icc_min
    data.frame(k = k, msep = weighted.mean(values$msep, values$internal_repeats), cv_complete = cv_ok,
               bootstrap_n = bs_n, angle_median = amed, angle_p95 = ap95, loading_min_median = loading,
               score_icc_min_median = icc, stable = stable,
               reason = if (!cv_ok || !complete_bs) "incomplete_validation" else if (any(samples$cuts_degenerate_group)) "cuts_near_degenerate_group" else if (stable) "stable" else "stability_threshold_not_met")
  })
  selection <- do.call(rbind, rows)
  wg_write(selection, file.path(request$output_dir, "selection.tsv"))
  # 未完成验证不能悄悄强制降到0维。
  if (any(!selection$cv_complete) || any(selection$reason %in% c("missing_bootstrap", "incomplete_validation"))) {
    return(wg_result(request, "waiting_input", "Candidate validation is incomplete", "selection.tsv"))
  }
  ranked <- selection[selection$stable & is.finite(selection$msep), , drop = FALSE]
  ranked <- ranked[order(ranked$msep, ranked$k), , drop = FALSE]
  if (!nrow(ranked)) return(wg_result(request, "not_estimable", "No valid component candidate", "selection.tsv"))
  selected <- as.integer(ranked$k[1])
  wg_json(list(selected_k = selected, selection_rule = "minimum_masked_msep_among_stable_candidates",
               candidate_upper_bound = max(candidates), selected_at_upper_bound = selected == max(candidates)),
          file.path(request$output_dir, "selection.json"))
  wg_result(request, "completed", artifacts = c("selection.tsv", "selection.json"), metrics = list(selected_k = selected))
}

wg_joint <- function(request) {
  w <- wg_read(request$inputs$wmh_scores)
  g <- wg_read(request$inputs$gmv_scores)
  keys <- union(w$subject_key, g$subject_key)
  ids <- rbind(wg_ids(w), wg_ids(g)[!g$subject_key %in% w$subject_key, , drop = FALSE])
  ids <- ids[match(keys, ids$subject_key), , drop = FALSE]
  wcols <- grep("^WMH_PC[0-9]+$", names(w), value = TRUE)
  gcols <- grep("^GMV_PC[0-9]+$", names(g), value = TRUE)
  if (!is.null(request$inputs$model)) {
    bundle <- readRDS(request$inputs$model)
    if (!identical(wcols, bundle$wmh_columns) || !identical(gcols, bundle$gmv_columns)) stop("Joint model component mismatch")
  } else {
    wg_assert_cohort(w, "C3")
    wg_assert_cohort(g, "C3")
    complete_keys <- intersect(w$subject_key, g$subject_key)
    wx <- wg_matrix(w[match(complete_keys, w$subject_key), , drop = FALSE], wcols)
    gx <- wg_matrix(g[match(complete_keys, g$subject_key), , drop = FALSE], gcols)
    complete <- rowSums(!is.finite(cbind(wx, gx))) == 0
    if (sum(complete) < 2L) return(wg_result(request, "not_estimable", "Joint scales need at least two complete C3 subjects"))
    block_scale <- function(x) if (!ncol(x)) 0 else sqrt(sum(apply(x[complete, , drop = FALSE], 2, stats::var)))
    bundle <- list(wmh_columns = wcols, gmv_columns = gcols, wmh_scale = block_scale(wx), gmv_scale = block_scale(gx),
                   training_n = sum(complete), definition = "trace_covariance_block_balance")
    if ((length(wcols) && (!is.finite(bundle$wmh_scale) || bundle$wmh_scale <= 0)) ||
        (length(gcols) && (!is.finite(bundle$gmv_scale) || bundle$gmv_scale <= 0))) stop("Invalid joint block variance")
    saveRDS(bundle, file.path(request$output_dir, "model.rds"))
  }
  out <- ids
  for (name in wcols) out[[paste0("JOINT_", name)]] <- wg_num(w[[name]][match(keys, w$subject_key)]) / bundle$wmh_scale
  for (name in gcols) out[[paste0("JOINT_", name)]] <- wg_num(g[[name]][match(keys, g$subject_key)]) / bundle$gmv_scale
  jcols <- grep("^JOINT_", names(out), value = TRUE)
  out$joint_valid <- if (length(jcols)) rowSums(!is.finite(wg_matrix(out, jcols))) == 0 else FALSE
  out$joint_dimension <- length(jcols)
  wg_write(out, file.path(request$output_dir, "scores.tsv"))
  artifacts <- c("scores.tsv", if (is.null(request$inputs$model)) "model.rds")
  wg_result(request, "completed", if (!length(jcols)) "Both modalities selected k=0; no joint spatial coordinates" else "",
            artifacts, nrow(out), list(dimension = length(jcols), wmh_scale = bundle$wmh_scale, gmv_scale = bundle$gmv_scale))
}

wg_acute_residualize <- function(request) {
  # RDS保存的ns基函数不会自动加载其S3方法所属命名空间。
  wg_require("splines")
  features <- wg_features(request)
  cfg <- wg_cfg(request, "acute")
  age_range <- as.numeric(unlist(wg_cfg(request, "age_range")))
  if (is.null(request$inputs$model) && !length(features)) {
    return(wg_result(request, "not_estimable", "No acute region passed configured coverage/prevalence checks"))
  }
  if (!is.null(request$inputs$model)) {
    bundle <- readRDS(request$inputs$model)
    if (is.null(bundle$age_range)) {
      return(wg_result(request, "not_estimable", "Frozen acute model lacks age-range provenance; refit required"))
    }
    age_range <- bundle$age_range
    data <- wg_read(request$inputs$data)
    features <- bundle$features
  } else {
    data <- wg_read(request$inputs$train)
    wg_assert_cohort(data, "C3")
    if (!"age" %in% names(data)) {
      return(wg_result(request, "waiting_input", "Acute C3 training table lacks mandatory age for eligibility", n_subjects = nrow(data)))
    }
    if (!"acute_volume_ml" %in% names(data)) stop("Missing acute_volume_ml")
    vol <- wg_num(data$acute_volume_ml)
    age <- wg_num(data$age)
    X <- wg_matrix(data, features)
    good <- is.finite(age) & age >= age_range[1] & age <= age_range[2] &
      is.finite(vol) & vol > 0 & rowSums(!is.finite(X)) == 0
    if (sum(good) <= as.integer(cfg$volume_spline_df) + 2L) return(wg_result(request, "not_estimable", "Insufficient acute-volume reference data", n_subjects = sum(good)))
    basis <- splines::ns(log1p(vol[good]), df = as.integer(cfg$volume_spline_df))
    fit_data <- as.data.frame(basis)
    colnames(fit_data) <- paste0("volume_ns", seq_len(ncol(basis)))
    fits <- lapply(seq_along(features), function(j) {
      d <- fit_data
      d$response <- X[good, j]
      stats::lm(stats::reformulate(names(fit_data), response = "response"), data = d)
    })
    bundle <- list(features = features, models = fits, basis = basis, basis_names = names(fit_data),
                   training_n = sum(good), age_range = age_range)
    saveRDS(bundle, file.path(request$output_dir, "model.rds"))
    if (!is.null(request$inputs$data)) data <- wg_read(request$inputs$data)
  }
  if (!"age" %in% names(data)) {
    return(wg_result(request, "waiting_input", "Acute projection table lacks mandatory age for eligibility", n_subjects = nrow(data)))
  }
  vol <- wg_num(data$acute_volume_ml)
  age <- wg_num(data$age)
  X <- wg_matrix(data, features)
  eligible_age <- is.finite(age) & age >= age_range[1] & age <= age_range[2]
  good <- eligible_age & is.finite(vol) & vol > 0 & rowSums(!is.finite(X)) == 0
  out <- wg_ids(data)
  out$age <- age
  out$acute_volume_ml <- vol
  residual_names <- paste0(features, "_resid")
  for (name in residual_names) out[[name]] <- NA_real_
  if (any(good)) {
    basis <- as.data.frame(stats::predict(bundle$basis, log1p(vol[good])))
    names(basis) <- bundle$basis_names
    for (j in seq_along(features)) out[[residual_names[j]]][good] <- X[good, j] - stats::predict(bundle$models[[j]], newdata = basis)
  }
  out$acute_residual_valid <- good
  out$acute_age_eligible <- eligible_age
  wg_write(out, file.path(request$output_dir, "residuals.tsv"))
  wg_json(list(features = as.list(residual_names)), file.path(request$output_dir, "features.json"))
  wg_result(request, "completed", artifacts = c("residuals.tsv", "features.json", if (is.null(request$inputs$model)) "model.rds"),
            n_subjects = nrow(data), metrics = list(valid_n = sum(good), age_ineligible_n = sum(!eligible_age)))
}

run_pca <- function(request) {
  dir.create(request$output_dir, recursive = TRUE, showWarnings = FALSE)
  switch(request$action, pca_reference = { request$k <- "max"; wg_pca_fit(request) },
         pca_fit = wg_pca_fit(request), pca_project = wg_pca_project(request), pca_cv = wg_pca_cv(request),
         pca_bootstrap = wg_pca_bootstrap(request), pca_select = wg_pca_select(request),
         joint = wg_joint(request), acute_residualize = wg_acute_residualize(request),
         stop("Unknown PCA action: ", request$action))
}
