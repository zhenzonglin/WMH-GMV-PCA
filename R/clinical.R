# 有序临床模型；调用 ordinal 求解器，禁止用结局重新选择影像成分。
clin_result <- function(req, status, reason, artifacts=character(), n_subjects=0L, extra=list()) {
  ans <- c(list(status=status, reason=reason, artifacts=as.list(unname(artifacts)), n_subjects=n_subjects), extra)
  jsonlite::write_json(ans, file.path(req$output_dir, "result.json"), auto_unbox=TRUE, pretty=TRUE, null="null", na="null")
  invisible(ans)
}
clin_table <- function(x, file) {
  utils::write.table(x, file, sep="\t", row.names=FALSE, quote=FALSE, na="")
}
clin_data <- function(path) {
  d <- utils::read.delim(path, colClasses="character", check.names=FALSE, na.strings=c("", "NA", "NaN"))
  keep <- c("subject_key", "participant_id", "cohort", "mrs_origin", "visit_date")
  for (n in setdiff(names(d), keep)) {
    if (!grepl("^(source_|mrs_origin_|visit_date_|index_date|death_date)", n))
      suppressWarnings(d[[n]] <- as.numeric(d[[n]]))
  }
  d
}
clin_ready <- function(fit, gradient_tolerance) {
  beta <- fit$beta
  if (any(!is.finite(beta))) stop("Non-finite or aliased regression coefficient; model not accepted")
  v <- stats::vcov(fit)
  if (any(!is.finite(v)) || any(diag(v) <= 0)) stop("Invalid covariance/Hessian; model not accepted")
  gradient <- if (!is.null(fit$gradient)) max(abs(fit$gradient)) else NA_real_
  convergence <- fit$optRes$convergence
  method<-if(!is.null(fit$control$method)) fit$control$method else "nlminb"
  success_codes<-if(identical(method,"ucminf")) c(1L,2L) else 0L
  if (!is.null(convergence) && length(convergence)==1L && !convergence %in% success_codes)
    stop("ordinal optimizer did not converge: ", convergence)
  if (is.finite(gradient) && gradient > gradient_tolerance) stop("ordinal maximum gradient exceeds configured tolerance")
  eigenvalues<-eigen(fit$Hessian,symmetric=TRUE,only.values=TRUE)$values
  list(max_gradient=gradient, log_likelihood=as.numeric(logLik(fit)), df=attr(logLik(fit), "df"),
       hessian_min_eigenvalue=min(eigenvalues),hessian_condition=max(eigenvalues)/min(eigenvalues),
       random_intercept_sd=if(!is.null(fit$ST)) as.numeric(unlist(fit$ST))[1] else NA_real_)
}
clin_effects <- function(fit, name, alpha) {
  b <- fit$beta
  se <- sqrt(diag(vcov(fit)))[names(b)]
  z <- b/se
  critical <- qnorm(1-alpha/2)
  data.frame(model=name, term=names(b), estimate=unname(b), se=unname(se),
    OR=exp(unname(b)), lower=exp(unname(b-critical*se)), upper=exp(unname(b+critical*se)),
    p=2*pnorm(-abs(unname(z))), q_secondary=NA_real_,
    interpretation="OR greater than 1 indicates a shift toward worse mRS", stringsAsFactors=FALSE)
}
clin_lrt <- function(smaller, larger, label) {
  ll0 <- logLik(smaller); ll1 <- logLik(larger)
  df <- attr(ll1,"df")-attr(ll0,"df")
  statistic <- 2*(as.numeric(ll1)-as.numeric(ll0))
  data.frame(test=label, statistic=if (df>0) statistic else NA_real_, df=df,
             p=if(df>0 && statistic>=0) pchisq(statistic,df,lower.tail=FALSE) else NA_real_,
             status=if(df>0) "tested" else "not_estimable_no_added_dimension")
}
clin_probabilities <- function(fit, d) {
  prediction_data <- d[, setdiff(names(d),c("mrs", "mrs_m12")), drop=FALSE]
  p <- predict(fit, newdata=prediction_data, type="prob")$fit
  if (is.null(dim(p))) p <- matrix(p,nrow=1)
  if (ncol(p)!=7L) stop("Prediction does not contain all seven mRS categories")
  if (any(!is.finite(p)) || max(abs(rowSums(p)-1))>1e-6) stop("Invalid ordinal probabilities")
  p
}
clin_prediction_frame <- function(p,d,model,label) {
  ans <- data.frame(subject_key=d$subject_key, model=model, validation=label,
                    observed_mrs=as.integer(d$mrs)-1L, p, check.names=FALSE)
  names(ans)[5:11] <- paste0("p_mrs",0:6)
  ans$p_independent <- rowSums(p[,1:3,drop=FALSE]); ans
}
clin_scores <- function(p,y,epsilon) {
  cum <- t(apply(p,1,cumsum))[,1:6,drop=FALSE]
  target <- outer(as.numeric(y),0:5,"<=")
  data.frame(metric=c("ranked_probability_score", "ordinal_log_loss"),
    value=c(mean(rowMeans((cum-target)^2)), -mean(log(pmax(epsilon,p[cbind(seq_along(y),as.numeric(y)+1)])))))
}
clin_calibration <- function(predictions,bins) {
  out <- list(); at <- 0L
  for (label in unique(predictions$validation)) for (model in unique(predictions$model)) {
    d <- predictions[predictions$validation==label & predictions$model==model,,drop=FALSE]
    if (!nrow(d)) next
    for (cut in 0:5) {
      probability <- rowSums(d[,paste0("p_mrs",0:cut),drop=FALSE])
      group <- pmin(bins, floor(probability*bins)+1L)
      for (bin in sort(unique(group))) {
        indices <- which(group==bin); at <- at+1L
        out[[at]] <- data.frame(validation=label,model=model,threshold=cut,bin=bin,n=length(indices),
          predicted=mean(probability[indices]),observed=mean(d$observed_mrs[indices]<=cut))
      }
    }
  }
  if(length(out)) do.call(rbind,out) else data.frame()
}
clin_timepoint_effects <- function(fit,components,visits,alpha) {
  # 线性对比只用于报告各访视的条件OR，不增加一组显著性检验。
  empty<-data.frame(visit=character(),component=character(),log_OR=numeric(),se=numeric(),
    OR=numeric(),lower=numeric(),upper=numeric(),estimand=character(),
    reference_visit=character(),contrast_terms=character(),stringsAsFactors=FALSE)
  if(!length(components)) return(empty)
  visits<-as.character(unlist(visits))
  beta<-fit$beta
  covariance<-vcov(fit)[names(beta),names(beta),drop=FALSE]
  critical<-qnorm(1-alpha/2)
  rows<-list()
  for(component in components) for(visit in visits) {
    if(!component %in% names(beta)) stop("Missing longitudinal main coefficient: ",component)
    terms<-component
    if(visit!=visits[1]) {
      possibilities<-c(paste0("visit",visit,":",component),paste0(component,":visit",visit))
      interaction<-intersect(possibilities,names(beta))
      if(length(interaction)!=1L) stop("Missing or ambiguous longitudinal interaction for ",component," at visit ",visit)
      terms<-c(terms,interaction)
    }
    contrast<-setNames(rep(0,length(beta)),names(beta));contrast[terms]<-1
    estimate<-as.numeric(sum(contrast*beta))
    variance<-as.numeric(t(contrast)%*%covariance%*%contrast)
    if(!is.finite(variance) || variance<=0) stop("Invalid timepoint contrast variance")
    se<-sqrt(variance)
    rows[[length(rows)+1L]]<-data.frame(visit=visit,component=component,log_OR=estimate,se=se,
      OR=exp(estimate),lower=exp(estimate-critical*se),upper=exp(estimate+critical*se),
      estimand="subject_specific_conditional_OR",reference_visit=visits[1],
      contrast_terms=paste(terms,collapse=" + "),stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}
clin_interaction_components <- function(columns, path, limit) {
  if (!length(columns)) return(character())
  if(is.null(path) || !file.exists(path)) stop("Component metadata required to avoid truncating a near-degenerate group")
  metadata <- utils::read.delim(path,stringsAsFactors=FALSE,check.names=FALSE)
  needed <- c("component","index","group_id","group_size","individually_identifiable")
  if(!all(needed %in% names(metadata)) || !all(columns %in% metadata$component)) stop("Incomplete component metadata")
  metadata <- metadata[metadata$component %in% columns,,drop=FALSE]
  candidates <- metadata[metadata$index<=limit,,drop=FALSE]
  keep <- vapply(candidates$group_id,function(group) {
    full <- metadata[metadata$group_id==group,,drop=FALSE]
    all(full$index<=limit) && nrow(full)==full$group_size[1]
  },logical(1))
  candidates$component[keep]
}
clin_clmm_control <- function(cfg) {
  opt<-cfg$optimizer
  if(is.null(opt)) stop("CLMM optimizer settings must be provided in configuration")
  if(identical(opt$method,"ucminf")) return(ordinal::clmm.control(method="ucminf",gradTol=opt$inner_gradient_tolerance,
    grtol=opt$outer_gradient_tolerance,xtol=opt$parameter_tolerance,maxeval=opt$max_evaluations,grad="central"))
  if(!identical(opt$method,"nlminb")) stop("Unsupported configured ordinal optimizer")
  ordinal::clmm.control(method="nlminb",gradTol=opt$inner_gradient_tolerance,
    rel.tol=opt$relative_tolerance,x.tol=opt$parameter_tolerance,
    iter.max=opt$max_iterations,eval.max=opt$max_evaluations)
}

run_clinical <- function(req) {
  if (!requireNamespace("ordinal",quietly=TRUE)) stop("Required mature backend is missing: ordinal")
  if (!requireNamespace("jsonlite",quietly=TRUE)) stop("Required package is missing: jsonlite")
  dir.create(req$output_dir,recursive=TRUE,showWarnings=FALSE)
  cfg <- req$config$analysis$clinical
  input <- req$inputs$data
  if (is.null(input) || !file.exists(input)) return(clin_result(req,"waiting_input","Clinical table is unavailable"))
  d <- clin_data(input)
  action <- req$action
  predict_only <- identical(action,"clinical_predict")
  part <- if (!is.null(req$analysis_part)) req$analysis_part else if (grepl("longitudinal$",action)) "longitudinal" else "single"
  b2 <- identical(action,"clinical_b2")
  cohort <- if (!is.null(req$cohort)) req$cohort else unique(d$cohort)[1]
  if (part=="longitudinal") {
    if (is.null(req$inputs$long) || !file.exists(req$inputs$long)) return(clin_result(req,"waiting_input","Longitudinal clinical table is unavailable"))
    d <- clin_data(req$inputs$long)
  }
  reference <- NULL
  if (!is.null(req$inputs$reference_models) && file.exists(req$inputs$reference_models)) reference <- readRDS(req$inputs$reference_models)
  if ((cohort=="GB" || predict_only) && is.null(reference)) return(clin_result(req,"waiting_input","Frozen C3 clinical specification and score scaling are required"))
  w <- grep("^WMH_PC[0-9]+$",names(d),value=TRUE)
  g <- grep("^GMV_PC[0-9]+$",names(d),value=TRUE)
  a <- grep("^ACUTE_PC[0-9]+$",names(d),value=TRUE)
  w <- sort(w); g <- sort(g); a <- sort(a)
  if (!is.null(reference)) {
    w <- reference$wmh_columns; g <- reference$gmv_columns; a <- reference$acute_columns
  }
  score_names <- c(w,g,if(b2) a else character())
  baseline <- unlist(cfg$covariates)
  structural <- c("icv_ml","wmh_total_ml","gmv_total_ml")
  outcome <- if(part=="longitudinal") "mrs" else paste0("mrs_m",max(unlist(cfg$visits)))
  required <- c("subject_key", baseline, structural,score_names,outcome,if(b2) "acute_volume_ml",if(part=="longitudinal") "visit")
  absent <- setdiff(required,names(d))
  if(length(absent)) return(clin_result(req,"waiting_input",paste("Missing fields:",paste(absent,collapse=", "))))
  if (b2 && (!length(a) || !length(c(w,g)))) return(clin_result(req,"not_estimable","Stable acute and chronic dimensions are required for location interactions"))
  # 所有嵌套模型使用同一完整病例集；不按各模型分别删缺失。
  eligible <- complete.cases(d[,required,drop=FALSE])
  eligible <- eligible & d$age>=req$config$analysis$age_range[[1]] & d$age<=req$config$analysis$age_range[[2]]
  if(isTRUE(cfg$exclude_conflicts) && "clinical_conflict" %in% names(d)) eligible <- eligible & !is.na(d$clinical_conflict) & d$clinical_conflict==0
  eligible <- eligible & d$icv_ml>0 & d$gmv_total_ml>0 & d$wmh_total_ml>=0
  if(b2) eligible <- eligible & d$acute_volume_ml>=0
  d <- d[which(eligible),,drop=FALSE]
  n <- length(unique(d$subject_key))
  if(n<if(predict_only) 1L else cfg$min_complete_cases) return(clin_result(req,"not_estimable","Too few complete patients",n_subjects=n))
  y <- d[[outcome]]
  if(any(!y %in% 0:6)) stop("Invalid mRS reached model backend")
  if(!predict_only && !all(0:6 %in% y)) return(clin_result(req,"not_estimable","One or more mRS categories absent; categories are not silently merged",n_subjects=n))
  d$mrs <- ordered(y,levels=0:6)
  d$visit <- if(part=="longitudinal") factor(d$visit,levels=unlist(cfg$visits)) else NULL
  scaling <- if(!is.null(reference)) reference$score_scaling else list()
  for (column in score_names) {
    if (is.null(reference)) {
      distinct <- !duplicated(d$subject_key)
      scaling[[column]] <- list(mean=mean(d[[column]][distinct]),sd=sd(d[[column]][distinct]))
    }
    item <- scaling[[column]]
    if(is.null(item) || !is.finite(item$sd) || item$sd<=0) return(clin_result(req,"not_estimable",paste("Non-varying score:",column),n_subjects=n))
    d[[column]] <- (d[[column]]-item$mean)/item$sd
  }
  d$log_icv <- log(d$icv_ml); d$log_wmh <- log1p(d$wmh_total_ml); d$log_gmv <- log(d$gmv_total_ml)
  if(b2) d$log_acute <- log1p(d$acute_volume_ml)
  # 可逆标准化改善阈值/斜率Hessian的数值尺度，不改变预设模型的拟合空间。
  if(!identical(cfg$continuous_covariate_standardization,"C3")) stop("Continuous covariate standardization must be configured as C3")
  covariate_scaling <- if(!is.null(reference)) reference$covariate_scaling else list()
  for(column in c("log_icv","log_wmh","log_gmv",if(b2) "log_acute")) {
    if(is.null(reference)) {
      distinct<-!duplicated(d$subject_key)
      covariate_scaling[[column]]<-list(mean=mean(d[[column]][distinct]),sd=sd(d[[column]][distinct]))
    }
    item<-covariate_scaling[[column]]
    if(is.null(item) || !is.finite(item$sd) || item$sd<=0)
      return(clin_result(req,"not_estimable",paste("Missing or non-varying frozen covariate scale:",column),n_subjects=n))
    d[[column]]<-(d[[column]]-item$mean)/item$sd
  }
  # 使用 C3 固定样条结点：GB 关联复现和预测单位保持相同。
  basis_specs <- if(!is.null(reference)) reference$basis_specs else list()
  basis_names <- list()
  for (variable in c("age","nihss",if(b2) "log_acute")) {
    df <- if(variable=="age") cfg$age_spline_df else if(variable=="nihss") cfg$nihss_spline_df else cfg$acute_volume_spline_df
    if(is.null(reference)) {
      basis <- splines::ns(d[[variable]],df=df)
      basis_specs[[variable]] <- list(knots=attr(basis,"knots"),Boundary.knots=attr(basis,"Boundary.knots"),intercept=FALSE)
    } else basis <- do.call(splines::ns,c(list(x=d[[variable]]),basis_specs[[variable]]))
    names <- paste0(variable,"_s",seq_len(ncol(basis)))
    d[names] <- as.data.frame(basis); basis_names[[variable]] <- names
  }
  if(predict_only) {
    label <- if(!is.null(req$validation_label)) req$validation_label else "frozen_model_prediction"
    predictions <- list();metrics <- list()
    for(name in names(reference$models)) {
      p<-clin_probabilities(reference$models[[name]],d)
      predictions[[name]]<-clin_prediction_frame(p,d,name,label)
      sc<-clin_scores(p,as.integer(d$mrs)-1L,cfg$prediction_epsilon)
      sc$model<-name;sc$validation<-label;metrics[[name]]<-sc
    }
    pred<-do.call(rbind,predictions)
    clin_table(pred,file.path(req$output_dir,"predictions.tsv"))
    clin_table(do.call(rbind,metrics),file.path(req$output_dir,"prediction_scores.tsv"))
    clin_table(clin_calibration(pred,cfg$calibration_bins),file.path(req$output_dir,"calibration.tsv"))
    return(clin_result(req,"completed","Frozen coefficients applied without heldout refitting",
      c("predictions.tsv","prediction_scores.tsv","calibration.tsv"),n))
  }
  base_terms <- c(basis_names$age,basis_names$nihss,setdiff(baseline,c("age","nihss")))
  total_terms <- c(base_terms,"log_icv","log_wmh","log_gmv")
  forms <- list()
  if(b2) {
    a_use <- clin_interaction_components(a,req$inputs$acute_component_metadata,cfg$max_acute_components)
    c_use <- c(clin_interaction_components(w,req$inputs$wmh_component_metadata,cfg$max_chronic_components),
               clin_interaction_components(g,req$inputs$gmv_component_metadata,cfg$max_chronic_components))
    # 完整近简并组可进入整体交互块，但不能单独赋予每个轴稳定的解剖含义。
    if(!length(a_use) || !length(c_use)) return(clin_result(req,"not_estimable","No identifiable stable interaction dimensions",n_subjects=n))
    location_terms <- as.vector(outer(a_use,c_use,paste,sep=":"))
    volume_terms <- paste("log_acute",c_use,sep=":")
    additive <- c(total_terms,basis_names$log_acute,a,w,g,volume_terms)
    forms <- list(additive=additive,location_interaction=c(additive,location_terms))
  } else if(part=="longitudinal") {
    terms <- c(total_terms,"visit",w,g)
    time_interactions<-if(length(c(w,g))) paste("visit",c(w,g),sep=":") else character()
    forms <- list(longitudinal_additive=terms,longitudinal_interaction=c(terms,time_interactions))
  } else forms <- list(M0=base_terms,M1=total_terms,MW=c(total_terms,w),MG=c(total_terms,g),MWG=c(total_terms,w,g))
  if(!is.null(reference)) forms <- reference$form_terms
  fits <- list(); effects <- list(); diagnostic <- list(); used <- list()
  for (name in names(forms)) {
    rhs <- paste(unique(forms[[name]]),collapse=" + ")
    f <- as.formula(paste("mrs ~",rhs))
    design <- model.matrix(f,d)
    if(qr(design)$rank<ncol(design)) return(clin_result(req,"not_estimable",paste("Rank deficient prespecified design:",name),n_subjects=n))
    if(n <= ncol(design)+cfg$minimum_residual_subjects) return(clin_result(req,"not_estimable",paste("Insufficient patients relative to model dimension:",name),n_subjects=n))
    if(part=="longitudinal") f <- as.formula(paste("mrs ~",rhs,"+ (1 | subject_key)"))
    identical_previous<-names(fits)[vapply(names(fits),function(previous)
      identical(unique(forms[[previous]]),unique(forms[[name]])),logical(1))]
    reused_from<-if(length(identical_previous)) identical_previous[1] else ""
    if(nzchar(reused_from)) fit<-fits[[reused_from]] else if(part=="longitudinal") {
      fit <- ordinal::clmm(f,data=d,link="logit",threshold="flexible",Hess=TRUE,nAGQ=cfg$nAGQ,control=clin_clmm_control(cfg))
    } else fit <- ordinal::clm(f,data=d,link="logit",threshold="flexible",Hess=TRUE)
    info <- tryCatch(clin_ready(fit,cfg$gradient_tolerance),error=function(e) {
      eigenvalues<-eigen(fit$Hessian,symmetric=TRUE,only.values=TRUE)$values
      clin_table(data.frame(model=name,optimizer=if(!is.null(fit$control$method)) fit$control$method else "Newton",
        optimizer_status=if(!is.null(fit$optRes$convergence)) fit$optRes$convergence else NA_integer_,
        max_gradient=if(!is.null(fit$gradient)) max(abs(fit$gradient)) else NA_real_,
        hessian_min_eigenvalue=min(eigenvalues),hessian_max_eigenvalue=max(eigenvalues),reason=conditionMessage(e)),
        file.path(req$output_dir,"diagnostic_failure.tsv"))
      saveRDS(fit,file.path(req$output_dir,"failed_model.rds"))
      stop(name,": ",conditionMessage(e),call.=FALSE)
    })
    fits[[name]] <- fit; used[[name]] <- f
    effects[[name]] <- clin_effects(fit,name,cfg$alpha)
    diagnostic[[name]] <- data.frame(model=name,n_subjects=n,n_records=nrow(d),parameters=info$df,
      reused_from=reused_from,
      log_likelihood=info$log_likelihood,max_gradient=info$max_gradient,
      hessian_min_eigenvalue=info$hessian_min_eigenvalue,hessian_condition=info$hessian_condition,
      random_intercept_sd=info$random_intercept_sd,
      estimand=if(part=="longitudinal") "subject_specific_conditional_OR" else "proportional_odds_OR")
  }
  effect <- do.call(rbind,effects); rownames(effect)<-NULL
  chosen <- if(b2) "location_interaction" else if(part=="longitudinal") "longitudinal_interaction" else "MWG"
  secondary <- effect$model==chosen & grepl(if(b2 || part=="longitudinal") ":" else "^(WMH_PC|GMV_PC)",effect$term)
  if(b2) secondary <- secondary & grepl("ACUTE_PC",effect$term)
  effect$q_secondary[secondary] <- p.adjust(effect$p[secondary],method="BH")
  tests <- if(b2) clin_lrt(fits$additive,fits$location_interaction,"acute_location_by_chronic_block") else if(part=="longitudinal")
    clin_lrt(fits$longitudinal_additive,fits$longitudinal_interaction,"visit_by_chronic_block") else
    rbind(clin_lrt(fits$M0,fits$M1,"global_volume_block"),clin_lrt(fits$M1,fits$MW,"WMH_spatial_block"),
      clin_lrt(fits$M1,fits$MG,"GMV_spatial_block"),clin_lrt(fits$M1,fits$MWG,"joint_spatial_block"))
  clin_table(effect,file.path(req$output_dir,"effects.tsv")); clin_table(tests,file.path(req$output_dir,"block_tests.tsv"))
  clin_table(do.call(rbind,diagnostic),file.path(req$output_dir,"diagnostics.tsv"))
  artifacts <- c("effects.tsv","block_tests.tsv","diagnostics.tsv","models.rds","analysis_subjects.tsv")
  clin_table(data.frame(subject_key=unique(d$subject_key)),file.path(req$output_dir,"analysis_subjects.tsv"))
  store <- list(models=fits,wmh_columns=w,gmv_columns=g,acute_columns=a,score_scaling=scaling,basis_specs=basis_specs,
                covariate_scaling=covariate_scaling,form_terms=forms,cohort=cohort,part=part,action=action,config=cfg,proportional_odds=list())
  if(part=="longitudinal") {
    sensitivity <- ordinal::clmm(used[[chosen]],data=d,link="logit",threshold="flexible",Hess=TRUE,nAGQ=cfg$validation_nAGQ,control=clin_clmm_control(cfg))
    clin_ready(sensitivity,cfg$gradient_tolerance)
    compare <- data.frame(term=names(fits[[chosen]]$beta),primary=fits[[chosen]]$beta,
                          higher_quadrature=sensitivity$beta,
                          absolute_difference=abs(fits[[chosen]]$beta-sensitivity$beta))
    clin_table(compare,file.path(req$output_dir,"quadrature_check.tsv"))
    observed <- as.data.frame(table(visit=d$visit,mrs=d$mrs))
    names(observed)[3] <- "n"
    clin_table(observed,file.path(req$output_dir,"observed_mrs_distribution.tsv"))
    timepoints<-clin_timepoint_effects(fits[[chosen]],c(w,g),cfg$visits,cfg$alpha)
    clin_table(timepoints,file.path(req$output_dir,"timepoint_effects.tsv"))
    artifacts <- c(artifacts,"quadrature_check.tsv","observed_mrs_distribution.tsv","timepoint_effects.tsv")
  } else {
    probs <- list(); scores <- list()
    for(name in names(fits)) {
      p <- clin_probabilities(fits[[name]],d)
      probs[[name]] <- clin_prediction_frame(p,d,name,paste0(cohort,"_association_apparent"))
      sc <- clin_scores(p,as.integer(d$mrs)-1L,cfg$prediction_epsilon)
      sc$model<-name; sc$validation<-paste0(cohort,"_association_apparent");scores[[name]]<-sc
    }
    if(!is.null(reference)) for(name in names(reference$models)) {
      p <- clin_probabilities(reference$models[[name]],d)
      probs[[paste0("frozen_",name)]] <- clin_prediction_frame(p,d,name,"GB_frozen_C3_external")
      sc <- clin_scores(p,as.integer(d$mrs)-1L,cfg$prediction_epsilon)
      sc$model<-name;sc$validation<-"GB_frozen_C3_external";scores[[paste0("frozen_",name)]]<-sc
    }
    pred <- do.call(rbind,probs)
    clin_table(pred,file.path(req$output_dir,"predictions.tsv")); clin_table(do.call(rbind,scores),file.path(req$output_dir,"prediction_scores.tsv"))
    bins <- cfg$calibration_bins
    clin_table(clin_calibration(pred,bins),file.path(req$output_dir,"calibration.tsv"))
    artifacts <- c(artifacts,"predictions.tsv","prediction_scores.tsv","calibration.tsv")
    nominal <- tryCatch(ordinal::nominal_test(fits[[chosen]]),error=function(e)e)
    if(inherits(nominal,"error")) {
      clin_table(data.frame(status="diagnostic_failed",reason=conditionMessage(nominal)),file.path(req$output_dir,"proportional_odds.tsv"))
    } else {
      tab <- as.data.frame(nominal);tab$term<-rownames(tab);rownames(tab)<-NULL
      clin_table(tab,file.path(req$output_dir,"proportional_odds.tsv"))
      pcol <- grep("Pr",names(tab),value=TRUE)[1]
      violates <- if(!is.na(pcol)) tab$term[!is.na(tab[[pcol]]) & tab[[pcol]]<cfg$nominal_alpha] else character()
      violates <- setdiff(violates,c("<none>","none"))
      store$proportional_odds <- list(violating_terms=violates)
      if(length(violates)) {
        partial <- tryCatch(ordinal::clm(used[[chosen]],nominal=as.formula(paste("~",paste(violates,collapse="+"))),
                              data=d,link="logit",Hess=TRUE),error=function(e)e)
        if(inherits(partial,"error")) clin_table(data.frame(status="not_estimable",reason=conditionMessage(partial)),file.path(req$output_dir,"partial_proportional_odds.tsv"))
        else {
          # 阈值特异效应保留原始系数，不误称为共同 OR。
          diag <- tryCatch({clin_ready(partial,cfg$gradient_tolerance);TRUE},error=function(e)e)
          if(inherits(diag,"error")) clin_table(data.frame(status="not_estimable",reason=conditionMessage(diag)),file.path(req$output_dir,"partial_proportional_odds.tsv"))
          else {saveRDS(partial,file.path(req$output_dir,"partial_proportional_odds.rds"));
            clin_table(data.frame(term=names(coef(partial)),estimate=coef(partial)),file.path(req$output_dir,"partial_proportional_odds.tsv"));artifacts<-c(artifacts,"partial_proportional_odds.rds")}
        }
        artifacts <- c(artifacts,"partial_proportional_odds.tsv")
      }
    }
    artifacts <- c(artifacts,"proportional_odds.tsv")
  }
  saveRDS(store,file.path(req$output_dir,"models.rds"))
  clin_result(req,"completed","Prespecified ordinal models fitted; diagnostics and apparent/external prediction labels retained",artifacts,n,
              list(n_records=nrow(d),n_complete=n,internal_validation="not_performed_in_this_task"))
}
