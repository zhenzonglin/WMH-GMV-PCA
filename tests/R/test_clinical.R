# Run from repository root: Rscript tests/R/test_clinical.R
source("R/clinical.R")
if (!requireNamespace("ordinal",quietly=TRUE)) stop("Install the pinned ordinal dependency before testing")
set.seed(20261030)
cfg <- list(analysis=list(age_range=c(50,75),clinical=list(visits=c(3,6,12),age_spline_df=3,nihss_spline_df=3,
    acute_volume_spline_df=3,max_acute_components=2,max_chronic_components=2,nAGQ=7,validation_nAGQ=11,
    alpha=0.05,nominal_alpha=0.01,prediction_epsilon=1e-8,exclude_conflicts=TRUE,
    min_complete_cases=30,calibration_bins=10,gradient_tolerance=0.001,minimum_residual_subjects=10,continuous_covariate_standardization="C3",
    covariates=c("age","sex","prestroke_mrs","nihss","prior_stroke","hypertension","diabetes","af","ivt","evt"),
    optimizer=list(method="nlminb",relative_tolerance=1e-10,parameter_tolerance=1e-8,
      inner_gradient_tolerance=1e-4,outer_gradient_tolerance=1e-4,max_iterations=1000,max_evaluations=2000))))
test_optimizer<-Sys.getenv("WMH_TEST_OPTIMIZER","ucminf")
cfg$analysis$clinical$optimizer$method<-test_optimizer
if(test_optimizer=="ucminf") cfg$analysis$clinical$optimizer$parameter_tolerance<-1e-12
test_inner<-Sys.getenv("WMH_TEST_INNER_GRADIENT","")
if(nzchar(test_inner)) cfg$analysis$clinical$optimizer$inner_gradient_tolerance<-as.numeric(test_inner)
test_outer<-Sys.getenv("WMH_TEST_OUTER_GRADIENT","")
if(nzchar(test_outer)) cfg$analysis$clinical$optimizer$outer_gradient_tolerance<-as.numeric(test_outer)
make_data <- function(n,cohort) {
  d <- data.frame(subject_key=paste0(cohort,":",sprintf("%05d",seq_len(n))),cohort=cohort,
    age=runif(n,50,75),sex=rbinom(n,1,.5),prestroke_mrs=sample(0:3,n,TRUE),nihss=rpois(n,6),
    prior_stroke=rbinom(n,1,.2),hypertension=rbinom(n,1,.6),diabetes=rbinom(n,1,.3),
    af=rbinom(n,1,.2),ivt=rbinom(n,1,.3),evt=rbinom(n,1,.2),icv_ml=rnorm(n,1400,100),
    wmh_total_ml=rgamma(n,2,.3),gmv_total_ml=rnorm(n,600,40),WMH_PC01=rnorm(n),GMV_PC01=rnorm(n),
    ACUTE_PC01=rnorm(n),acute_volume_ml=rgamma(n,2,.1),clinical_conflict=0)
  d$mrs_m12 <- as.integer(cut(.6*d$WMH_PC01+.2*d$GMV_PC01+rlogis(n),
                              breaks=c(-Inf,-2,-1,-.3,.3,1,2,Inf),labels=FALSE))-1L
  d
}
work <- Sys.getenv("WMH_TEST_OUTPUT",tempfile("wmh-clinical-test-"))
dir.create(work,recursive=TRUE,showWarnings=FALSE)
d <- make_data(500,"C3")
input <- file.path(work,"C3.tsv");clin_table(d,input)
req <- list(action="clinical_b1_single",config=cfg,inputs=list(data=input),output_dir=file.path(work,"c3_b1"),cohort="C3")
result <- run_clinical(req)
stopifnot(result$status=="completed")
model <- readRDS(file.path(req$output_dir,"models.rds"))
effects <- read.delim(file.path(req$output_dir,"effects.tsv"))
stopifnot(effects$OR[effects$model=="MWG" & effects$term=="WMH_PC01"]>1)
# 可逆总量标准化不改变拟合空间：独立重拟合未标准化设计并比较全部类别概率。
scaled_design<-model$models$MWG$model
unscaled_design<-scaled_design
for(column in names(model$covariate_scaling)) {
  item<-model$covariate_scaling[[column]]
  unscaled_design[[column]]<-unscaled_design[[column]]*item$sd+item$mean
}
unscaled_fit<-ordinal::clm(formula(model$models$MWG),data=unscaled_design,link="logit",Hess=TRUE)
scaled_probability<-predict(model$models$MWG,newdata=scaled_design[,setdiff(names(scaled_design),"mrs"),drop=FALSE],type="prob")$fit
unscaled_probability<-predict(unscaled_fit,newdata=unscaled_design[,setdiff(names(unscaled_design),"mrs"),drop=FALSE],type="prob")$fit
equivalence_error<-max(abs(scaled_probability-unscaled_probability))
cat("Covariate reparameterization maximum probability difference:",equivalence_error,"\n")
stopifnot(equivalence_error<1e-7)
gb <- make_data(400,"GB");gb_input<-file.path(work,"GB.tsv");clin_table(gb,gb_input)
gb_req<-req;gb_req$cohort<-"GB";gb_req$inputs<-list(data=gb_input,reference_models=file.path(req$output_dir,"models.rds"));gb_req$output_dir<-file.path(work,"gb_b1")
stopifnot(run_clinical(gb_req)$status=="completed")
gb_model <- readRDS(file.path(gb_req$output_dir,"models.rds"))
stopifnot(identical(model$score_scaling,gb_model$score_scaling),identical(model$basis_specs,gb_model$basis_specs),identical(model$covariate_scaling,gb_model$covariate_scaling))
p <- read.delim(file.path(gb_req$output_dir,"predictions.tsv"))
stopifnot("GB_frozen_C3_external" %in% p$validation,max(abs(rowSums(p[paste0("p_mrs",0:6)])-1))<1e-6)
b2<-req;b2$action<-"clinical_b2";b2$output_dir<-file.path(work,"c3_b2")
for(modality in c("wmh","gmv","acute")) {
  metadata_path<-file.path(work,paste0(modality,"_components.tsv"))
  clin_table(data.frame(component=paste0(toupper(modality),"_PC01"),index=1,group_id=1,
              group_size=1,individually_identifiable=TRUE),metadata_path)
  b2$inputs[[paste0(modality,"_component_metadata")]]<-metadata_path
}
stopifnot(run_clinical(b2)$status=="completed")
blocks<-read.delim(file.path(b2$output_dir,"block_tests.tsv"))
stopifnot(blocks$test=="acute_location_by_chronic_block",blocks$df==2)
near_path<-file.path(work,"near.tsv")
clin_table(data.frame(component=c("WMH_PC01","WMH_PC02","WMH_PC03"),index=1:3,group_id=c(1,2,2),
  group_size=c(1,2,2),individually_identifiable=c(TRUE,FALSE,FALSE)),near_path)
stopifnot(identical(clin_interaction_components(c("WMH_PC01","WMH_PC02","WMH_PC03"),near_path,2),"WMH_PC01"))
patient_intercept <- rnorm(180,0,1.2)
long <- do.call(rbind,lapply(c(3,6,12),function(v) {
  z<-d[seq_len(180),];z$visit<-v
  z$mrs<-as.integer(cut(patient_intercept+.6*z$WMH_PC01+.2*z$GMV_PC01-.025*v+rlogis(nrow(z)),
                         breaks=c(-Inf,-2,-1,-.3,.3,1,2,Inf),labels=FALSE))-1L
  z
}))
long_path<-file.path(work,"long.tsv");clin_table(long,long_path)
long_req<-req;long_req$action<-"clinical_b1_longitudinal";long_req$inputs$long<-long_path
long_req$output_dir<-file.path(work,"c3_long")
stopifnot(run_clinical(long_req)$status=="completed")
stopifnot(file.exists(file.path(long_req$output_dir,"quadrature_check.tsv")))
timepoints<-read.delim(file.path(long_req$output_dir,"timepoint_effects.tsv"))
long_fit<-readRDS(file.path(long_req$output_dir,"models.rds"))$models$longitudinal_interaction
beta<-long_fit$beta;covariance<-vcov(long_fit)
at3<-timepoints[timepoints$component=="WMH_PC01" & timepoints$visit==3,]
at6<-timepoints[timepoints$component=="WMH_PC01" & timepoints$visit==6,]
interaction<-intersect(c("visit6:WMH_PC01","WMH_PC01:visit6"),names(beta))
expected_variance<-covariance["WMH_PC01","WMH_PC01"]+covariance[interaction,interaction]+2*covariance["WMH_PC01",interaction]
stopifnot(abs(at3$log_OR-beta[["WMH_PC01"]])<1e-10,
          abs(at6$log_OR-beta[["WMH_PC01"]]-beta[[interaction]])<1e-10,
          abs(at6$se^2-expected_variance)<1e-10,
          abs(log(at6$lower)-(at6$log_OR-qnorm(1-cfg$analysis$clinical$alpha/2)*sqrt(expected_variance)))<1e-10,
          !any(c("p","q_secondary") %in% names(timepoints)),
          nrow(clin_timepoint_effects(long_fit,character(),c(3,6,12),.05))==0L)
# K=0 是合法完整流程：没有伪造 visit:，并复用相同加性拟合。
zero_long<-long[,!grepl("^(WMH_PC|GMV_PC)",names(long)),drop=FALSE]
zero_long_path<-file.path(work,"zero_long.tsv");clin_table(zero_long,zero_long_path)
zero_req<-long_req;zero_req$inputs$long<-zero_long_path;zero_req$output_dir<-file.path(work,"c3_zero_long")
stopifnot(run_clinical(zero_req)$status=="completed")
zero_models<-readRDS(file.path(zero_req$output_dir,"models.rds"))
zero_tests<-read.delim(file.path(zero_req$output_dir,"block_tests.tsv"))
zero_diagnostics<-read.delim(file.path(zero_req$output_dir,"diagnostics.tsv"))
zero_timepoints<-read.delim(file.path(zero_req$output_dir,"timepoint_effects.tsv"))
stopifnot(identical(zero_models$models$longitudinal_additive,zero_models$models$longitudinal_interaction),
  zero_tests$status=="not_estimable_no_added_dimension",zero_tests$df==0,is.na(zero_tests$p),
  zero_diagnostics$reused_from[zero_diagnostics$model=="longitudinal_interaction"]=="longitudinal_additive",
  nrow(zero_timepoints)==0L)
# 缺少变量与缺少结局等级应标记不可执行，不能换一个更简单模型。
bad<-d;bad$af<-NULL;bad_path<-file.path(work,"bad.tsv");clin_table(bad,bad_path)
bad_req<-req;bad_req$inputs$data<-bad_path;bad_req$output_dir<-file.path(work,"bad")
stopifnot(run_clinical(bad_req)$status=="waiting_input")
emptycat<-d;emptycat$mrs_m12[emptycat$mrs_m12==6]<-5
empty_path<-file.path(work,"emptycat.tsv");clin_table(emptycat,empty_path)
bad_req$inputs$data<-empty_path;bad_req$output_dir<-file.path(work,"emptycat")
stopifnot(run_clinical(bad_req)$status=="not_estimable")
cat("Clinical R single-point, interaction, frozen-GB and failure-mode tests passed.\n")
cat("Test artifacts: ",work,"\n",sep="")
