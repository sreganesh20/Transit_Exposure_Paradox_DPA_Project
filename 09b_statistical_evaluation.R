# Step 9B: prediction-only statistical evaluation. Run from repository root.
suppressPackageStartupMessages(library(data.table))
input <- Sys.getenv('ML_INPUT_FILE','data/processed/07a_coverage_aware_workloads_2024.csv')
src <- Sys.getenv('EVAL_9A_DIR','outputs/09a_ml_framework_validation')
out <- Sys.getenv('EVAL_9B_DIR','outputs/09b_statistical_evaluation_validation')
B <- as.integer(Sys.getenv('EVAL_BOOT_REPS','5000'))
stopifnot(!is.na(B), B >= 100L, file.exists(input), dir.exists(src))
stopifnot(normalizePath(src) != normalizePath(out,mustWork=FALSE))
if(dir.exists(out) && length(list.files(out))) stop('Choose an empty output directory to avoid stale exports.')
read9 <- function(n) fread(file.path(src,paste0(n,'.csv')))
manifest <- read9('export_manifest')
for(i in seq_len(nrow(manifest))) {
 x <- fread(file.path(src,manifest$file[i])); stopifnot(nrow(x)==manifest$rows[i],ncol(x)==manifest$columns[i])
}
meta <- read9('run_metadata'); stopifnot(nrow(meta)==1L,unname(tools::md5sum(input))==meta$input_md5)
d <- fread(input); p <- read9('final_predictions'); models <- c('LM-B','LM-A','GAM-B','GAM-A','RF-B','RF-A')
stopifnot(!anyDuplicated(d$query_id), !anyDuplicated(p[,.(query_id,model)]),setequal(p$model,models))
e <- d[split %in% c('UCI_spatial_holdout','CTA_transfer')]
stopifnot(setequal(p$query_id,e$query_id),nrow(p)==6L*nrow(e),all(p[, .N,by=query_id]$N==6L))
stopifnot(setequal(e$split,c('UCI_spatial_holdout','CTA_transfer')))
idx <- match(p$query_id,e$query_id)
stopifnot(all(p$workload==e$split[idx]),all(p$coverage_status==e$coverage_status[idx]),all(p$eligible_full_coverage==e$eligible_full_coverage[idx]),all(p$actual_count==e$crimes_2024[idx]))
close <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-9,check.attributes=FALSE)))
stopifnot(all(is.finite(p$predicted_log)),all(p$predicted_log>=0),all(is.finite(p$predicted_count)))
close(p$actual_log,log1p(p$actual_count)); close(p$predicted_count,expm1(p$predicted_log))
p[, `:=`(station_id=e$station_id[idx],area_km2=e$area_km2[idx],halfwidth=e$query_halfwidth_m[idx],block=paste(floor((e$center_lat[idx]-41)/.020),floor((e$center_lon[idx]+88)/.025),sep='_'))]
stopifnot(all(p[workload=='UCI_spatial_holdout',block]==p[workload=='UCI_spatial_holdout',spatial_block]))
p[,unused:=FALSE]; if('is_poc_query' %in% names(e)) p[,unused:=workload=='UCI_spatial_holdout' & !e$is_poc_query[idx]]
if ('is_poc_query' %in% names(e)) {
  expected_unused <- unique(p[unused==TRUE & eligible_full_coverage==TRUE,query_id])
  if (length(expected_unused) || file.exists(file.path(src,'unused_holdout_ids.csv'))) {
    ids <- read9('unused_holdout_ids')
    stopifnot(!anyDuplicated(ids$query_id),setequal(ids$query_id,expected_unused))
  }
}
s <- read9('support_diagnostics'); stopifnot(!anyDuplicated(s$query_id),setequal(s$query_id,e$query_id)); p[,outside:=s$augmented_outside[match(query_id,s$query_id)]]; stopifnot(!anyNA(p$outside))
metric <- function(a,l,c) { z<-log1p(a); r2<-function(x,y) if(sum((x-mean(x))^2)==0) NA_real_ else 1-sum((x-y)^2)/sum((x-mean(x))^2); list(n=length(a),log_mae=mean(abs(z-l)),log_rmse=sqrt(mean((z-l)^2)),log_r2=r2(z,l),count_mae=mean(abs(a-c)),count_rmse=sqrt(mean((a-c)^2)),count_r2=r2(a,c),relative_mae=if(mean(a)>0) mean(abs(a-c))/mean(a) else NA_real_,count_bias=mean(c-a),log_bias=mean(l-z)) }
q <- p[eligible_full_coverage==TRUE]; primary <- q[,metric(actual_count,predicted_log,predicted_count),by=.(workload,model)]
old <- read9('final_metrics'); setorder(old,workload,model); setorder(primary,workload,model); stopifnot(nrow(old)==nrow(primary)); for(n in names(old)) close(old[[n]],primary[[n]])
p[, `:=`(family=sub('-[AB]$','',model),config=substr(model,nchar(model),nchar(model)),log_error=abs(actual_log-predicted_log),count_error=abs(actual_count-predicted_count))]
paired <- dcast(p,workload+query_id+family+block+station_id+eligible_full_coverage+coverage_status+area_km2+halfwidth+outside+unused~config,value.var=c('log_error','count_error'))
paired[, `:=`(gain=log_error_B-log_error_A,count_gain=count_error_B-count_error_A)]
old <- read9('paired_errors'); chk <- paired[eligible_full_coverage==TRUE]; setorder(chk,workload,query_id,family); setorder(old,workload,query_id,model_family); close(chk$workload,old$workload); close(chk$query_id,old$query_id); close(chk$family,old$model_family); close(chk$block,old$spatial_block); close(chk$log_error_B,old$baseline); close(chk$log_error_A,old$augmented); close(chk$gain,old$error_improvement)
# Common strata across baseline and augmented models; partial coverage descriptive only.
strata <- function(x) {
 a<-x[eligible_full_coverage==TRUE]; a[,group:='eligible']; b<-copy(x); b[,group:=paste0('coverage:',coverage_status)]; c<-copy(a); c[,group:=ifelse(outside,'support:outside','support:inside')]; f<-copy(a[workload=='CTA_transfer']); f[,group:=paste0('halfwidth_m:',halfwidth)]; g<-copy(a[workload=='UCI_spatial_holdout']); g[,group:=as.character(cut(area_km2,c(-Inf,.1,1,10,Inf),include.lowest=TRUE))]; g[,group:=paste0('area_km2:',group)]; h<-copy(a[unused==TRUE]); h[,group:='previously_unused']; rbindlist(list(a,b,c,f,g,h),use.names=TRUE)
}
metrics <- strata(p)[,metric(actual_count,predicted_log,predicted_count),by=.(workload,group,model)]
comparisons <- strata(paired)[,.(n=.N,baseline_log_mae=mean(log_error_B),augmented_log_mae=mean(log_error_A),log_gain=mean(gain),log_gain_pct=if(mean(log_error_B)>0) 100*mean(gain)/mean(log_error_B) else NA_real_,count_gain=mean(count_gain),count_gain_pct=if(mean(count_error_B)>0) 100*mean(count_gain)/mean(count_error_B) else NA_real_,improved_pct=100*mean(gain>0),tied_pct=100*mean(gain==0)),by=.(workload,group,family)]
# Resample whole clusters; intervals conditional on fitted models, exploratory.
set.seed(meta$seed)
boot <- function(x,cluster,scheme) {
 x<-copy(x); x[,cluster_id:=get(cluster)]; stopifnot(!anyNA(x$cluster_id)); x[, { a<-.SD[,.(n=.N,g=mean(gain)),by=cluster_id]; k<-nrow(a); draws<-replicate(B,{j<-sample.int(k,k,replace=TRUE);c(weighted.mean(a$g[j],a$n[j]),mean(a$g[j]))}); list(n_clusters=k,n_queries=sum(a$n),query_weighted=weighted.mean(a$g,a$n),block_weighted=mean(a$g),query_lower=if(k>1) quantile(draws[1,],.025) else NA_real_,query_upper=if(k>1) quantile(draws[1,],.975) else NA_real_,equal_lower=if(k>1) quantile(draws[2,],.025) else NA_real_,equal_upper=if(k>1) quantile(draws[2,],.975) else NA_real_,scheme=scheme) },by=.(workload,family)]
}
intervals <- rbindlist(list(boot(chk,'block','geographic_block'),boot(chk[workload=='CTA_transfer'],'station_id','station_cluster')))
# Compare existing block estimates, not independently seeded interval endpoints.
old <- read9('bootstrap_results'); v<-intervals[scheme=='geographic_block']; setorder(v,workload,family); setorder(old,workload,model_family); close(v$workload,old$workload); close(v$family,old$model_family); close(v$query_weighted,old$observed_improvement); close(v$n_clusters,old$n_blocks)
stopifnot(nrow(primary)==12L,nrow(intervals)==9L,
          all(is.finite(paired$gain)),all(is.finite(paired$count_gain)),
          nrow(chk)==3L*uniqueN(q$query_id))
script_arg <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
script_md5 <- if(length(script_arg)==1L) unname(tools::md5sum(sub('^--file=','',script_arg))) else NA_character_
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tables<-list(primary_metrics=primary,stratified_metrics=metrics,paired_comparisons=comparisons,cluster_intervals=intervals,run_metadata=data.table(input_md5=meta$input_md5,source_dir=src,bootstrap_reps=B,seed=meta$seed,script_md5=script_md5,R_version=R.version.string,data_table_version=as.character(packageVersion("data.table")),source_manifest_md5=unname(tools::md5sum(file.path(src,"export_manifest.csv")))))
for(n in names(tables)) fwrite(tables[[n]],file.path(out,paste0(n,'.csv')))
png(file.path(out, 'mobility_improvement.png'),
    width = 1200, height = 800, res = 120)

z <- comparisons[group == 'eligible']
m <- dcast(z, family ~ workload, value.var = 'log_gain_pct')
workloads <- c('CTA_transfer','UCI_spatial_holdout')
stopifnot(setequal(names(m)[-1],workloads))
vals <- t(as.matrix(m[, ..workloads]))
limits <- range(c(0,vals))
span <- max(diff(limits),1)

par(mar = c(5, 5, 5, 2))

barplot(vals,
        beside = TRUE,
        names.arg = m$family,
        ylim = limits + c(-0.08,0.30)*span,
        ylab = 'Log MAE improvement (%)',
        col = c('steelblue', 'darkorange'))

abline(h = 0)

legend('top', inset = 0.02,
       legend = c('CTA transfer', 'UCI spatial holdout'),
       fill = c('steelblue', 'darkorange'),
       horiz = TRUE, bty = 'n')

dev.off()
# Generate descriptive findings from the current run, not hardcoded results.
findings <- comparisons[group=='eligible']
setorder(findings,workload,family)
findings_text <- with(findings,sprintf('%s / %s: log MAE gain %.6f (%.2f%%); count MAE gain %.3f (%.2f%%); %.2f%% of queries improved in absolute log error.',workload,family,log_gain,log_gain_pct,count_gain,count_gain_pct,improved_pct))
best <- primary[order(log_mae),.SD[1L],by=workload]
best_text <- with(best,sprintf('%s: lowest observed log MAE is %s (%.6f).',workload,model,log_mae))
interval_text <- with(intervals,sprintf('%s / %s / %s: %d clusters; query-weighted gain %.6f, exploratory 95%% percentile interval [%.6f, %.6f]; equal-cluster gain %.6f, interval [%.6f, %.6f].',workload,family,scheme,n_clusters,query_weighted,query_lower,query_upper,block_weighted,equal_lower,equal_upper))
unused_text <- if(any(findings$workload=='UCI_spatial_holdout') && any(comparisons$group=='previously_unused')) {
  u <- comparisons[group=='previously_unused']
  c('Previously unused eligible holdout is a subset of the full holdout, not independent replication.',with(u,sprintf('%s: %d queries; log MAE gain %.2f%%.',family,n,log_gain_pct)))
} else 'No previously unused eligible holdout subset reported.'
writeLines(c('Primary analysis uses eligible full coverage only. Positive gains favor augmentation. Relative MAE = count MAE / mean observed count. Bias = predicted minus observed.',
 'Coverage strata outside full coverage are descriptive, not comparable full-city outcomes. Support flags are marginal percentile bounds, not joint support.',
 'Intervals are exploratory 95% percentile cluster bootstrap intervals conditional on the fitted models. They do not capture model fitting uncertainty. Spatial overlap between clusters remains possible. Few clusters make intervals especially uncertain. No formal significance or causal claims.',
 'Area bands are fixed descriptive thresholds in square kilometers. Station clustering is an alternative sensitivity analysis, not independent confirmation.',
 '', 'Observed primary comparisons:',findings_text,best_text,'','Cluster sensitivity:',interval_text,'',unused_text),file.path(out,'interpretation.txt'))
cat('9B checks passed; outputs:',out,'\n')
