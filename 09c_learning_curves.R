# Step 9C - fixed-hyperparameter, nested training-size sensitivity.
# Run from repository root. Source 9A helpers only (no tuning or 9A rerun).
# Configuration: ML_INPUT_FILE, LC_9A_DIR, LC_OUTPUT_DIR, LC_SIZES,
# LC_REPS, LC_REUSE_FULL, GAM_ENGINE, GAM_NUM_THREADS, RF_NUM_THREADS.
# Checkpoints are deliberately separate for smoke and full-scale runs.
main_9c <- function() {
  stopifnot(file.exists('09a_model_helpers.R'), file.exists('09c_learning_curves.R'))
  source('09a_model_helpers.R', local = TRUE)
  input <- normalizePath(Sys.getenv('ML_INPUT_FILE', 'data/processed/08_full_workload_2024.csv'), mustWork = TRUE)
  src <- normalizePath(Sys.getenv('LC_9A_DIR', 'outputs/09a_ml_framework_full'), mustWork = TRUE)
  out <- normalizePath(Sys.getenv('LC_OUTPUT_DIR', 'outputs/09c_learning_curves_full'), mustWork = FALSE)
  stopifnot(src != out)
  read9 <- function(n) fread(file.path(src, paste0(n, '.csv')))
  hash <- function(path) { h <- unname(tools::md5sum(path)); stopifnot(!anyNA(h)); h }
  same <- function(a, b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-9,check.attributes=FALSE)))
  manifest <- read9('export_manifest'); meta <- read9('run_metadata')
  stopifnot(nrow(meta)==1L, hash(input)==meta$input_md5, meta$seed==SEED, meta$n_folds==N_FOLDS,
            !anyDuplicated(manifest$file), all(file.exists(file.path(src,manifest$file))))
  needed <- c('selected_hyperparameters','spatial_folds','final_predictions','final_metrics','final_fit_summary',
              'gam_baseline_tuning','gam_augmented_tuning','rf_baseline_tuning','rf_augmented_tuning','run_metadata')
  stopifnot(all(paste0(needed,'.csv') %in% manifest$file))
  d <- prepare_features(fread(input))
  train <- d[split=='UCI_train' & eligible_full_coverage==TRUE]
  hold <- d[split=='UCI_spatial_holdout' & eligible_full_coverage==TRUE]
  cta <- d[split=='CTA_transfer' & eligible_full_coverage==TRUE]
  N <- nrow(train); blocks <- sort(unique(train$spatial_block)); K <- length(blocks)
  stopifnot(N==meta$training_queries, nrow(hold)>0L,nrow(cta)>0L,
            !anyNA(train$spatial_block), !anyNA(hold$spatial_block), K>=N_FOLDS,
            !length(intersect(train$spatial_block,hold$spatial_block)))
  evals <- list(UCI_spatial_holdout=hold, CTA_transfer=cta)
  if ('is_poc_query' %in% names(hold)) {
    stopifnot(!anyNA(hold$is_poc_query))
    u <- hold[is_poc_query==FALSE]
    if(nrow(u)) {
      ids <- read9('unused_holdout_ids')
      stopifnot(!anyDuplicated(ids$query_id),setequal(ids$query_id,u$query_id))
      evals$UCI_previously_unused <- u
    }
  }
  rm(d)
  folds <- read9('spatial_folds')
  cnt <- train[,.(n_queries=.N),by=spatial_block]
  stopifnot(!anyDuplicated(folds$spatial_block),setequal(folds$spatial_block,cnt$spatial_block),
            setequal(folds$fold_id,seq_len(N_FOLDS)))
  same(folds$n_queries,cnt$n_queries[match(folds$spatial_block,cnt$spatial_block)])
  models <- c('LM-B','LM-A','GAM-B','GAM-A','RF-B','RF-A')
  hyper <- read9('selected_hyperparameters')
  stopifnot(nrow(hyper)==6L,setequal(hyper$model,models),!anyDuplicated(hyper$model))
  hyper <- hyper[match(models,model)]
  for(id in models) {
    hp <- hyper[model==id]
    tag <- if(endsWith(id,'-B')) 'baseline' else 'augmented'
    if(startsWith(id,'GAM')) {
      t <- read9(paste0('gam_',tag,'_tuning'))
      setorder(t,log_mae,k)
      stopifnot(hp$k==t$k[1L],hp$k %in% GAM_K_GRID,hp$gam_engine==GAM_ENGINE)
    }
    if(startsWith(id,'RF')) {
      t <- read9(paste0('rf_',tag,'_tuning'))
      setorder(t,log_mae,-min_node_size,mtry)
      stopifnot(hp$mtry==t$mtry[1L],hp$min_node_size==t$min_node_size[1L],
                hp$num_trees==RF_NUM_TREES, hp$num_trees==t$num_trees[1L])
    }
  }
  stopifnot(GAM_ENGINE %in% c('bam','gam')) # fixed engine across all sample sizes
  # Filter default sizes to the available workload; explicit invalid sizes fail.
  default_sizes <- c(565L,1000L,5000L,20000L,50000L,100000L)
  raw_sizes <- Sys.getenv('LC_SIZES','')
  sizes <- if(nzchar(raw_sizes)) suppressWarnings(as.numeric(strsplit(raw_sizes,',',fixed=TRUE)[[1L]])) else default_sizes[default_sizes<N]
  stopifnot(length(sizes)>0L,all(is.finite(sizes)),all(sizes==floor(sizes)),all(sizes>=K),all(sizes<N))
  sizes <- sort(unique(as.integer(sizes)))
  reps <- suppressWarnings(as.numeric(Sys.getenv('LC_REPS','3')))
  stopifnot(length(reps)==1L,is.finite(reps),reps==floor(reps),reps>=1,reps<=100)
  reps <- as.integer(reps)
  reuse_flag <- Sys.getenv('LC_REUSE_FULL','TRUE'); stopifnot(reuse_flag %in% c('TRUE','FALSE'))
  reuse_full <- reuse_flag=='TRUE'
  seeds <- as.integer(SEED+seq_len(reps)-1L)
  # Verify the full-size reference using the received predictions and 9A metrics.
  ref <- read9('final_predictions'); ref <- ref[eligible_full_coverage==TRUE]
  eval_ids <- c(hold$query_id,cta$query_id)
  stopifnot(nrow(ref)==6L*length(eval_ids),setequal(ref$query_id,eval_ids),setequal(ref$model,models),
            !anyDuplicated(ref[,.(query_id,model)]),all(is.finite(ref$predicted_log)),all(ref$predicted_log>=0),
            all(is.finite(ref$predicted_count)))
  same(ref$predicted_count,expm1(ref$predicted_log))
  fit_summary <- read9('final_fit_summary')
  stopifnot(nrow(fit_summary)==6L,setequal(fit_summary$model,models),all(fit_summary$training_queries==N))
  get_reference <- function(id,work) {
    ev <- evals[[work]]
    src_work <- if(work=='UCI_previously_unused') 'UCI_spatial_holdout' else work
    p <- ref[model==id & workload==src_work & query_id %in% ev$query_id]
    stopifnot(nrow(p)==nrow(ev),setequal(p$query_id,ev$query_id))
    p <- p[match(ev$query_id,query_id)]
    same(p$actual_count,ev$crimes_2024); same(p$actual_log,ev$log_crimes)
    evaluate_predictions(ev$crimes_2024,p$predicted_log)
  }
  reference <- rbindlist(lapply(models,function(id) rbindlist(lapply(names(evals),function(w) {
    m <- get_reference(id,w); m[,`:=`(model=id,workload=w)]; m
  }))))
  original <- read9('final_metrics'); chk <- reference[workload!='UCI_previously_unused']
  setorder(chk,workload,model);setorder(original,workload,model)
  stopifnot(nrow(chk)==nrow(original)); for(nm in names(original)) same(chk[[nm]],original[[nm]])
  if('UCI_previously_unused' %in% names(evals)) {
    old_u <- read9('unused_holdout_metrics'); new_u <- reference[workload=='UCI_previously_unused'];setorder(old_u,model);setorder(new_u,model)
    stopifnot(nrow(old_u)==nrow(new_u));for(nm in setdiff(names(old_u),'workload')) same(old_u[[nm]],new_u[[nm]])
  }
  rm(ref,chk,original)
  # Checkpoint signatures include data, exact source code, 9A files, runtime and design.
  sig_files <- c(input,'09a_model_helpers.R','09c_learning_curves.R',file.path(src,c('export_manifest.csv',paste0(needed,'.csv'))))
  if('UCI_previously_unused' %in% names(evals)) sig_files <- c(sig_files,file.path(src,c('unused_holdout_ids.csv','unused_holdout_metrics.csv')))
  signature <- list(version='9c-v1',files=setNames(hash(sig_files),sig_files),sizes=sizes,reps=reps,seeds=seeds,reuse_full=reuse_full,
                    engine=GAM_ENGINE,gam_threads=GAM_NUM_THREADS,rf_threads=RF_NUM_THREADS,trees=RF_NUM_TREES,model_seed=SEED,
                    R=R.version.string,RNG=RNGkind(),packages=sapply(c('data.table','mgcv','ranger'),function(p)as.character(packageVersion(p))))
  ckdir <- file.path(out,'checkpoints'); cfg <- file.path(ckdir,'configuration.rds')
  atomic_rds <- function(x,path) {
    tmp <- tempfile(pattern='.tmp-',tmpdir=dirname(path)); on.exit(unlink(tmp),add=TRUE)
    saveRDS(x,tmp); stopifnot(!file.exists(path),file.rename(tmp,path))
  }
  if(file.exists(cfg)) {
    if(!identical(readRDS(cfg),signature)) stop('Checkpoint configuration differs; choose a new LC_OUTPUT_DIR.')
  } else {
    if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop('Output directory is not empty; choose a new LC_OUTPUT_DIR.')
    dir.create(ckdir,recursive=TRUE,showWarnings=FALSE); atomic_rds(signature,cfg)
  }
  nested <- function(seed) {
    set.seed(seed)
    groups <- split(seq_len(N),train$spatial_block)
    anchors <- vapply(groups,function(i)i[sample.int(length(i),1L)],integer(1))
    rest <- setdiff(seq_len(N),anchors)
    c(unname(anchors),rest[sample.int(length(rest))])
  }
  sample_hash <- function(ids) {
    f <- tempfile();on.exit(unlink(f),add=TRUE);writeLines(as.character(ids),f,useBytes=TRUE);hash(f)
  }
  metric_names <- c('log_mae','log_rmse','log_r2','count_mae','count_rmse','count_r2')
  evaluate_fit <- function(id,subtrain) {
    hp <- hyper[model==id]
    features <- if(endsWith(id,'-B')) baseline_features else augmented_features
    fd <- as.data.frame(subtrain[,c(features,target),with=FALSE])
    family <- sub('-[BA]$','',id)
    fit_seconds <- system.time({
      fit <- switch(family,LM=fit_lm(fd,features),GAM=fit_gam(fd,features,hp$k),RF=fit_rf(fd,features,hp$mtry,hp$min_node_size))
    })[['elapsed']]
    predict_fn <- switch(family,LM=predict_lm,GAM=predict_gam,RF=predict_rf)
    pred_seconds <- system.time({
      preds <- lapply(evals[c('UCI_spatial_holdout','CTA_transfer')],function(ev) {
        lp <- as.numeric(predict_fn(fit,as.data.frame(ev[,features,with=FALSE])))
        stopifnot(length(lp)==nrow(ev),all(is.finite(lp)))
        lp
      })
      if('UCI_previously_unused' %in% names(evals)) {
        ix <- match(evals$UCI_previously_unused$query_id,evals$UCI_spatial_holdout$query_id)
        stopifnot(!anyNA(ix))
        preds$UCI_previously_unused <- preds$UCI_spatial_holdout[ix]
      }
      res <- rbindlist(lapply(names(evals),function(work) {
        m <- evaluate_predictions(evals[[work]]$crimes_2024,preds[[work]])
        m[,workload:=work];m
      }))
    })[['elapsed']]
    res[,`:=`(fit_seconds=fit_seconds,prediction_evaluation_seconds=pred_seconds)]
    res
  }
  check_fit <- function(obj,key) {
    stopifnot(identical(obj$key,key),is.data.table(obj$metrics),nrow(obj$metrics)==length(evals),
              setequal(obj$metrics$workload,names(evals)),!anyDuplicated(obj$metrics$workload),
              all(obj$metrics$n==vapply(evals,nrow,integer(1))[obj$metrics$workload]),
              all(vapply(obj$metrics[,..metric_names],function(x)all(is.finite(x)),logical(1))),
              all(is.finite(obj$metrics$fit_seconds)),all(obj$metrics$fit_seconds>=0),
              all(is.finite(obj$metrics$prediction_evaluation_seconds)),all(obj$metrics$prediction_evaluation_seconds>=0))
  }
  results <- list(); sampling <- list(); block_sampling <- list(); timings <- list()
  fitted_now <- 0L; resumed <- 0L; total <- reps*length(sizes)*length(models); job <- 0L
  for(r in seq_len(reps)) {
    ord <- nested(seeds[r]); stopifnot(length(ord)==N,!anyDuplicated(ord),identical(ord,nested(seeds[r])))
    for(size_n in sizes) {
      subtrain <- train[ord[seq_len(size_n)]]
      stopifnot(nrow(subtrain)==size_n,!anyDuplicated(subtrain$query_id),uniqueN(subtrain$spatial_block)==K)
      sh <- sample_hash(subtrain$query_id)
      sampling[[length(sampling)+1L]] <- data.table(repetition=r,sampling_seed=seeds[r],training_size=size_n,spatial_blocks=K,sample_md5=sh)
      bc <- subtrain[,.(sample_queries=.N),by=spatial_block]
      bc <- merge(bc,folds,by='spatial_block',all.x=TRUE);bc[,`:=`(repetition=r,training_size=size_n)]
      block_sampling[[length(block_sampling)+1L]] <- bc
      for(id in models) {
        job <- job+1L
        key <- list(repetition=r,seed=seeds[r],size=size_n,model=id,sample_md5=sh)
        path <- file.path(ckdir,sprintf('r%03d_n%07d_%s.rds',r,size_n,id))
        if(file.exists(path)) {
          obj <- readRDS(path);check_fit(obj,key);resumed <- resumed+1L
          cat(sprintf('[%d/%d] RESUMED r=%d n=%d %s\n',job,total,r,size_n,id))
        } else {
          cat(sprintf('[%d/%d] FIT r=%d n=%d %s\n',job,total,r,size_n,id));flush.console()
          obj <- list(key=key,metrics=evaluate_fit(id,subtrain));check_fit(obj,key)
          atomic_rds(obj,path);fitted_now <- fitted_now+1L
          cat(sprintf('  saved; fit %.2f s, evaluation %.2f s\n',obj$metrics$fit_seconds[1],obj$metrics$prediction_evaluation_seconds[1]))
          gc(verbose=FALSE)
        }
        m <- copy(obj$metrics);m[,`:=`(model=id,repetition=r,sampling_seed=seeds[r],training_size=size_n,sample_md5=sh,origin='nested_sample')]
        results[[length(results)+1L]] <- m
        timings[[length(timings)+1L]] <- m[1L,.(model,repetition,training_size,fit_seconds,prediction_evaluation_seconds)]
      }
    }
  }
  if(reuse_full) {
    reference[,`:=`(repetition=0L,sampling_seed=NA_integer_,training_size=N,sample_md5=sample_hash(train$query_id),origin='9A_full_reference',fit_seconds=NA_real_,prediction_evaluation_seconds=NA_real_)]
    results[[length(results)+1L]] <- reference
  }
  raw <- rbindlist(results,use.names=TRUE);setorder(raw,workload,model,training_size,repetition)
  stopifnot(nrow(raw)==(total+(if(reuse_full) 6L else 0L))*length(evals),!anyDuplicated(raw[,.(workload,model,training_size,repetition)]))
  summary <- raw[,c(list(n_repetitions=.N),setNames(lapply(.SD,mean),paste0(metric_names,'_mean')),setNames(lapply(.SD,sd),paste0(metric_names,'_sd'))),by=.(workload,model,training_size,origin),.SDcols=metric_names]
  raw[,`:=`(family=sub('-[BA]$','',model),configuration=substr(model,nchar(model),nchar(model)))]
  gains <- dcast(raw,workload+family+training_size+repetition+origin+sample_md5~configuration,value.var=c('log_mae','count_mae'))
  stopifnot(!anyNA(gains$log_mae_A),!anyNA(gains$log_mae_B),all(gains$log_mae_B>0))
  gains[,`:=`(log_gain=log_mae_B-log_mae_A,log_gain_pct=100*(log_mae_B-log_mae_A)/log_mae_B,count_gain=count_mae_B-count_mae_A)]
  gain_summary <- gains[,.(n_repetitions=.N,mean_log_gain=mean(log_gain),sd_log_gain=sd(log_gain),mean_log_gain_pct=mean(log_gain_pct),sd_log_gain_pct=sd(log_gain_pct),mean_count_gain=mean(count_gain)),by=.(workload,family,training_size,origin)]
  tables <- list(learning_curve_runs=raw,learning_curve_summary=summary,paired_mobility_gains=gains,mobility_gain_summary=gain_summary,
                 sampling_summary=rbindlist(sampling),sampling_by_block=rbindlist(block_sampling),fit_runtimes=rbindlist(timings),
                 selected_hyperparameters=hyper,run_metadata=data.table(input_md5=meta$input_md5,source_manifest_md5=hash(file.path(src,'export_manifest.csv')),helper_md5=hash('09a_model_helpers.R'),script_md5=hash('09c_learning_curves.R'),seed=SEED,repetitions=reps,sizes=paste(sizes,collapse=','),reuse_full=reuse_full,gam_engine=GAM_ENGINE,gam_threads=GAM_NUM_THREADS,rf_threads=RF_NUM_THREADS,R_version=R.version.string,mgcv_version=as.character(packageVersion('mgcv')),ranger_version=as.character(packageVersion('ranger')),data_table_version=as.character(packageVersion('data.table'))))
  for(nm in names(tables)) fwrite(tables[[nm]],file.path(out,paste0(nm,'.csv')))
  cols <- c(LM='#377EB8',GAM='#B65C00',RF='#16836A')
  plot_curve <- function(work) {
    z <- summary[workload==work]; ylim <- range(z$log_mae_mean,z$log_mae_mean+ifelse(is.na(z$log_mae_sd),0,z$log_mae_sd),z$log_mae_mean-ifelse(is.na(z$log_mae_sd),0,z$log_mae_sd));ylim[1]<-max(0,ylim[1])
    png(file.path(out,paste0('learning_curve_',work,'.png')),width=1400,height=900,res=140)
    on.exit(dev.off(),add=TRUE);par(mar=c(5,5,4,2))
    plot(NA,xlim=range(z$training_size),ylim=ylim,log='x',xaxt='n',xlab='Eligible training queries (log scale)',ylab='Mean log MAE',main=gsub('_',' ',work))
    axis(1,at=sort(unique(z$training_size)),labels=format(sort(unique(z$training_size)),big.mark=',',scientific=FALSE,trim=TRUE),las=2,cex.axis=.75)
    for(id in models) {
      q <- z[model==id];setorder(q,training_size);co<-cols[sub('-[BA]$','',id)];ls<-if(endsWith(id,'-B'))2 else 1
      lines(q$training_size,q$log_mae_mean,col=co,lty=ls,lwd=2);points(q$training_size,q$log_mae_mean,pch=ifelse(q$origin=='9A_full_reference',18,16),col=co)
      q <- q[!is.na(log_mae_sd) & log_mae_sd>0];if(nrow(q)) arrows(q$training_size,pmax(0,q$log_mae_mean-q$log_mae_sd),q$training_size,q$log_mae_mean+q$log_mae_sd,angle=90,code=3,length=.035,col=co)
    }
    legend('topright',legend=models,col=cols[sub('-[BA]$','',models)],lty=ifelse(endsWith(models,'-B'),2,1),lwd=2,bty='n',cex=.85)
    mtext('Bars: +/- 1 SD across nested samples; diamond: single 9A reference (no SD).',side=3,line=.1,cex=.75)
  }
  for(work in names(evals)) plot_curve(work)
  writeLines(c('Fixed full-training-selected hyperparameters; no retuning at smaller sizes.',
    'One randomly selected anchor per training block, followed by a random permutation of remaining queries; rare blocks are overrepresented at small sizes.',
    'Each repetition shares its nested subsets across model configurations. RF fitting seed is fixed by the 9A helper; variability reflects training subsampling.',
    'Evaluation workloads are fixed; no holdout or CTA tuning. Previously unused holdout is a subset, not independent replication.',
    'The full endpoint reuses the single validated 9A fit per configuration, produced on the teammate\'s original R environment. It is not a repeated estimate; SD is NA.',
    'Error bars are descriptive sample-to-sample SD, not confidence intervals. GAM effective basis size may be limited by available distinct predictor values.',
    'Training blocks are defined by query centers; overlapping query footprints may remain dependent. No causal or significance claims.'),file.path(out,'methodology.txt'))
  cat(sprintf('\n9C CHECKS PASSED | new fits: %d | resumed: %d | reused full-size models: %d\n',fitted_now,resumed,if(reuse_full)6L else 0L))
  cat('Outputs:',out,'\n')
  invisible(list(new_fits=fitted_now,resumed=resumed,summary=summary))
}
main_9c()
