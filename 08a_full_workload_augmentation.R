# CSP 571 | Full rectangular workload augmentation (Step 8A)
# Run from repository root. Default: benchmark only. No raw downloads.
# Full run: Sys.setenv(DPA08_MODE='full'); source('08a_full_workload_augmentation.R')
# Optional: DPA08_BATCH_SIZE=1000 (default), DPA08_BENCH_N=1000.
# Uses same projected polygon definition as Steps 4B and 7A.
# Checkpoints are written under data/raw (ignored by Git).

run08 <- function() {
  if (!requireNamespace('sf', quietly=TRUE)) stop('Install sf: install.packages("sf")')
  library(sf)
  old_s2 <- sf_use_s2(); on.exit(sf_use_s2(old_s2), add=TRUE)
  sf_use_s2(FALSE)
  mode <- Sys.getenv('DPA08_MODE', 'benchmark')
  batch_size <- as.integer(Sys.getenv('DPA08_BATCH_SIZE', '1000'))
  bench_n <- as.integer(Sys.getenv('DPA08_BENCH_N', '1000'))
  check <- function(x, msg) if (!isTRUE(x)) stop(msg, call.=FALSE)
  check(mode %in% c('benchmark','full'), 'DPA08_MODE must be benchmark or full.')
  check(is.finite(batch_size) && batch_size >= 100L && batch_size <= 5000L, 'Batch size must be 100–5000.')
  check(is.finite(bench_n) && bench_n >= 100L && bench_n <= 5000L, 'Benchmark size must be 100–5000.')
  out <- 'outputs'; proc <- 'data/processed'; cache <- 'data/raw/08a_checkpoints'
  for (p in c(out, proc, cache)) dir.create(p, recursive=TRUE, showWarnings=FALSE)
  uci_files <- list.files('data/raw/uci', pattern='^Range-Queries-Aggregates\\.csv$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
  check(length(uci_files)==1L, 'Expected one raw Range-Queries-Aggregates.csv. Run Step 1 if missing.')
  paths <- c(uci=uci_files[1], crime=file.path(proc,'chicago_crimes_2024_clean.csv'),
             station=file.path(proc,'station_crime_ridership_2024.csv'),
             pilot=file.path(proc,'07a_coverage_aware_workloads_2024.csv'),
             split=file.path(out,'05_uci_spatial_split.csv'),
             boundary='data/Boundaries_-_City_20261121.geojson')
  check(all(file.exists(paths)), paste('Missing inputs:',paste(paths[!file.exists(paths)],collapse=', ')))
  read <- function(p) read.csv(p, stringsAsFactors=FALSE, check.names=FALSE)
  t0 <- proc.time()[['elapsed']]
  u <- read(paths['uci']); names(u) <- tolower(trimws(names(u)))
  req <- c('x','y','x_range','y_range','count','sum_','avg')
  check(all(req %in% names(u)), 'UCI schema differs from Step 4B.')
  for (v in req) u[[v]] <- suppressWarnings(as.numeric(u[[v]]))
  u$source_row_number <- seq_len(nrow(u))
  valid <- is.finite(u$x)&is.finite(u$y)&is.finite(u$x_range)&is.finite(u$y_range)&u$x_range>0&u$y_range>0
  check(nrow(u)==200000L, 'Expected 200,000 raw UCI records.')
  check(all(valid), paste('Invalid rectangle geometry on',sum(!valid),'rows. Audit before continuing; no silent dropping.'))
  u$query_id <- paste0('UCI_ROW_',u$source_row_number)
  pilot <- read(paths['pilot']); sp <- read(paths['split'])
  pu <- pilot[pilot$query_source=='UCI_original_rectangle',,drop=FALSE]
  cta <- pilot[pilot$query_source=='CTA_centered_rectangle',,drop=FALSE]
  check(nrow(pu)==800L && nrow(cta)==426L && !anyDuplicated(pilot$query_id), 'Unexpected pilot workload.')
  check(!anyDuplicated(sp$query_id) && setequal(sp$query_id,pu$query_id), 'Invalid original split IDs.')
  held <- unique(sp$spatial_block[sp$split=='UCI_spatial_holdout'])
  trained <- unique(sp$spatial_block[sp$split=='UCI_train'])
  check(!length(intersect(held,trained)), 'Original split has shared blocks.')
  crime <- read(paths['crime']); stations <- read(paths['station'])
  stations <- stations[stations$query_radius_m==500,,drop=FALSE]
  check(nrow(stations)==142L && !anyDuplicated(stations$station_id), 'Expected 142 unique active stations.')
  check(nrow(crime)==256583L && all(is.finite(crime$latitude)&is.finite(crime$longitude)), 'Crime inputs differ from pilot.')
  check(all(is.finite(stations$latitude)&is.finite(stations$longitude)&is.finite(stations$total_entries_2024)&stations$total_entries_2024>=0), 'Invalid station features.')
  project_points <- function(x,lon,lat) st_transform(st_geometry(st_as_sf(x,coords=c(lon,lat),crs=4326,remove=FALSE)),26916)
  crime_g <- project_points(crime,'longitude','latitude'); station_g <- project_points(stations,'longitude','latitude')
  types <- toupper(trimws(as.character(crime$primary_type)))
  entries <- as.numeric(stations$total_entries_2024)
  city <- st_read(paths['boundary'],quiet=TRUE)
  city <- st_make_valid(st_transform(st_union(st_make_valid(city)),26916))
  check(length(city)==1L && all(st_is_valid(city)) && !any(st_is_empty(city)), 'Invalid city boundary.')
  classify <- function(p) ifelse(p>=99.99,'fully covered',ifelse(p<=0.0001,'outside Chicago','partially covered'))
  rect <- function(x,y,dx,dy) st_polygon(list(matrix(c(x-dx,y-dy,x+dx,y-dy,x+dx,y+dy,x-dx,y+dy,x-dx,y-dy),ncol=2,byrow=TRUE)))
  ft <- 1200/3937
  extra_features <- function(d, width, height, center_g) {
    d$width_m <- width; d$height_m <- height
    d$aspect_ratio <- width/height
    d$perimeter_m <- 2*(width+height)
    nearest <- st_nearest_feature(center_g,station_g)
    d$nearest_station_distance_m <- as.numeric(st_distance(center_g,station_g[nearest],by_element=TRUE))
    d$station_density_per_km2 <- d$stations_inside_query/d$area_km2
    d$ridership_density_per_km2 <- d$total_cta_entries_at_stations_inside_query_2024/d$area_km2
    d
  }
  # Batched GEOS: only this batch's point-membership lists exist at any time.
  # Do not substitute a bounding box in UTM: projected UCI rectangles may rotate.
  process_batch <- function(idx) {
    q <- u[idx,,drop=FALSE]
    centers <- st_geometry(st_as_sf(q,coords=c('x','y'),crs=3435,remove=FALSE))
    ll <- st_coordinates(st_transform(centers,4326))
    check(all(is.finite(ll)) && all(ll[,1]> -88.2 & ll[,1]< -87.2 & ll[,2]>41.5 & ll[,2]<42.2), 'Unexpected UCI center coordinates; inspect inferred CRS.')
    geom <- st_transform(st_sfc(lapply(seq_len(nrow(q)),function(i) rect(q$x[i],q$y[i],q$x_range[i],q$y_range[i])),crs=3435),26916)
    area <- as.numeric(st_area(geom)); check(all(is.finite(area)&area>0),'Invalid query areas.')
    tm <- proc.time()[['elapsed']]
    inc <- st_intersects(geom,crime_g,sparse=TRUE)
    cnt <- lengths(inc)
    count_type <- function(k) vapply(inc,function(j) sum(types[j]==k,na.rm=TRUE),integer(1))
    theft <- count_type('THEFT'); robbery <- count_type('ROBBERY'); battery <- count_type('BATTERY'); assault <- count_type('ASSAULT')
    hit_n <- sum(as.double(cnt)); rm(inc)
    crime_sec <- proc.time()[['elapsed']]-tm
    tm <- proc.time()[['elapsed']]
    stn <- st_intersects(geom,station_g,sparse=TRUE)
    stn_n <- lengths(stn); stn_ent <- vapply(stn,function(j) sum(entries[j]),numeric(1)); rm(stn)
    # Exact coverage for boundary-straddling queries; fully contained polygons
    # have their entire area inside. Disjoint polygons have zero inside area.
    full <- lengths(st_covered_by(geom,city,sparse=TRUE))>0L
    touches <- lengths(st_intersects(geom,city,sparse=TRUE))>0L
    inside <- numeric(length(idx)); inside[full] <- area[full]
    partial <- which(touches & !full)
    if (length(partial)) {
      inter <- suppressWarnings(st_intersection(st_sf(i=partial,geometry=geom[partial]),city))
      if (nrow(inter)) {
        sums <- tapply(as.numeric(st_area(inter)),inter$i,sum)
        inside[as.integer(names(sums))] <- as.numeric(sums)
      }
    }
    check(all(is.finite(inside)&inside>=0&inside<=area+pmax(1e-4,1e-8*area)), 'Invalid boundary overlap.')
    pct <- pmin(100,pmax(0,100*inside/area)); coverage <- classify(pct)
    block <- paste(floor((ll[,2]-41)/0.020),floor((ll[,1]+88)/0.025),sep='_')
    d <- data.frame(query_id=q$query_id,query_source='UCI_original_rectangle',geometry_type='rectangle',
      query_halfwidth_m=NA_real_,station_id=NA_character_,station_name=NA_character_,area_km2=area/1e6,
      original_uci_count=q$count,original_uci_sum_arrests=q$sum_,original_uci_avg_beat=q$avg,
      crimes_2024=cnt,theft_2024=theft,robbery_2024=robbery,battery_2024=battery,assault_2024=assault,
      stations_inside_query=stn_n,total_cta_entries_at_stations_inside_query_2024=stn_ent,
      center_lon=ll[,1],center_lat=ll[,2],crimes_per_km2_2024=cnt/(area/1e6),
      split=ifelse(block %in% held,'UCI_spatial_holdout','UCI_train'),spatial_block=block,
      pct_inside_chicago=pct,coverage_status=coverage,inside_city_area_km2=inside/1e6,
      zero_crime_query=cnt==0,eligible_full_coverage=coverage=='fully covered',
      stringsAsFactors=FALSE)
    d <- extra_features(d,2*q$x_range*ft,2*q$y_range*ft,st_transform(centers,26916))
    d$is_poc_query <- d$query_id %in% pu$query_id
    d$source_row_number <- q$source_row_number
    list(data=d,crime_seconds=crime_sec,other_seconds=proc.time()[['elapsed']]-tm,point_memberships=hit_n)
  }
  # Validate all 800 original UCI rows against the committed pilot output.
  pilot_idx <- match(pu$query_id,u$query_id)
  check(!anyNA(pilot_idx),'Pilot IDs absent from original file.')
  cat('Checking 800 pilot queries against Step 7A...\n')
  pilot_batches <- split(pilot_idx,ceiling(seq_along(pilot_idx)/batch_size))
  pilot_results <- lapply(pilot_batches,process_batch)
  actual <- do.call(rbind,lapply(pilot_results,`[[`,'data'))
  actual <- actual[match(pu$query_id,actual$query_id),,drop=FALSE]
  exact_cols <- c('crimes_2024','theft_2024','robbery_2024','battery_2024','assault_2024','stations_inside_query','total_cta_entries_at_stations_inside_query_2024','split','spatial_block','coverage_status','zero_crime_query','eligible_full_coverage')
  validation <- do.call(rbind,lapply(exact_cols,function(k) data.frame(field=k,mismatches=sum(is.na(actual[[k]])!=is.na(pu[[k]]) | (!is.na(actual[[k]]) & !is.na(pu[[k]]) & actual[[k]]!=pu[[k]])))))
  for (k in c('area_km2','pct_inside_chicago','center_lon','center_lat')) {
    tol <- if(k=='pct_inside_chicago') 1e-5 else if(k=='area_km2') 1e-8+pmax(abs(pu[[k]]),1)*1e-6 else 1e-8
    validation <- rbind(validation,data.frame(field=k,mismatches=sum(abs(actual[[k]]-pu[[k]])>tol)))
  }
  write.csv(validation,file.path(out,'08a_pilot_validation.csv'),row.names=FALSE)
  print(validation,row.names=FALSE)
  check(all(validation$mismatches==0),'Pilot mismatch. Stop and inspect outputs/08a_pilot_validation.csv before full processing.')
  # Independent random sample for realistic processing speed (not a model split).
  set.seed(5712026)
  bench_idx <- sample(setdiff(seq_len(nrow(u)),pilot_idx),bench_n)
  cat('Benchmarking',bench_n,'additional UCI queries...\n')
  b0 <- proc.time()[['elapsed']]
  bb <- lapply(split(bench_idx,ceiling(seq_along(bench_idx)/batch_size)),process_batch)
  seconds <- proc.time()[['elapsed']]-b0
  bench <- data.frame(n_queries=bench_n,batch_size=batch_size,elapsed_seconds=seconds,
    crime_seconds=sum(vapply(bb,`[[`,numeric(1),'crime_seconds')),
    other_seconds=sum(vapply(bb,`[[`,numeric(1),'other_seconds')),
    point_memberships=sum(vapply(bb,`[[`,numeric(1),'point_memberships')),
    rough_full_seconds=seconds*nrow(u)/bench_n)
  write.csv(bench,file.path(out,'08a_benchmark.csv'),row.names=FALSE)
  print(bench,row.names=FALSE)
  cat('Pilot validation PASSED. Rough full-run estimate (not a guarantee):',round(bench$rough_full_seconds/60,1),'minutes.\n')
  if(mode=='benchmark') {cat('BENCHMARK COMPLETE. No full-scale dataset written.\n'); return(invisible(bench))}
  rm(bb,pilot_results,actual); gc()
  # Only reuse checkpoints from the same input bytes, batch size and code version.
  signature <- list(version='08a-v1',input_md5=unname(tools::md5sum(paths)),batch_size=batch_size,sf=as.character(packageVersion('sf')),external=sf_extSoftVersion())
  manifest <- file.path(cache,'manifest.rds')
  if(file.exists(manifest)) check(identical(readRDS(manifest),signature),'Checkpoint inputs/settings changed. Rename data/raw/08a_checkpoints before restarting.') else saveRDS(signature,manifest)
  groups <- split(seq_len(nrow(u)),ceiling(seq_len(nrow(u))/batch_size))
  files <- file.path(cache,sprintf('batch_%05d.rds',seq_along(groups)))
  for(i in seq_along(groups)) {
    if(file.exists(files[i])) {
      z <- readRDS(files[i]); check(identical(z$data$query_id,u$query_id[groups[[i]]]),'Checkpoint query IDs mismatch.')
    } else {
      z <- process_batch(groups[[i]])
      tmp <- paste0(files[i],'.tmp'); saveRDS(z,tmp,compress=FALSE)
      check(file.rename(tmp,files[i]),'Could not finalize checkpoint.')
    }
    cat(sprintf('Batch %d/%d | %d queries | %.1fs crime + %.1fs other\n',i,length(groups),nrow(z$data),z$crime_seconds,z$other_seconds))
    rm(z); if(i%%10L==0L) gc()
  }
  # Assemble only aggregated rows, never the individual incident join.
  all_u <- do.call(rbind,lapply(files,function(p) readRDS(p)$data))
  cta_centers <- project_points(cta,'center_lon','center_lat')
  cta <- extra_features(cta,2*cta$query_halfwidth_m,2*cta$query_halfwidth_m,cta_centers)
  cta$is_poc_query <- TRUE; cta$source_row_number <- NA_integer_
  check(setequal(names(all_u),names(cta)),'CTA/UCI output schema mismatch.')
  all <- rbind(all_u,cta[,names(all_u),drop=FALSE])
  check(nrow(all)==200426L && !anyNA(all$query_id) && !anyDuplicated(all$query_id),'Final workload ID/count mismatch.')
  check(setequal(all_u$query_id,u$query_id),'Full workload missed UCI rows.')
  check(!length(intersect(all_u$spatial_block[all_u$split=='UCI_train'],all_u$spatial_block[all_u$split=='UCI_spatial_holdout'])),'Training/holdout block leakage.')
  check(all(is.finite(all$area_km2)&all$area_km2>0&is.finite(all$crimes_2024)&all$crimes_2024>=0),'Invalid primary features/target.')
  check(all(all$eligible_full_coverage==(all$coverage_status=='fully covered')),'Coverage eligibility mismatch.')
  check(all(all$zero_crime_query==(all$crimes_2024==0)),'Zero flag mismatch.')
  # Original test blocks are preserved. New blocks go to training; new rows
  # within old holdout blocks are separately identifiable by is_poc_query.
  summary <- aggregate(rep(1L,nrow(all)),list(split=all$split,coverage_status=all$coverage_status),sum); names(summary)[3] <- 'n_queries'
  write.csv(summary,file.path(out,'08a_full_coverage_summary.csv'),row.names=FALSE)
  dest <- file.path(proc,'08_full_workload_2024.csv'); tmp <- paste0(dest,'.tmp')
  write.csv(all,tmp,row.names=FALSE); if(file.exists(dest)) unlink(dest)
  check(file.rename(tmp,dest),'Could not finalize full dataset.')
  lines <- c('STEP 8A: FULL RECTANGULAR WORKLOAD AUGMENTATION',
    paste('Total rows:',nrow(all)),paste('Original UCI rows:',nrow(all_u)),paste('CTA rows:',nrow(cta)),
    paste('Fully covered UCI training rows:',sum(all$split=='UCI_train' & all$eligible_full_coverage)),
    paste('Fully covered UCI holdout rows:',sum(all$split=='UCI_spatial_holdout' & all$eligible_full_coverage)),
    paste('Previously unused fully covered UCI holdout rows:',sum(all$split=='UCI_spatial_holdout' & all$eligible_full_coverage & !all$is_poc_query)),
    paste('Fully covered CTA rows:',sum(all$split=='CTA_transfer' & all$eligible_full_coverage)),
    paste('Training blocks:',length(unique(all$spatial_block[all$split=='UCI_train' & all$eligible_full_coverage]))),
    paste('Elapsed seconds this invocation:',round(proc.time()[['elapsed']]-t0,1)),
    'All 800 pilot rows passed exact count/coverage-class/split checks and numeric tolerance checks.',
    'Original UCI CRS inferred EPSG:3435; same polygon reconstruction as pilot, UTM EPSG:26916 for overlay.',
    'New center blocks assigned to training; original held-out blocks never enter training. Overlapping footprints may share incidents.',
    'Width/height are native projected dimensions converted to meters; nearest station distance is center-to-station in UTM.',
    'Coverage is a jurisdiction proxy, not completeness; retain all rows including zeros.',
    'Historical UCI aggregates are provenance only; never use them as primary 2024 model features.',
    '',capture.output(print(summary,row.names=FALSE)), '',capture.output(sessionInfo()))
  writeLines(lines,file.path(out,'08_processing_summary.txt'))
  cat(paste(lines[1:10],collapse='\n'),'\nFULL WORKLOAD COMPLETE\n')
  invisible(all)
}
run08()
