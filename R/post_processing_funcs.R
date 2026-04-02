themes <- ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 12, angle = 45, vjust = 0.5),
    axis.text.y = ggplot2::element_text(size = 16, face = "bold", angle = 0),
    axis.title.x = ggplot2::element_text(size = 20, face = "bold", hjust = 1),
    axis.title.y = ggplot2::element_text(size = 20, face = "bold", hjust = 1),
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 18),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 20, face = "bold"),
    legend.text = ggplot2::element_text(colour = "black", size = 20),
    strip.text.y = ggplot2::element_text(size = 20, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 20, face = "bold"),
    strip.background = ggplot2::element_rect(colour = "black", fill = "white"),
    legend.background = ggplot2::element_rect(
      fill = "white", linewidth = 0.5, linetype = "solid", colour = "black"
    )
  )

### ECLIPSE + LATITUDE LONGITUDE ALTITUDE ON ISS###
process_llae <- function(lla_pathname, match_date1, match_date2){
  lla_headers <- c("year", "month", "day", "julian_day", "hour", "minute",
                   "second", "lat", "long", "alt")
  
  lla_files <- base::list.files(lla_pathname, pattern = "SST_.*\\.txt", full.names = TRUE)
  lla_data <- list()
  for (file in lla_files) {
    file_name <- tools::file_path_sans_ext(basename(file))
    lla_data[[file_name]] <- data.table::fread(file, header = FALSE, col.names = lla_headers)
  }
  lla <- dplyr::bind_rows(lla_data, .id = "source_file") %>%
    dplyr::mutate(
      datetime = lubridate::make_datetime(year, month, day, hour, minute, second),
      utc_time_lla = as.POSIXct(datetime, tz = "UTC")
    ) %>%
    dplyr::filter(dplyr::between(utc_time_lla, match_date1, match_date2)) %>%
    dplyr::select(lat, long, alt, utc_time_lla, julian_day)
  return(lla)
}

### MAX CURRENT + ENERGY JOIN ###
process_max <- function(current, energy){
  current %>%
    dplyr::rename_with(~ stringr::str_replace(., "^V", "bin"), dplyr::starts_with("V")) %>%
    dplyr::rename(utc_time_epee = timestamp) %>%
    dplyr::select(utc_time_epee, dplyr::all_of(bin_select)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(bin_select), names_to = "energy_bin_org", values_to = "current_value_org") %>%
    dplyr::full_join(
      energy %>% dplyr::rename_with(~ stringr::str_replace(., "^V", "bin"), dplyr::starts_with("V")) %>%
        dplyr::rename(utc_time_epee = timestamp) %>%
        dplyr::select(utc_time_epee, dplyr::all_of(bin_select)) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(bin_select), names_to = "energy_bin_org", values_to = "energy_value_org"),
      dplyr::join_by(utc_time_epee, energy_bin_org)
    ) %>%
    dplyr::slice_max(order_by = current_value_org, by = utc_time_epee, with_ties = FALSE) %>%
    tidyr::drop_na(utc_time_epee) %>%
    dplyr::mutate(
      noise_floor_org = ifelse(current_value_org < 0.15, "noise", "valid"),
      scp_epee_org = (energy_value_org - ((0.5 * oxygen_mass * velocity^2) / charge))
    )
}

### SV ###
process_sv <- function(calendar_dates, sv_pathname, match_date1, match_date2, oxygen_mass, velocity, charge){
  sv_list <- foreach::foreach(day = seq_along(calendar_dates), .errorhandling = 'pass') %do% {
    day_ymd <- calendar_dates[day]
    file_path <- paste0(sv_pathname, "/", day_ymd, "_df.RData")
    print(file_path)
    loaded_name <- base::load(file_path)
    print(paste("Processing:", file_path, "→", loaded_name))
    df <- base::get(loaded_name)
    return(df)
  }
  str(sv_list)
  # Combine all loaded data frames
  sv <- dplyr::bind_rows(sv_list) %>%
    dplyr::rename(
      current_value_sv = current,
      energy_bin_sv = energy_bin,
      energy_value_sv = energy
    ) %>%
    dplyr::mutate(
      utc_time_sv = lubridate::ymd_hms(timestamp),
      scp_epee_sv = (energy_value_sv - ((0.5 * oxygen_mass * velocity^2) / charge))
    ) %>%
    dplyr::filter(dplyr::between(utc_time_sv, match_date1, match_date2)) %>%
    dplyr::select(-timestamp)
}

### FPMU ###
process_fpmu <- function(fpmu_pathname, match_date1, match_date2){
  newcolnames <- c("year", "month", "day", "hour", "min", "sec", "msec",
                   "lat", "long", "alt", "ne_wlp", "te_wlp", "ne_pip", "ne_gis",
                   "vf_wlp", "vf_fpp")
  files <- base::list.files(path = fpmu_pathname, recursive = TRUE, pattern = "*.txt", full.names = TRUE)
  df_list <- list()
  for (f in seq_along(files)) {
    df <- utils::read.table(files[[f]], stringsAsFactors = FALSE)
    file_id <- as.character(f)
    df <- stats::setNames(df, newcolnames)
    df$id <- file_id
    df_list[[file_id]] <- df
  }
  
  fpmu <- dplyr::bind_rows(df_list, .id = "source_file") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ replace(.x, .x == -999, NA))) %>%
    dplyr::distinct(year, month, day, hour, min, sec, msec, .keep_all = TRUE) %>%
    tidyr::unite("ymd", year:day, sep = "-") %>%
    tidyr::unite("hms", hour:sec, sep = ":") %>%
    tidyr::unite("time_full", ymd, hms, sep = " ") %>%
    dplyr::mutate(time_fpmu = lubridate::ymd_hms(time_full)) %>%
    dplyr::select(c(time_fpmu, lat, long, alt, ne_wlp, ne_pip, vf_wlp, vf_fpp)) %>%
    dplyr::filter(dplyr::between(time_fpmu, match_date1, match_date2)) %>%
    dplyr::rename(
      alt_fpmu = alt,
      lat_fpmu = lat,
      lon_fpmu = long
    )
}

### FULL JOIN ###
merge_data <- function(llae, epee_max, sv, fpmu, tol_sec = 1, prefer = "later") {
  prefer <- base::match.arg(prefer)
  
  full <- dplyr::full_join(epee_max, sv, dplyr::join_by(utc_time_epee == utc_time_sv), keep = TRUE) %>%
    dplyr::full_join(fpmu, dplyr::join_by(utc_time_epee == time_fpmu), keep = TRUE) %>%
    tidyr::drop_na(energy_bin_org) %>%
    dplyr::mutate(.rid = dplyr::row_number())
  
  # 1) First: exact timestamp match
  full_exact <- full %>%
    dplyr::left_join(
      (llae %>% dplyr::select(utc_time_lla, lat, long, alt)),
      by = c("utc_time_epee" = "utc_time_lla")
    )
  
  # 2) Rows that still need coords (all three missing)
  need <- full_exact %>%
    dplyr::filter(dplyr::if_all(c(lat, long, alt), is.na)) %>%
    dplyr::select(.rid, utc_time_epee)
  
  if (nrow(need) == 0L) {
    return(full_exact %>% dplyr::select(-.rid))  # everything matched exactly
  }
  
  # 3) Within ± tol_sec window, generate candidates and compute closeness
  cand <- fuzzyjoin::difference_inner_join(
    need %>% dplyr::transmute(.rid, epee_sec = as.numeric(utc_time_epee)),
    (llae %>% dplyr::select(utc_time_lla, lat, long, alt)) %>%
      dplyr::transmute(lla_sec = as.numeric(utc_time_lla), lat, long, alt),
    by = c("epee_sec" = "lla_sec"),
    max_dist = tol_sec,
    distance_col = ".abs_diff"
  ) %>%
    dplyr::mutate(.signed = epee_sec - lla_sec)  # negative => llae is later/equal than epee
  
  # 4) Pick the single nearest candidate per row, with a clear tie-break rule
  picked <- cand %>%
    dplyr::arrange(.rid, abs(.signed), .signed > 0) %>%
    dplyr::group_by(.rid) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(.rid, lat, long, alt)
  
  # 5) Patch the missing coords and clean up
  full_exact %>%
    dplyr::left_join(picked, by = ".rid", suffix = c("", "_fill")) %>%
    dplyr::mutate(
      lat  = dplyr::coalesce(lat,  lat_fill),
      long = dplyr::coalesce(long, long_fill),
      alt  = dplyr::coalesce(alt,  alt_fill)
      ) %>%
    dplyr::select(-.rid, -tidyselect::ends_with("_fill"))
}

### CURRENT AND ENERGY###
rawepee_process <- function(calendar_dates, rawdata_pathname, julian_days, file_level = "L2"){
  datalist <- foreach::foreach(day = seq_along(calendar_dates), .combine = 'list', .errorhandling = 'pass') %do% {
    day_ymd <- calendar_dates[day]
    day_julian <- julian_days[day]
    
    nc_file <- paste0(rawdata_pathname, '/STP-H9_EPEE_', day_ymd, '_', day_julian, '_L2.nc')
    print(nc_file)
    ncin <- ncdf4::nc_open(nc_file)
    if (file_level == "L1") {
      raw_time <- ncdf4::ncvar_get(ncin, "TYPE1_PKT_DATA_TIME_ARRAY")
      energy   <- ncdf4::ncvar_get(ncin, "TYPE1_PKT_SWEEP_VOLTAGE_ARRAY")
      current  <- ncdf4::ncvar_get(ncin, "TYPE1_PKT_SWEEP_DATA_ARRAY")
    } else if (file_level == "L2") {
      raw_time <- ncdf4::ncvar_get(ncin, "EPEETime")
      energy   <- ncdf4::ncvar_get(ncin, "EPEEEnergy")
      current  <- ncdf4::ncvar_get(ncin, "EPEECurrent")
    } else {
      stop("Unknown file_level: must be 'L1' or 'L2'")
    }
    ncdf4::nc_close(ncin)
    
    epee_epoch <- as.POSIXct("1980-01-06 00:00:00", tz = "UTC")
    time_vec <- epee_epoch + raw_time
    
    energy_df <- as.data.frame(energy) %>%
      dplyr::mutate(timestamp = time_vec,
                    timeind = base::seq_along(time_vec),
                    timemin = get_timesec(timestamp) / 60)
    
    current_df <- as.data.frame(current) %>%
      dplyr::mutate(timestamp = time_vec,
                    timeind = base::seq_along(time_vec),
                    timemin = get_timesec(timestamp) / 60)
    
    list(energy_df = energy_df, current_df = current_df)
  }
  
  energy_all  <- dplyr::bind_rows(lapply(datalist, `[[`, "energy_df"))
  current_all <- dplyr::bind_rows(lapply(datalist, `[[`, "current_df"))
  
  return(list(energy = energy_all, current = current_all))
}

estimate_noise_floor <- function(current_full, energy_full, calendar_dates, calendar_date_choice, quarter_day){
  df_curr_all <- current_full %>%
    dplyr::filter(lubridate::as_date(timestamp) == lubridate::as_date(calendar_dates[calendar_date_choice]))
  df_enrg_all <- energy_full %>%
    dplyr::filter(lubridate::as_date(timestamp) == lubridate::as_date(calendar_dates[calendar_date_choice]))
  mat_curr_all <- as.matrix(df_curr_all[, ind_bins])
  mat_enrg_all <- as.matrix(df_enrg_all[, ind_bins])
  
  n_time <- nrow(mat_curr_all)
  chunks <- base::split(1:n_time, base::cut(1:n_time, breaks = 4, labels = FALSE))
  
  time_idx <- chunks[[quarter_day]]
  df_curr <- df_curr_all[time_idx, ]
  df_enrg <- df_enrg_all[time_idx, ]
  mat_curr <- mat_curr_all[time_idx, ]
  mat_enrg <- mat_enrg_all[time_idx, ]
  
  id_vars <- c("timestamp","timeind","timemin")
  df_long <- reshape2::melt(df_curr, id.vars = id_vars)
  df_long$energy <- reshape2::melt(df_enrg, id.vars = id_vars)$value
  df_long$energy_bin <- as.numeric(df_long$variable)
  energy_bins <- levels(df_long$variable)
  
  df_max <- get_max_df(df_curr, df_enrg, energy_bins)
  
  all_bg <- max(df_long$value) <= 0.18
  df_long$predmu <- df_long$value
  if (!all_bg) {
    svec_curr_inputs <- as.matrix(df_long[, c('timemin','energy')])
    fit_curr <- fit_scaled(y = df_long$value, inputs = svec_curr_inputs, nug = NULL, n.est = 5e4)
    df_long$predmu <- predictions_scaled(fit_curr, svec_curr_inputs, m = 100, joint = FALSE, predvar = FALSE)
  }
  
  trust_enrg <- 2:97
  mat_smoothcurr <- unmelt(df_long, id_vars, value.var = 'predmu')
  
  bg_timeinds <- get_noisefloor_timeinds(df_long, df_max, pcnt_bad = 1, bad_bins = 50:100)
  if (length(bg_timeinds) > 0) {
    bg_matinds <- sapply(bg_timeinds, function(x) which(x == df_max$timeind))
  } else {
    sum_abs_diffs <- apply(mat_smoothcurr[, 1:30], 1, function(row) {
      diffs <- diff(row)
      sum(abs(diffs))
    })
    srtd_chngs <- sort(sum_abs_diffs, index.return = TRUE)$ix
    bg_matinds <- srtd_chngs[1:round(0.001 * length(srtd_chngs))]
  }
  
  floor_profile <- apply(mat_curr[bg_matinds, 1:100], 2, mean)
  floor_profile_all <- matrix(NA, nrow = length(bg_matinds), ncol = length(floor_profile))
  
  for (bb in 1:length(bg_matinds)) {
    tmp_bg <- unlist(mat_smoothcurr[bg_matinds[bb], 1:100])
    tmp <- get_optim_noisefloor(
      data_profile = tmp_bg,
      richards_energy_bins = 1:97,
      gaussian_energy_bins = 2:20,
      return_gaussian = TRUE,
      parabolic_adjustment = TRUE,
      verbose = FALSE
    )
    floor_profile_all[bb, ] <- tmp$floor_profile
    floor_profile_all[bb, trust_enrg] <- tmp$floor_profile[trust_enrg]
  }
  
  bb_kp <- which.min(apply(floor_profile_all[, 2:20], 1, sum))
  floor_profile[trust_enrg] <- floor_profile_all[bb_kp, trust_enrg]
  return(list(
    df_long = df_long,
    noisefloor_timeinds = df_max$timeind[bg_matinds],
    floor_profile = floor_profile,
    tmp = tmp,
    tmp_bg = tmp_bg
  ))
}
