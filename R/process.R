### VARIABLES ####
all_dates      <- seq.Date(from = as.Date("2023-09-28"), to = as.Date("2023-09-29"), by = "day")
calendar_dates <- format(all_dates, "%Y%m%d")
julian_days    <- format(all_dates, "%j")
ind_bins       <- 1:100
enrg_cut       <- 40
alpha          <- 0.05

# output dirs
process_fig_dir <- here::here("figures", "process")
log_dir         <- here::here("logfiles")

# Detect available cores and leave some free for other processes
num_cores <- max(1, parallel::detectCores() - 2)
cl <- parallel::makeCluster(num_cores)  # PSOCK on mac

# Export the project root so workers can source R/ files without file_path
project_root <- normalizePath(here::here(), mustWork = TRUE)

parallel::clusterExport(
  cl,
  varlist = c(
    "project_root",
    # data/paths
    "sv_pathname", "rawdata_pathname", "process_fig_dir", "log_dir",
    # config / indices
    "calendar_dates", "julian_days", "ind_bins", "enrg_cut"
  ),
  envir = environment()
)

# Source functions and activate renv ONCE per worker (no library() calls needed if you use pkg::)
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(project_root, "renv", "activate.R"))
  source(file.path(project_root, "R", "process_funcs.R"))
  source(file.path(project_root, "R", "post_processing_funcs.R"))
}))

doParallel::registerDoParallel(cl)

# write print messages and errors to a logfile
foreach::foreach(
  day = seq_along(calendar_dates),
  .combine = 'c',
  .errorhandling = 'pass'
) %dopar% {
  tryCatch({
    day_ymd    <- calendar_dates[day]
    day_julian <- julian_days[day]
    
    # log
    log_file <- file.path(log_dir, sprintf("logs_%s_%s.txt", day_ymd, day_julian))
    log_conn <- file(log_file, open = "w")
    on.exit(try(close(log_conn), silent = TRUE), add = TRUE)
    writeLines(paste(Sys.time(), "Starting", day_ymd), log_conn)
    
    # bring in netcdf L2 file and process
    nc_file   <- file.path(rawdata_pathname, sprintf("STP-H9_EPEE_%s_%s_L2.nc", day_ymd, day_julian))
    data_list <- load_netcdf_data(nc_file, log_conn)
    if (is.null(data_list)) {
      writeLines(paste(Sys.time(), "Skipping", day_ymd, "- load_netcdf_data returned NULL"), log_conn)
      return(NA_character_)
    }
    
    df_enrg_all  <- data_list$energy_df
    df_curr_all  <- data_list$current_df
    mat_curr_all <- as.matrix(df_curr_all[, ind_bins])
    mat_enrg_all <- as.matrix(df_enrg_all[, ind_bins])
    
    # split day into 4 chunks
    df_out <- df_long_all <- df_max_all <- NULL
    n_time <- nrow(mat_curr_all)
    chunks <- split(seq_len(n_time), cut(seq_len(n_time), breaks = 4, labels = FALSE))
    
    out_list  <- vector("list", length(chunks))
    max_list  <- vector("list", length(chunks))
    long_list <- vector("list", length(chunks))
    
    for (kk in seq_along(chunks)) {
      time_idx <- chunks[[kk]]
      df_curr  <- df_curr_all[time_idx, ]
      df_enrg  <- df_enrg_all[time_idx, ]
      mat_curr <- mat_curr_all[time_idx, ]
      mat_enrg <- mat_enrg_all[time_idx, ]
      writeLines(paste(Sys.time(), "Matrix conversion done for", day_ymd, "-chunk", kk), log_conn)
      
      # long form (use data.table::melt explicitly)
      id_vars  <- c("timestamp", "timeind", "timemin")
      df_long  <- reshape2::melt(df_curr, id.vars = id_vars)
      df_long$energy     <- reshape2::melt(df_enrg, id.vars = id_vars)$value
      df_long$energy_bin <- as.numeric(df_long$variable)
      energy_bins        <- levels(df_long$variable)
      
      # STEP 0a: max for each time step (raw)
      df_max <- get_max_df(df_curr, df_enrg, energy_bins)
      
      # STEP 1a: GP fit if not all background
      all_bg         <- max(df_long$value) <= 0.18
      df_long$predmu <- df_long$value
      if (!all_bg) {
        svec_curr_inputs <- as.matrix(df_long[, c("timemin", "energy")])
        
        #BLOCK FOR HOLD-OUT PREDICTIONS###
        #######
        set.seed(123)
        n_rows <- nrow(df_long)
        idx_te <- sample.int(n_rows, size = floor(0.10 * n_rows))
        idx_tr <- setdiff(seq_len(n_rows), idx_te)
        
        # fit on training only
        fit_tr <- fit_scaled(y = df_long$value[idx_tr], inputs = svec_curr_inputs[idx_tr, , drop = FALSE],
                             ms = c(30), nu = 1.5, nug = NULL, n.est = min(5000, length(idx_tr)), scale = "parms")
        
        # predict variances on holdout
        pred_hold <- predictions_scaled(fit = fit_tr, locs_pred = svec_curr_inputs[idx_te, , drop = FALSE],
                                        m = 100, joint = FALSE, predvar = TRUE, scale = "parms")
        sd_hold <- sqrt(pred_hold$vars)
        
        # robust signal amplitude (5th-95th percentile)
        amp_robust <- diff(quantile(df_long$value, probs = c(0.05, 0.95), na.rm = TRUE))
        median_sd_pct_chunk <- 100 * median(sd_hold, na.rm = TRUE) / amp_robust
        p90_sd_pct_chunk    <- 100 * quantile(sd_hold, 0.90, na.rm = TRUE) / amp_robust
        
        # write to log (log_conn exists in your scope)
        writeLines(sprintf("Chunk %d holdout median SD = %.4g%% ; p90 = %.4g%%", kk, median_sd_pct_chunk, p90_sd_pct_chunk),
                   log_conn)
        
        # full-chunk fit for production
        fit_curr <- fit_scaled(y = df_long$value, inputs = svec_curr_inputs, ms = c(30),
                               nu = 1.5, nug = NULL, n.est = min(5e4, nrow(svec_curr_inputs)), scale = "parms")
        #####
        
        
        
        ##fit_curr <- fit_scaled(y = df_long$value, inputs = svec_curr_inputs, nug = NULL, n.est = 5e4)
        writeLines(paste(Sys.time(), "Current fit for", day_ymd), log_conn)
        
        ind_bigcurr      <- (df_long$value > 0.15)
        timeinds_runall  <- setdiff(unique(df_long$timeind), unique(df_long$timeind[ind_bigcurr]))
        ind_timerun      <- df_long$timeind %in% timeinds_runall
        ind_lowbin       <- (df_long$energy <= (enrg_cut + 10))
        inds             <- (ind_lowbin | ind_timerun | ind_bigcurr)
        df_long$predmu[inds] <- predictions_scaled(fit_curr, svec_curr_inputs[inds, ], m = 100, joint = FALSE, predvar = FALSE)
      }
      writeLines(paste(Sys.time(), "SVecchia fit for", day_ymd), log_conn)
      
      # STEP 1b: noise floor estimate
      trust_enrg      <- 2:97
      mat_smoothcurr  <- unmelt(df_long, id_vars, value.var = "predmu")
      bg_timeinds     <- get_noisefloor_timeinds(df_long, df_max, pcnt_bad = 1, bad_bins = 50:100)
      if (length(bg_timeinds) > 0) {
        bg_matinds <- sapply(bg_timeinds, function(x) which(x == df_max$timeind))
      } else {
        sum_abs_diffs <- apply(mat_smoothcurr[, 1:30], 1, function(row) sum(abs(diff(row))))
        srtd_chngs    <- sort(sum_abs_diffs, index.return = TRUE)$ix
        bg_matinds    <- srtd_chngs[1:round(0.001 * length(srtd_chngs))]
      }
      floor_profile     <- apply(mat_curr[bg_matinds, 1:100], 2, mean)
      floor_profile_all <- matrix(NA, nrow = length(bg_matinds), ncol = length(floor_profile))
      writeLines(paste("Floor profile all rows:", ncol(floor_profile_all)), log_conn)
      
      for (bb in seq_along(bg_matinds)) {
        tmp_bg <- unlist(mat_smoothcurr[bg_matinds[bb], 1:100])
        tmp <- get_optim_noisefloor(
          data_profile = tmp_bg,
          richards_energy_bins = 1:97,
          gaussian_energy_bins = 2:20,
          parabolic_adjustment = TRUE, verbose = FALSE
        )
        floor_profile_all[bb, ]          <- tmp$floor_profile
        floor_profile_all[bb, trust_enrg] <- tmp$floor_profile[trust_enrg]
      }
      writeLines(paste("Floor profile rows from get_optim_noisefloor:", ncol(floor_profile)), log_conn)
      bb_kp <- which.min(apply(floor_profile_all[, 2:20, drop = FALSE], 1, sum))
      floor_profile[trust_enrg] <- floor_profile_all[bb_kp, trust_enrg]
      
      # diagnostic plot (explicit grDevices/graphics)
      grDevices::png(file.path(process_fig_dir, sprintf("%s_%d.png", day_ymd, kk)), width = 1000, height = 400)
      graphics::par(mfrow = c(1, 2))
      graphics::matplot(t(mat_curr[bg_matinds, ]), type = "l", xlab = "energy bin", ylab = "current")
      graphics::lines(floor_profile, lwd = 2, col = "red")
      graphics::matplot(t(mat_curr[bg_matinds, ]) - floor_profile, type = "l", xlab = "energy bin", ylab = "current - noise floor")
      graphics::abline(h = 0, col = "red", lwd = 2)
      graphics::par(mfrow = c(1, 1))
      grDevices::dev.off()
      writeLines(paste(Sys.time(), "Noise floor plot generated for", day_ymd), log_conn)
      
      # STEP 1c: remove noise floor
      df_floor <- df_curr
      df_floor[, ind_bins] <- matrix(rep(floor_profile, nrow(df_curr)), nrow = nrow(df_curr), ncol = length(floor_profile), byrow = TRUE)
      writeLines(paste("df_floor:", paste(dim(df_floor), collapse = "x")), log_conn)
      
      df_long$floor           <- reshape2::melt(df_floor, id.vars = id_vars)$value
      df_long$predmu_nofloor  <- df_long$predmu - df_long$floor
      
      df_smoothcurrent_nofloor <- unmelt(df_long, id_vars, value.var = "predmu_nofloor")
      df_max_nofloor           <- get_max_df(df_smoothcurrent_nofloor, df_enrg, energy_bins, bin_inds = 1:97)
      
      # STEP 3: build output slice
      df_tmp <- df_max_nofloor[, c("current", "energy_bin", "energy", "timeind", "timestamp")]
      df_tmp$current_omitted     <- NA_real_
      df_tmp$energy_omitted      <- NA_real_
      df_tmp$energy_bin_omitted  <- NA_real_
      
      out_list[[kk]]  <- df_tmp
      max_list[[kk]]  <- df_max
      long_list[[kk]] <- df_long
      
      cum_rows <- sum(vapply(out_list[seq_len(kk)], nrow, integer(1), USE.NAMES = FALSE))
      writeLines(paste("Combined df_out rows for", day_ymd, ":", cum_rows), log_conn)
    }
    
    # combine chunks
    df_out      <- dplyr::bind_rows(out_list)
    df_max_all  <- dplyr::bind_rows(max_list)
    df_long_all <- dplyr::bind_rows(long_list)
    writeLines(paste("Final combined df_out rows for", day_ymd, ":", nrow(df_out)), log_conn)
    
    # Postprocess
    df_tmp <- df_out
    writeLines(paste(Sys.time(), "Step 1.0: Checking df_max_all structure"), log_conn)
    writeLines(utils::capture.output(str(df_max_all)), log_conn)
    
    if (!"energy_bin" %in% names(df_max_all)) stop("Column 'energy_bin' missing from df_max_all")
    if (!is.numeric(df_max_all$energy_bin))   stop("Column 'energy_bin' is not numeric")
    
    writeLines(paste(Sys.time(), "Step 1.1: Passed checks, creating bad_enrg"), log_conn)
    df_max_all$bad_enrg <- df_max_all$energy_bin > 19
    
    writeLines("Step 2: Running rle and computing bad_run_span", log_conn)
    rle_bad <- base::rle(1 * df_max_all$bad_enrg)
    df_max_all$bad_run_span <- base::inverse.rle(list(
      lengths = rle_bad$lengths,
      values  = ifelse(rle_bad$values == 1, rle_bad$lengths, 0)
    ))
    
    writeLines("Step 3: Calculating bad_run_pcnt and removing spans", log_conn)
    df_max_all$bad_run_pcnt <- 100 * df_max_all$bad_run_span / nrow(df_max_all)
    indices_to_remove1 <- base::which(df_max_all$bad_run_pcnt >= 5)
    
    writeLines("Step 4: Running loess on df_tmp", log_conn)
    subset_df <- base::subset(df_tmp, current > 0.021 & energy < enrg_cut)
    writeLines(paste("Loess subset rows:", nrow(subset_df)), log_conn)
    
    tmplo  <- stats::loess(energy_bin ~ timeind, data = subset_df, span = 0.0125)
    predlo <- stats::predict(tmplo, newdata = df_tmp$timeind)
    
    writeLines("Step 5: Computing difflo and indices_to_remove2", log_conn)
    difflo <- df_tmp$energy_bin - predlo
    indices_to_remove2 <- base::which(difflo < -2.5)
    
    writeLines("Step 6: Filtering for unphysical energy values", log_conn)
    indices_to_remove3 <- base::which(df_tmp$energy > (enrg_cut + 5))
    
    indices_to_remove <- sort(unique(c(indices_to_remove1, indices_to_remove2, indices_to_remove3)))
    writeLines(paste("Indices to remove:", length(indices_to_remove)), log_conn)
    
    if (length(indices_to_remove) > 0) {
      writeLines("Step 7: Setting omitted columns", log_conn)
      stopifnot(all(c("current_omitted", "energy_omitted", "energy_bin_omitted") %in% colnames(df_tmp)))
      df_tmp$current_omitted[indices_to_remove]    <- df_tmp$current[indices_to_remove]
      df_tmp$energy_omitted[indices_to_remove]     <- df_tmp$energy[indices_to_remove]
      df_tmp$energy_bin_omitted[indices_to_remove] <- df_tmp$energy_bin[indices_to_remove]
      df_tmp$current[indices_to_remove]    <- NA_real_
      df_tmp$energy[indices_to_remove]     <- NA_real_
      df_tmp$energy_bin[indices_to_remove] <- NA_real_
    }
    writeLines(paste(Sys.time(), "Postprocess complete for", day_ymd), log_conn)
    
    # ggplot diagnostics (explicit ggplot2::)
    plot_map(df_long_all, df_curr_all, do_log = TRUE) +
      ggplot2::geom_point(data = subset(df_max_all, energy <= enrg_cut),
                          ggplot2::aes(x = timeind, y = energy_bin),
                          inherit.aes = FALSE, alpha = 0.3) +
      ggplot2::geom_point(data = subset(df_max_all, energy > enrg_cut),
                          ggplot2::aes(x = timeind, y = energy_bin),
                          inherit.aes = FALSE, alpha = 0.3, color = "red")
    ggplot2::ggsave(file.path(process_fig_dir, sprintf("%s_raw.png", day_ymd)),
                    width = 12, height = 9, units = "in", dpi = 300)
    
    plot_energy(df_tmp, df_curr_all, do_log = TRUE, enrg_cut = enrg_cut + 5)
    ggplot2::ggsave(file.path(process_fig_dir, sprintf("%s_energy.png", day_ymd)),
                    width = 8, height = 4, units = "in", dpi = 300)
    
    plot_current(df_tmp, df_curr_all, do_log = FALSE, enrg_cut = enrg_cut + 5)
    ggplot2::ggsave(file.path(process_fig_dir, sprintf("%s_current.png", day_ymd)),
                    width = 8, height = 4, units = "in", dpi = 300)
    
    plot_map(df_long_all, df_curr_all, do_log = TRUE) +
      ggplot2::geom_point(data = df_tmp, ggplot2::aes(x = timeind, y = energy_bin),
                          inherit.aes = FALSE, alpha = 0.3) +
      ggplot2::geom_point(data = df_tmp, ggplot2::aes(x = timeind, y = energy_bin_omitted),
                          inherit.aes = FALSE, alpha = 0.3, color = "red")
    ggplot2::ggsave(file.path(process_fig_dir, sprintf("%s_map.png", day_ymd)),
                    width = 12, height = 9, units = "in", dpi = 300)
    
    # Save data frame
    df_tmp$energy_fitted <- NULL
    save(df_tmp, file = file.path(sv_pathname, sprintf("%s_df.RData", day_ymd)))
    
    return(day_ymd)
    
  }, error = function(e) {
    # Append to the per-day log
    log_file <- file.path(log_dir, sprintf("logs_%s_%s.txt", day_ymd, day_julian))
    err_conn <- file(log_file, open = "a")
    on.exit(try(close(err_conn), silent = TRUE), add = TRUE)
    writeLines(sprintf("%s Error in %s: %s",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       day_ymd, conditionMessage(e)), err_conn)
    NA_character_
  })
}

parallel::stopCluster(cl)



