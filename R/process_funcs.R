#########
######   scaled Vecchia approximation (Katzfuss, Guinness, Lawrence)  #######

#' fit parameters using scaled Vecchia, assuming matern covariance
#' (docstring unchanged)

fit_scaled = function(y, inputs, ms = c(30), trend = 'pre', X, nu = 1.5, nug = 0, scale = 'parms',
                      var.ini, ranges.ini, select = Inf, print.level = 0, max.it = 32, tol.dec = 4,
                      n.est = min(5e3, nrow(inputs)), find.vcf = TRUE, vcf.scorefun = ls) {
  
  ## dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  
  ## specify trend covariates
  if (missing(X)) {
    if (trend == 'zero') {
      X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
    } else if (trend == 'intercept') {
      X = as.matrix(rep(1, n))
    } else if (trend == 'linear') {
      X = cbind(rep(1, n), inputs)
    } else if (trend == 'pre') {
      X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
      beta = mean(y)
      y = y - beta
    } else stop('invalid trend option specified')
  } else trend = 'X'
  
  ## default variance parameter
  if (missing(var.ini)) {
    cur.var = summary(stats::lm(y ~ X - 1))$sigma^2
  } else cur.var = var.ini
  
  ## default range parameters
  input.ranges = apply(inputs, 2, function(x) diff(range(x)))
  if (missing(ranges.ini)) cur.ranges = .2 * input.ranges else cur.ranges = ranges.ini
  active = rep(TRUE, d)
  
  ## fixed nugget?
  if (is.null(nug)) {
    fix.nug = FALSE; nug = .01 * stats::var(y)
  } else fix.nug = TRUE
  
  ## smoothness: fixed? bessel?
  if (is.null(nu)) {
    covfun = 'matern_scaledim'
    cur.oth = c(3.5, nug)
    fix.nu = FALSE
  } else if (nu %in% (.5 + (1:4))) {
    covfun = paste0("matern", nu * 10, "_scaledim")
    cur.oth = nug
    fix.nu = FALSE
  } else {
    covfun = 'matern_scaledim'
    cur.oth = c(nu, nug)
    fix.nu = TRUE
  }
  
  ## only use subsample for estimation?
  if (n.est < n) {
    ind.est = sample(1:n, n.est)
    y.full = y; inputs.full = inputs; X.full = X
    y = y[ind.est]; inputs = inputs[ind.est, , drop = FALSE]; X = X[ind.est, , drop = FALSE]
  }
  
  ## decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y) - 1))
  
  ### for increasing m
  for (i.m in 1:length(ms)) {
    m = ms[i.m]
    if (i.m < length(ms)) { tol = 10^(-tol.dec - 2) } else { tol = 10^(-tol.dec) }
    
    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while (conv == FALSE & maxit <= max.it) {
      
      if (print.level > 0) {
        print(paste0('m=', m, ', maxit=', maxit)); print(cur.ranges)
      }
      
      ## check for inactive input dims (large range params)
      active = (cur.ranges < input.ranges * select)
      if (sum(active, na.rm = TRUE) == 0) stop('all inputs inactive. increase select?')
      cur.ranges[!active] = Inf
      
      ## specify how to scale input dimensions
      cur.ranges[!active] = Inf
      
      ## order and condition based on current params
      if (scale == 'parms') {
        scales = 1 / cur.ranges
      } else if (scale == 'ranges') {
        scales = 1 / input.ranges
      } else if (scale == 'none') {
        scales = 1
      } else stop(paste0('invalid argument scale=', scale))
      
      ## order and condition based on current params
      ord = GPvecchia::order_maxmin_exact(t(t(inputs) * scales))
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      NNarray = GpGp::find_ordered_nn(t(t(inputs.ord) * scales), m)
      
      ## starting and fixed parameters
      cur.parms = c(cur.var, cur.ranges[active], cur.oth)
      fixed = NULL
      if (fix.nu)  fixed = c(fixed, length(cur.parms) - 1)
      if (fix.nug) fixed = c(fixed, length(cur.parms))
      
      ## fisher scoring
      fit = GpGp::fit_model(
        y.ord, inputs.ord[, active, drop = FALSE], X.ord,
        NNarray = NNarray, m_seq = m, convtol = tol,
        start_parms = cur.parms, max_iter = maxit,
        covfun_name = covfun, silent = (print.level < 2),
        reorder = FALSE, fixed_parms = fixed
      )
      cur.var = fit$covparms[1]
      cur.ranges[active] = fit$covparms[1 + (1:sum(active))]
      cur.oth = fit$covparms[-(1:(1 + sum(active)))]
      conv = fit$conv
      maxit = maxit * 2
    }
  }
  
  ### prepare fit object for subsequent prediction
  fit$covparms = c(cur.var, cur.ranges, cur.oth)
  fit$trend = trend
  if (n.est < n) {
    fit$y = y.full
    fit$locs = inputs.full
    fit$X = X.full
  } else {
    fit$locs = inputs.ord
  }
  if (trend == 'zero') {
    fit$X = as.matrix(rep(0, n))
  } else if (trend == 'pre') {
    fit$betahat = beta
    fit$y = fit$y + beta
    fit$trend = 'intercept'
    fit$X = as.matrix(rep(1, n))
  }
  
  ### find variance correction factor, if requested
  if (find.vcf) {
    fit$vcf = fit_vcf(fit, scale = scale, scorefun = vcf.scorefun)
  } else fit$vcf = 1
  
  return(fit)
}


#######   prediction   ########

predictions_scaled <- function(fit, locs_pred, m = 100, joint = TRUE, nsims = 0,
                               predvar = FALSE, X_pred, scale = 'parms') {
  
  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  vcf = if (is.null(fit$vcf)) 1 else fit$vcf
  
  # specify trend if missing
  if (missing(X_pred)) {
    if (fit$trend == 'zero') {
      X_pred = as.matrix(rep(0, n_pred))
    } else if (fit$trend == 'intercept') {
      X_pred = as.matrix(rep(1, n_pred))
    } else if (fit$trend == 'linear') {
      X_pred = cbind(rep(1, n_pred), locs_pred)
    } else stop('X_pred must be specified')
  }
  
  # specify how to scale input dimensions
  if (scale == 'parms') {
    scales = 1 / covparms[1 + (1:ncol(locs_obs))]
  } else if (scale == 'ranges') {
    scales = 1 / apply(locs_obs, 2, function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=', scale))
  
  if (joint) {  # joint predictions
    
    # get orderings
    temp = order_maxmin_pred(t(t(locs_obs) * scales), t(t(locs_pred) * scales))
    ord1 = temp$ord
    ord2 = temp$ord_pred
    
    # reorder stuff
    X_obs  <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1, , drop = FALSE]
    Xord_pred <- X_pred[ord2, , drop = FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta
    
    # put all coordinates together
    locs_all <- rbind(locs_obs[ord1, , drop = FALSE], locs_pred[ord2, , drop = FALSE])
    
    # get nearest neighbor array
    sm = if (n_pred < 1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(t(t(locs_all) * scales), m,
                                        fix.first = n_obs, searchmult = sm)
    NNarray_pred = NNarray_all[-(1:n_obs), -1]
    
    means = numeric(length = n_pred)
    if (nsims > 0) samples = array(dim = c(n_pred, nsims))
    
    # make predictions sequentially
    for (i in 1:n_pred) {
      
      # NN conditioning sets
      NN      = sort(NNarray_pred[i, ])
      NN_obs  = NN[NN <= n_obs]
      NN_pred = NN[NN >  n_obs] - n_obs
      
      # (co-)variances
      K  = get(covfun_name)(covparms, locs_all[c(NN, i + n_obs), ])
      cl = t(chol(K))
      
      # prediction
      y.NN      = y[NN_obs]
      means[i]  = cl[m + 1, 1:m] %*% base::forwardsolve(cl[1:m, 1:m], c(y.NN, means[NN_pred]))
      if (nsims > 0) { # conditional simulation
        pred.var = cl[m + 1, m + 1]^2 * vcf
        for (s in 1:nsims) {
          pm = cl[m + 1, 1:m] %*% base::forwardsolve(cl[1:m, 1:m], c(y.NN, samples[NN_pred, s]))
          samples[i, s] = stats::rnorm(1, pm, sqrt(pred.var))
        }
      }
    }
    
    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if (nsims == 0) {
      preds = means
    } else {
      samples[ord2, ] = samples + c(Xord_pred %*% beta)
      preds = list(means = means, samples = samples)
    }
    
  } else {  # separate predictions
    
    if (nsims > 0) stop('cannot produce joint samples when joint=FALSE')
    
    y  = y_obs - X_obs %*% beta
    
    # find the NNs
    m = min(m, nrow(locs_obs))
    NNarray = FNN::get.knnx(t(t(locs_obs) * scales),
                            t(t(locs_pred) * scales), m)$nn.index
    
    means = vars = numeric(length = n_pred)
    for (i in 1:n_pred) {
      NN = NNarray[i, ]
      K  = get(covfun_name)(covparms, rbind(locs_obs[NN, ], locs_pred[i, ]))
      cl = t(chol(K))
      means[i] = cl[m + 1, 1:m] %*% base::forwardsolve(cl[1:m, 1:m], y[NN])
      vars[i]  = cl[m + 1, m + 1]^2 * vcf
    }
    means = means + c(X_pred %*% beta)
    
    if (predvar == FALSE) {
      preds = means
    } else {
      preds = list(means = means, vars = vars)
    }
  }
  
  return(preds)
}


#######   obs.pred maxmin ordering   ########
order_maxmin_pred <- function(locs, locs_pred, refine = FALSE) {
  
  ord <- 1:nrow(locs)
  ord_pred <- GPvecchia::order_maxmin_exact(locs_pred)
  
  if (refine) {
    locs_all = rbind(locs, locs_pred)
    n <- nrow(locs)
    m <- min(round(sqrt(n)), 200)
    n_pred <- nrow(locs_pred)
    
    NN <- FNN::get.knn(locs_all, k = m)$nn.index
    index_in_position <- c(ord, n + ord_pred, rep(NA, n_pred))
    position_of_index <- order(index_in_position[1:(n + n_pred)])
    
    curlen <- n + n_pred
    nmoved <- 0
    for (j in (n + 1):(n + 2 * n_pred)) {
      nneigh <- round(min(m, 1 * (n + n_pred) / (j - nmoved + 1)))
      neighbors <- NN[index_in_position[j], 1:nneigh]
      if (min(position_of_index[neighbors], na.rm = TRUE) < j) {
        nmoved <- nmoved + 1
        curlen <- curlen + 1
        position_of_index[index_in_position[j]] <- curlen
        index_in_position[curlen] <- index_in_position[j]
        index_in_position[j] <- NA
      }
    }
    ord_pred <- index_in_position[!is.na(index_in_position)][(n + 1):(n + n_pred)] - n
  }
  return(list(ord = ord, ord_pred = ord_pred))
}


#######   find NN for prediction locations   ########
find_ordered_nn_pred <- function(locs, m, fix.first = 0, searchmult = 2) {
  
  if (is.null(ncol(locs))) locs <- as.matrix(locs)
  
  n <- nrow(locs)
  m <- min(m, n - 1)
  mult <- 2
  
  ee <- min(apply(locs, 2, stats::sd))
  locs <- locs + matrix(ee * 1e-6 * stats::rnorm(n * ncol(locs)), n, ncol(locs))
  
  NNarray <- matrix(NA, n, m + 1)
  
  maxval <- min(mult * m + 1, n)
  if (fix.first <= maxval) {
    NNarray[1:maxval, ] <- GpGp::find_ordered_nn_brute(locs[1:maxval, , drop = FALSE], m)
  } else {
    maxval = fix.first
    NNarray[1:(m + 1), ] <- GpGp::find_ordered_nn_brute(locs[1:(m + 1), , drop = FALSE], m)
    NNarray[1:maxval, 1] = 1:maxval
    NNarray[(m + 1):maxval, 1 + (1:m)] = matrix(rep(1:m, maxval - m), byrow = TRUE, ncol = m)
  }
  query_inds <- min(maxval + 1, n):n
  data_inds  <- 1:n
  msearch <- m
  
  while (length(query_inds) > 0) {
    msearch <- min(max(query_inds), round(searchmult * msearch))
    data_inds <- 1:min(max(query_inds), n)
    NN <- FNN::get.knnx(locs[data_inds, , drop = FALSE], locs[query_inds, , drop = FALSE], msearch)$nn.index
    less_than_k <- t(sapply(1:nrow(NN), function(k) NN[k, ] <= query_inds[k]))
    sum_less_than_k <- apply(less_than_k, 1, sum)
    ind_less_than_k <- which(sum_less_than_k >= m + 1)
    
    NN_m <- t(sapply(ind_less_than_k, function(k) NN[k, ][less_than_k[k, ]][1:(m + 1)]))
    NNarray[query_inds[ind_less_than_k], ] <- NN_m
    query_inds <- query_inds[-ind_less_than_k]
  }
  
  return(NNarray)
}


##########   line search for variance correction factor   ###
fit_vcf = function(fit, m.pred = 140, n.test = min(1e3, round(nrow(fit$locs) / 5)),
                   scale = 'parms', scorefun = ls) {
  
  fitsearch = fit
  inds.test = sample(1:nrow(fit$locs), n.test)
  fitsearch$y   = fit$y[-inds.test]
  fitsearch$locs = fit$locs[-inds.test, , drop = FALSE]
  fitsearch$X   = fit$X[-inds.test, , drop = FALSE]
  
  preds = predictions_scaled(
    fitsearch, locs_pred = fit$locs[inds.test, , drop = FALSE],
    m = m.pred, joint = FALSE, predvar = TRUE, scale = scale,
    X_pred = fit$X[inds.test, , drop = FALSE]
  )
  
  y.test = fit$y[inds.test]
  objfun = function(vcf) scorefun(y.test, preds$means, preds$vars * vcf)
  vcf = stats::optimize(objfun, c(1e-6, 1e6))$minimum
  
  return(vcf)
}

### log score
ls = function(dat, mu, sig2) -mean(stats::dnorm(dat, mu, sqrt(sig2), log = TRUE))


######### Helper functions #########

get_timesec = function(time_vec){
  tdiff = c(0, time_vec[2:length(time_vec)] - time_vec[1:(length(time_vec) - 1)])
  return(cumsum(tdiff))
}

load_netcdf_data <- function(nc_file, log_conn, file_level = "L2") {
  if (!file.exists(nc_file)) {
    writeLines(paste(Sys.time(), "File not found:", nc_file), log_conn)
    return(NULL)
  }
  # Load NetCDF
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
  
  energy_df <- as.data.frame(energy) |>
    dplyr::mutate(
      timestamp = time_vec,
      timeind   = 1:length(time_vec),
      timemin   = get_timesec(timestamp) / 60
    )
  
  current_df <- as.data.frame(current) |>
    dplyr::mutate(
      timestamp = time_vec,
      timeind   = 1:length(time_vec),
      timemin   = get_timesec(timestamp) / 60
    )
  
  return(list(
    energy_df  = energy_df,
    current_df = current_df
  ))
}

# Define a function to read Rds, CSV, or TXT files
setup_data = function(file_name){
  file_extension <- tools::file_ext(file_name)
  if (file_extension == "Rds" | file_extension == "rds") {
    df <- readRDS(file_name)
  } else if (file_extension == "csv" | file_extension == "txt") {
    df <- utils::read.csv(file_name, stringsAsFactors = FALSE, header = TRUE)
  } else {
    stop("Unsupported file type. Please provide an Rds, CSV, or TXT file.")
  }
  bin_columns <- grep("^bin", names(df), value = TRUE)
  df <- df[rowSums(is.na(df[, bin_columns])) < length(bin_columns), ]
  df$timeind = 1:nrow(df)
  df$timemin = get_timesec(as.POSIXct(df$timestamp)) / 60
  return(df)
}

unmelt <- function(df_long, id_vars, value.var = "value") {
  # keep only needed columns first
  dt <- data.table::as.data.table(df_long)[, c(id_vars, "variable", value.var), with = FALSE]
  # build a **formula**, not a character string
  form <- as.formula(paste(paste(id_vars, collapse = " + "), "~ variable"))
  # wide cast
  df <- data.table::dcast(dt, form, value.var = value.var)
  # reorder columns: bins first, then id_vars
  non_id_vars <- setdiff(names(df), id_vars)
  df <- df[, c(non_id_vars, id_vars), with = FALSE]
  # (optional) return a plain data.frame to avoid data.table semantics later
  df <- as.data.frame(df)
  return(df)
}


######### Main algorithm functions #########
get_max_df = function(df_curr, df_enrg, energy_bins, bin_inds = 1:length(energy_bins),
                      df_ll = NULL, df_sd = NULL){
  mat = as.matrix(df_curr[, bin_inds])
  df_max = data.frame(
    variable  = factor(rep(NA, nrow(df_curr)), levels = energy_bins),
    current   = NA,
    energy_bin = NA,
    energy    = NA,
    timeind   = df_curr$timeind,
    timestamp = df_curr$timestamp,
    timemin   = df_curr$timemin
  )
  if (!is.null(df_ll)) within_floor = rep(NA, nrow(mat))
  if (!is.null(df_sd)) predsd = rep(NA, nrow(mat))
  
  for (ii in 1:nrow(mat)) {
    bin_ind = which.max(mat[ii, ])
    df_max$energy_bin[ii] = bin_ind
    df_max$energy[ii]     = df_enrg[ii, bin_ind]
    df_max$variable[ii]   = energy_bins[bin_ind]
    df_max$current[ii]    = mat[ii, bin_ind]
    if (!is.null(df_ll)) {
      within_floor[ii] = (df_ll[ii, bin_ind] < 0)
    }
    if (!is.null(df_sd)) {
      predsd[ii] = df_sd[ii, bin_ind]
    }
  }
  if (!is.null(df_ll)) df_max$within_floor = within_floor
  if (!is.null(df_sd)) df_max$predsd = predsd
  return(df_max)
}


get_decr = function(df_long, bad_timeinds, energy_bin_max = 20) {
  tot_dec = rep(NA, length(bad_timeinds))
  for (bb in seq_along(bad_timeinds)) {
    bind  = bad_timeinds[bb]
    tmpdf = subset(df_long, timeind == bind)
    valid = tmpdf$energy_bin <= energy_bin_max
    tmpdiff = diff(tmpdf$predmu)
    valid_diffs = valid[-1]
    tot_dec[bb] = sum(tmpdiff[valid_diffs & tmpdiff < 0])
  }
  return(abs(tot_dec))
}

get_noisefloor_timeinds = function(df_long, df_max, pcnt_bad = 10, bad_bins = 50:100){
  is_bad_time <- sapply(df_max$variable, function(x)
    x %in% c(paste0('bin', bad_bins), paste0('V', bad_bins))
  )
  bad_timeinds = df_max$timeind[is_bad_time]
  tot_dec = get_decr(df_long, bad_timeinds)
  
  keep_going = TRUE; window_width = 0; desired_l = round((pcnt_bad / 100) * sum(is_bad_time))
  noisefloor_timeinds_old = c()
  while (keep_going) {
    window_width = window_width + 1
    roll_wind = zoo::rollsum(tot_dec, k = 2 * window_width + 1, fill = NA, align = "center")
    noisefloor_timeinds = bad_timeinds[roll_wind == 0]; noisefloor_timeinds = noisefloor_timeinds[!is.na(noisefloor_timeinds)]
    stop_cond_1 = (length(noisefloor_timeinds) < desired_l)
    stop_cond_2 = sum(is.na(roll_wind)) == length(roll_wind)
    if (stop_cond_1 | stop_cond_2) {
      keep_going = FALSE
      if (stop_cond_2) { noisefloor_timeinds = noisefloor_timeinds_old }
    } else {
      noisefloor_timeinds_old = noisefloor_timeinds
    }
  }
  return(noisefloor_timeinds)
}

get_background_timeinds = function(df_max, energy_floor = 50, alpha = 0.5, get_all = FALSE, pad = 20){
  vals = subset(df_max, energy > energy_floor)$timeind
  df = data.frame(val = vals)
  
  factoextra::fviz_nbclust(df, stats::kmeans, method = "silhouette") -> optimal_silhouette_plot
  optimal_k <- as.numeric(optimal_silhouette_plot$data$clusters[which.max(optimal_silhouette_plot$data$y)])
  
  kmeans_result <- stats::kmeans(df$val, centers = as.numeric(optimal_k), nstart = 25)
  df$cluster <- as.factor(kmeans_result$cluster)
  
  if (get_all) {
    meds = df |>
      dplyr::group_by(cluster) |>
      dplyr::summarize(
        med = stats::median(val),
        q1  = stats::quantile(val, 0.1),
        q2  = stats::quantile(val, 0.9),
        .groups = "drop"
      )
    df$resid = NA_real_
    for (ii in 1:nrow(df)) {
      df$resid[ii] = df$val[ii] - as.numeric(meds[as.numeric(df$cluster[ii]), 'med'])
    }
    resid_cutoff = max(meds$q2 - meds$q1)
    df_sub = subset(df, abs(resid) < resid_cutoff)
    cluster_summary = df_sub |>
      dplyr::group_by(cluster) |>
      dplyr::summarize(
        V1 = min(val),
        V2 = max(val),
        .groups = "drop"
      )
    mxtimeind = max(df_max$timeind)
    pads = (pad / 100) * (cluster_summary$V2 - cluster_summary$V1)
    cluster_summary$V1 = sapply(cluster_summary$V1 - pads, function(x) max(x, 1))
    cluster_summary$V2 = sapply(cluster_summary$V2 + pads, function(x) min(x, mxtimeind))
  } else {
    cluster_sizes <- df |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(size = dplyr::n(), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(size))
    
    if (nrow(cluster_sizes) >= 2) {
      largest_clusters <- cluster_sizes$cluster[1:2]
      df <- df |>
        dplyr::filter(cluster %in% largest_clusters)
      optimal_k = 2
    }
    
    cluster_summary <- df |>
      dplyr::group_by(cluster) |>
      dplyr::summarize(
        V1 = stats::quantile(val, 0.5 - alpha / 2),
        V2 = stats::quantile(val, 0.5 + alpha / 2),
        .groups = "drop"
      )
  }
  
  inds_bad = c()
  for (jj in 1:optimal_k) {
    tmp_rng = round(unname(unlist(cluster_summary[jj, -1])))
    inds_bad = c(inds_bad, tmp_rng[1]:tmp_rng[2])
  }
  return(inds_bad)
}

# Fit an optimal parametric function to get the noise floor.
get_optim_noisefloor <- function(data_profile, richards_energy_bins = 1:97, gaussian_energy_bins = 1:20,
                                 parabolic_adjustment = TRUE, return_gaussian = FALSE,
                                 max_iter = 1000, min_iter = 100, tol = 1e-6, verbose = FALSE) {
  
  f_richards <- function(x, A, k, x0, nu, A0) {
    A0 + A / (1 + exp(-k * (x - x0)))^(1 / nu)
  }
  f_gaussian <- function(x, amp, mu, sigma) {
    amp * exp(-((x - mu)^2) / (2 * sigma^2))
  }
  f_parabola <- function(x, p, q, r) {
    p * x^2 + q * x + r
  }
  
  sse_richards <- function(par, x, data) {
    sum((data - f_richards(x, par[1], par[2], par[3], par[4], par[5]))^2)
  }
  sse_gaussian <- function(par, x, resid) {
    sum((resid - f_gaussian(x, par[1], par[2], par[3]))^2)
  }
  sse_parabola <- function(par, x, resid) {
    sum((resid - f_parabola(x, par[1], par[2], par[3]))^2)
  }
  
  richards_par <- c(A = 0.14, k = 0.1, x0 = 10, nu = 1, A0 = 0)
  
  prof_diff = diff(data_profile)
  neg_bins = which(prof_diff < 0 & abs(prof_diff) > 0.003) + 1; neg_bins = neg_bins[neg_bins < 25]
  if (length(neg_bins) > 0) {
    contaminated_bins = sort(c(neg_bins, (min(neg_bins) - 1):(min(neg_bins) - 1 - length(neg_bins))))
    contaminated_bins = contaminated_bins[contaminated_bins > min(richards_energy_bins)]
    if (length(contaminated_bins) > 2) {
      contaminated_bins = contaminated_bins[2:(length(contaminated_bins) - 1)]
    }
  } else { contaminated_bins = c() }
  
  res_richards <- stats::optim(
    par = richards_par, fn = sse_richards,
    x = setdiff(richards_energy_bins, contaminated_bins),
    data = data_profile[setdiff(richards_energy_bins, contaminated_bins)],
    method = "L-BFGS-B",
    lower = c(0, 0, 0, 0, -Inf),
    upper = c(Inf, 1, Inf, Inf, Inf)
  )
  richards_par <- res_richards$par
  
  gaussian_par_init <- list(
    c(amp = 0.00, mu = 10, sigma = 2),
    c(amp = 0.01, mu = 5,  sigma = 2),
    c(amp = 0.01, mu = 10, sigma = 2),
    c(amp = 0.01, mu = 15, sigma = 2)
  )
  gaussian_started <- FALSE
  
  parabola_par <- c(p = 0, q = 0, r = 0)
  parabola_fit <- rep(0, 100)
  
  prev_richards_par <- richards_par
  prev_parabola_par <- parabola_par
  prev_gaussian_par <- c(amp = 0, mu = 10, sigma = 3)
  parabola_started <- FALSE
  
  for (iter in 1:max_iter) {
    richards_fit <- f_richards(1:100, richards_par[1], richards_par[2], richards_par[3], richards_par[4], richards_par[5])
    parabola_fit <- if (parabola_started) f_parabola(1:100, parabola_par[1], parabola_par[2], parabola_par[3]) else rep(0, 100)
    
    residuals_for_gaussian <- data_profile - richards_fit - parabola_fit
    if (iter == 1) {
      best_sse <- Inf
      for (ll in 1:length(gaussian_par_init)) {
        res_gaussian <- stats::optim(
          par = gaussian_par_init[[ll]],
          fn = sse_gaussian,
          x = gaussian_energy_bins,
          resid = residuals_for_gaussian[gaussian_energy_bins],
          method = "L-BFGS-B",
          lower = c(0, min(gaussian_energy_bins), 1e-2),
          upper = c(Inf, max(gaussian_energy_bins), 3)
        )
        tmp_gaussian_par <- res_gaussian$par
        tmp_sse <- sse_gaussian(tmp_gaussian_par, gaussian_energy_bins, residuals_for_gaussian[gaussian_energy_bins])
        if (tmp_sse < best_sse) {
          gaussian_par <- tmp_gaussian_par
          best_sse <- tmp_sse
        }
      }
    } else {
      res_gaussian <- stats::optim(
        par = gaussian_par,
        fn = sse_gaussian,
        x = gaussian_energy_bins,
        resid = residuals_for_gaussian[gaussian_energy_bins],
        method = "L-BFGS-B",
        lower = c(0, min(gaussian_energy_bins), 1e-2),
        upper = c(Inf, max(gaussian_energy_bins), 3)
      )
      gaussian_par <- res_gaussian$par
    }
    
    amp <- gaussian_par[1]
    if (iter <= 10 && amp > 0) {
      mu <- gaussian_par[2]
      sigma <- gaussian_par[3]
      mask_region <- (richards_energy_bins < (mu - 3 * sigma)) | (richards_energy_bins > (mu + 3 * sigma))
      masked_bins <- richards_energy_bins[mask_region]
      if (length(masked_bins) < 5) masked_bins <- richards_energy_bins
    } else {
      masked_bins <- richards_energy_bins
    }
    
    adjusted_data <- data_profile - f_gaussian(1:100, gaussian_par[1], gaussian_par[2], gaussian_par[3]) - parabola_fit
    res_richards <- stats::optim(
      par = richards_par,
      fn = sse_richards,
      x = masked_bins,
      data = adjusted_data[masked_bins],
      method = "L-BFGS-B",
      lower = c(0, 0, 0, 0, -Inf),
      upper = c(Inf, 1, Inf, Inf, Inf)
    )
    richards_par <- res_richards$par
    
    delta_richards <- sum(abs(richards_par - prev_richards_par))
    delta_gaussian <- sum(abs(gaussian_par - prev_gaussian_par))
    
    if (!parabola_started && delta_richards < tol && delta_gaussian < tol) {
      if (verbose) cat(sprintf("Parabola activated at iteration %d.\n", iter))
      parabola_started <- TRUE
    }
    
    if (parabola_started && parabolic_adjustment) {
      residuals_for_parabola <- data_profile - richards_fit - f_gaussian(1:100, gaussian_par[1], gaussian_par[2], gaussian_par[3])
      res_parabola <- stats::optim(
        par = parabola_par,
        fn = sse_parabola,
        x = masked_bins,
        resid = residuals_for_parabola[masked_bins],
        method = "L-BFGS-B"
      )
      parabola_par <- res_parabola$par
    }
    
    if (parabola_started) {
      delta_parabola <- sum(abs(parabola_par - prev_parabola_par))
      if (verbose) {
        cat(sprintf("Iteration %d:\n", iter))
        cat("  Richards params:", richards_par, "\n")
        cat("  Gaussian params:", gaussian_par, "\n")
        cat("  Parabola params:", parabola_par, "\n")
        cat(sprintf("  ΔRichards: %.6f, ΔGaussian: %.6f, ΔParabola: %.6f\n\n",
                    delta_richards, delta_gaussian, delta_parabola))
      }
      if (delta_richards < tol && delta_gaussian < tol && delta_parabola < tol) {
        if (verbose) cat("Convergence achieved.\n")
        break
      }
    }
    
    prev_richards_par <- richards_par
    prev_gaussian_par <- gaussian_par
    prev_parabola_par <- parabola_par
  }
  
  final_richards_fit <- f_richards(1:100, richards_par[1], richards_par[2], richards_par[3], richards_par[4], richards_par[5])
  final_parabola_fit <- if (parabolic_adjustment && parabola_started) f_parabola(1:100, parabola_par[1], parabola_par[2], parabola_par[3]) else rep(0, 100)
  
  floor_profile <- final_richards_fit + final_parabola_fit
  print(class(floor_profile))
  
  if (return_gaussian) {
    final_gaussian_fit <- f_gaussian(1:100, gaussian_par[1], gaussian_par[2], gaussian_par[3])
    return(list(
      floor_profile     = floor_profile,
      richards_profile  = final_richards_fit,
      gaussian_profile  = final_gaussian_fit,
      parabola_profile  = final_parabola_fit,
      richards_params   = richards_par,
      gaussian_params   = gaussian_par,
      parabola_params   = parabola_par
    ))
  } else {
    return(list(
      floor_profile     = floor_profile,
      parabola_profile  = final_parabola_fit,
      richards_profile  = final_richards_fit,
      richards_params   = richards_par,
      parabola_params   = parabola_par
    ))
  }
}

get_omit_windowinds = function(df_out, threshold_curr = 0.02, desired_window_size = 3000){
  df_out <- df_out |>
    dplyr::mutate(is_near_zero = ifelse(current <= threshold_curr, 1, 0))
  
  runs <- base::rle(df_out$is_near_zero)
  runs_df <- data.frame(lengths = runs$lengths, values = runs$values)
  
  runs_df <- runs_df |>
    dplyr::mutate(
      start_index = cumsum(dplyr::lag(lengths, default = 0)) + 1,
      end_index   = cumsum(lengths)
    ) |>
    dplyr::filter(values == 1 & lengths > desired_window_size)
  
  if (nrow(runs_df) > 0) {
    indices_to_remove <- c(unlist(apply(runs_df, 1, function(x) seq(x["start_index"], x["end_index"]))))
  } else { indices_to_remove <- c() }
  return(indices_to_remove)
}


######### Plotting functions #########

plot_map = function(df_long, df, type = 'current', do_log = FALSE, df_max = NULL,
                    inds_vis = round(seq(1, nrow(df), length.out = 4)),
                    save_fnm = "", width = 15, height = 10, do_fittednofloor = FALSE){
  
  if (do_log) {
    trans  = function(x){ log10(x) }
    fillab = paste0('log10(', type, ')')
  } else {
    trans  = function(x){ x }
    fillab = type
  }
  
  if (!do_fittednofloor) {
    g = ggplot2::ggplot(df_long, ggplot2::aes(x = timeind, y = variable, fill = trans(value)))
  } else {
    g = ggplot2::ggplot(df_long, ggplot2::aes(x = timeind, y = variable, fill = trans(predmu_nofloor)))
  }
  
  g = g +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(y = 'energy bin', fill = fillab, x = 'time') +
    ggplot2::scale_x_continuous(
      breaks = df$timeind[inds_vis],
      labels = df$timestamp[inds_vis]
    )
  
  if (!is.null(df_max)) {
    g = g + ggplot2::geom_point(
      data = df_max,
      ggplot2::aes(x = timeind, y = variable),
      color = 'red', inherit.aes = FALSE, alpha = 0.3
    )
  }
  
  if (!save_fnm == "") {
    ggplot2::ggsave(g, filename = save_fnm, width = width, height = height)
  }
  return(g)
}

plot_energy = function(df_max, df, do_log = FALSE, save_fnm = "", width = 15, height = 6,
                       inds_vis = round(seq(1, nrow(df), length.out = 4)), enrg_cut = NULL){
  
  if (do_log) {
    trans  = function(x){ log10(x) }
    fillab = paste0('log10(current)')
  } else {
    trans  = function(x){ x }
    fillab = 'current'
  }
  
  if (!is.null(enrg_cut)) { df_max$within_floor = df_max$energy > enrg_cut }
  
  if ('within_floor' %in% colnames(df_max)) {
    g = ggplot2::ggplot(
      subset(df_max, !within_floor),
      ggplot2::aes(x = timeind, y = energy, color = trans(current))
    ) +
      ggplot2::geom_point(data = subset(df_max, within_floor), alpha = 0.3, color = 'red') +
      ggplot2::geom_point(alpha = 0.5)
  } else {
    g = ggplot2::ggplot(
      df_max,
      ggplot2::aes(x = timeind, y = energy, color = trans(current))
    ) +
      ggplot2::geom_point(alpha = 0.5)
  }
  
  g = g +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(y = 'energy at peak current', color = fillab, x = 'time') +
    ggplot2::scale_x_continuous(
      breaks = df$timeind[inds_vis],
      labels = df$timestamp[inds_vis]
    )
  
  if (!save_fnm == "") {
    ggplot2::ggsave(g, filename = save_fnm, width = width, height = height)
  }
  return(g)
}

plot_current = function(df_max, df, do_log = FALSE, save_fnm = "", width = 15, height = 6,
                        inds_vis = round(seq(1, nrow(df), length.out = 4)), enrg_cut = NULL){
  
  if (do_log) {
    trans  = function(x){ log10(x) }
    fillab = paste0('log10(peak current)')
  } else {
    trans  = function(x){ x }
    fillab = 'peak current'
  }
  
  if (!is.null(enrg_cut)) { df_max$within_floor = df_max$energy > enrg_cut }
  
  if ('within_floor' %in% colnames(df_max)) {
    g = ggplot2::ggplot(
      subset(df_max, !within_floor),
      ggplot2::aes(x = timeind, y = trans(current), color = energy)
    ) +
      ggplot2::geom_point(data = subset(df_max, within_floor), alpha = 0.3, color = 'red') +
      ggplot2::geom_point(alpha = 0.5)
  } else {
    g = ggplot2::ggplot(
      df_max,
      ggplot2::aes(x = timeind, y = trans(current), color = energy)
    ) +
      ggplot2::geom_point(alpha = 0.5)
  }
  
  g = g +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(y = fillab, color = 'energy', x = 'time') +
    ggplot2::scale_x_continuous(
      breaks = df$timeind[inds_vis],
      labels = df$timestamp[inds_vis]
    )
  
  if (!save_fnm == "") {
    ggplot2::ggsave(g, filename = save_fnm, width = width, height = height)
  }
  return(g)
}

