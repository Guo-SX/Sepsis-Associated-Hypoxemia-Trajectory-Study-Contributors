
























pkgs <- c("dplyr","tibble","purrr","gbmt","ggplot2","tidyr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

stopifnot(exists("imp_out"))
stopifnot(all(c("day","stay_id_code") %in% names(imp_out)))

dat0 <- imp_out %>%
  mutate(
    day = as.integer(day),
    stay_id_code = as.integer(stay_id_code)
  ) %>%
  arrange(stay_id_code, day)

if (!("ne_total_mcg" %in% names(dat0))) {
  stop("imp_out 缺少 ne_total_mcg 列（你说有肾上腺素用量）。如果列名不同请改脚本。")
}
dat0 <- dat0 %>%
  mutate(ne_total_mcg_log = log1p(ne_total_mcg))

pad_vectors <- function(vector_list, row_names) {
  max_len <- max(sapply(vector_list, length))
  out <- as.data.frame(do.call(rbind, lapply(vector_list, function(x) {
    c(x, rep(NA, max_len - length(x)))
  })))
  rownames(out) <- row_names
  colnames(out) <- paste0("Class", seq_len(max_len))
  out
}

calculate_OCC <- function(appa, prior) {
  (appa / (1 - appa)) / (prior / (1 - prior))
}

get_first_nonnull <- function(x, default = NA) {
  if (!is.null(x)) return(x)
  default
}

get_ic_value <- function(ic, key) {
  if (is.null(ic)) return(NA_real_)
  nm <- names(ic)
  idx <- which(tolower(nm) == tolower(key))
  if (length(idx) == 1) return(as.numeric(ic[[idx]]))
  NA_real_
}

get_npm <- function(m) {
  npm <- NA_real_
  if (!is.null(m$npm) && length(m$npm) > 0) npm <- suppressWarnings(as.numeric(m$npm))
  if (is.na(npm) && !is.null(m$np) && length(m$np) > 0) npm <- suppressWarnings(as.numeric(m$np))
  npm
}

get_loglik <- function(m) {
  ll <- NA_real_
  if (!is.null(m$loglik) && length(m$loglik) > 0) ll <- suppressWarnings(as.numeric(m$loglik))
  if (is.na(ll) && !is.null(m$logLik) && length(m$logLik) > 0) ll <- suppressWarnings(as.numeric(m$logLik))
  if (is.na(ll)) {
    aic <- get_ic_value(m$ic, "AIC")
    npm <- get_npm(m)
    if (!is.na(aic) && !is.na(npm)) ll <- -(aic - 2*npm)/2
  }
  ll
}

entropy01 <- function(m) {
  P <- NULL
  if (!is.null(m$pprob)) P <- as.matrix(m$pprob)
  if (is.null(P) && !is.null(m$posterior)) P <- as.matrix(m$posterior)
  if (is.null(P) && !is.null(m$post)) P <- as.matrix(m$post)
  if (is.null(P)) return(NA_real_)

  rs <- rowSums(P)
  if (any(rs == 0) || any(!is.finite(rs))) return(NA_real_)
  P <- P / rs
  P[P <= 0] <- NA_real_

  n <- nrow(P); G <- ncol(P)
  if (n == 0 || G <= 1) return(NA_real_)

  ent_raw  <- -rowSums(P * log(P), na.rm = TRUE)
  ent_mean <- mean(ent_raw, na.rm = TRUE)

  1 - (ent_mean / log(G))
}

class_pct6 <- function(assign, G) {
  out <- rep(NA_real_, 6)
  if (!is.null(assign) && length(assign) > 0 && !is.na(G) && G >= 1) {
    tb <- prop.table(table(assign)) * 100
    for (g in seq_len(min(G, 6))) {
      if (as.character(g) %in% names(tb)) out[g] <- as.numeric(tb[[as.character(g)]])
    }
  }
  setNames(out, paste0("%class", 1:6))
}

fit_one_multistart <- function(dat, varNames,
                               id_col = "stay_id_code",
                               time_col = "day",
                               d, ng,
                               n_start = 20, seed_base = 202601,
                               scaling = 0, verbose = TRUE) {

  fits <- vector("list", n_start)
  ics  <- vector("list", n_start)

  for (i in seq_len(n_start)) {
    set.seed(seed_base + i)
    if (verbose) {
      message(sprintf("[Fit] d=%d, K=%d, start=%d/%d (seed=%d)",
                      d, ng, i, n_start, seed_base + i))
    }

    fits[[i]] <- gbmt(
      x.names = varNames,
      unit    = id_col,
      time    = time_col,
      d       = d,
      ng      = ng,
      data    = dat,
      scaling = scaling
    )
    ics[[i]] <- fits[[i]]$ic
  }

  IC_all <- as.data.frame(do.call(rbind, ics))
  rownames(IC_all) <- paste0("start_", seq_len(n_start))

  if ("BIC" %in% colnames(IC_all)) {
    best_i <- which.min(IC_all$BIC)
  } else if ("bic" %in% colnames(IC_all)) {
    best_i <- which.min(IC_all$bic)
  } else {
    best_i <- which.min(IC_all[, 1])
  }

  list(
    best_model = fits[[best_i]],
    best_start = best_i,
    IC_all     = IC_all,
    all_models = fits
  )
}

fit_grid <- function(dat, varNames,
                     id_col = "stay_id_code", time_col = "day",
                     d_vec = c(1,2,3), ng_vec = 1:6,
                     n_start = 20, seed_base = 202601,
                     scaling = 0, verbose = TRUE) {

  grid <- list()
  summary_rows <- list()
  seed_cursor <- seed_base

  for (d in d_vec) {
    for (ng in ng_vec) {

      ms <- fit_one_multistart(
        dat = dat, varNames = varNames,
        id_col = id_col, time_col = time_col,
        d = d, ng = ng,
        n_start = n_start,
        seed_base = seed_cursor,
        scaling = scaling,
        verbose = verbose
      )

      key <- paste0("d", d, "_K", ng)
      grid[[key]] <- ms

      ic_best <- ms$best_model$ic
      ic_best_row <- as.data.frame(t(ic_best))
      ic_best_row$d <- d
      ic_best_row$ng <- ng
      ic_best_row$best_start <- ms$best_start
      ic_best_row$key <- key
      summary_rows[[key]] <- ic_best_row

      seed_cursor <- seed_cursor + n_start + 100
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows) %>%
    relocate(key, d, ng, best_start)

  list(grid = grid, summary = summary_df)
}

extract_quality2 <- function(m, dat_use, id_col = "stay_id_code") {
  appa  <- get_first_nonnull(m$appa, NA)
  prior <- get_first_nonnull(m$prior, NA)
  occ   <- if (all(is.na(appa)) || all(is.na(prior))) NA else calculate_OCC(appa, prior)
  prop  <- prop.table(table(m$assign))

  n_id <- length(unique(dat_use[[id_col]]))
  K    <- length(unique(m$assign))

  post_mat <- NULL
  if (!is.null(m$pprob)) post_mat <- m$pprob
  if (is.null(post_mat) && !is.null(m$posterior)) post_mat <- m$posterior
  if (is.null(post_mat) && !is.null(m$post)) post_mat <- m$post

  list(
    n_id = n_id, K = K,
    appa = appa, prior = prior, occ = occ, prop = prop,
    loglik = get_loglik(m),
    entropy = entropy01(m),
    AIC = get_ic_value(m$ic, "AIC"),
    BIC = get_ic_value(m$ic, "BIC"),
    SABIC = dplyr::coalesce(get_ic_value(m$ic, "SSBIC"),
                            get_ic_value(m$ic, "SABIC"),
                            get_ic_value(m$ic, "ssbic"),
                            get_ic_value(m$ic, "sabic")),
    posterior = post_mat
  )
}

save_group_and_means_strict <- function(dat_use, best_model, varNames,
                                        out_prefix,
                                        id_col = "stay_id_code",
                                        time_col = "day") {

  ids <- sort(unique(dat_use[[id_col]]))
  if (length(ids) != length(best_model$assign)) {
    stop(sprintf("ID数量(%d) 与 best_model$assign长度(%d) 不一致：请确认 unit/id_col 完全一致。",
                 length(ids), length(best_model$assign)))
  }

  group_df <- tibble::tibble(
    !!id_col := ids,
    group = as.integer(best_model$assign)
  )
  write.csv(group_df, paste0(out_prefix, "_group.csv"), row.names = FALSE)

  dat2 <- dat_use %>%
    mutate(!!id_col := as.integer(.data[[id_col]]))

  group_df2 <- group_df %>%
    mutate(!!id_col := as.integer(.data[[id_col]]))

  dat_with_group <- dat2 %>%
    left_join(group_df2, by = id_col)

  if (any(is.na(dat_with_group$group))) {
    stop("join 后出现 group=NA：通常是ID类型/取值不一致导致。")
  }

  group_day_means <- dat_with_group %>%
    group_by(group, .data[[time_col]]) %>%
    summarise(
      n_obs = dplyr::n(),
      across(all_of(varNames), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(group, .data[[time_col]])

  write.csv(group_day_means, paste0(out_prefix, "_group_day_means.csv"), row.names = FALSE)

  list(group_df = group_df, dat_with_group = dat_with_group, group_day_means = group_day_means)
}

run_one_scenario <- function(dat0,
                             scenario_name,
                             lactate_col,
                             sf_col,
                             pf_col,
                             include_ne = TRUE,
                             ne_col = "ne_total_mcg_log",
                             id_col = "stay_id_code",
                             time_col = "day",
                             d_vec = c(1,2,3),
                             ng_vec = 1:6,
                             n_start = 20,
                             seed_base = 202601,
                             scaling = 0,
                             verbose = TRUE,
                             thr_appa = 0.70,
                             thr_occ  = 5,
                             thr_prop = 0.05) {

  varNames <- c(lactate_col, sf_col, pf_col)
  if (include_ne) varNames <- c(varNames, ne_col)

  need_cols <- c(id_col, time_col, varNames)
  miss <- setdiff(need_cols, names(dat0))
  if (length(miss) > 0) stop("缺少这些列：", paste(miss, collapse = ", "))

  dat_use <- dat0 %>%
    select(all_of(need_cols)) %>%
    arrange(.data[[id_col]], .data[[time_col]]) %>%
    as.data.frame()

  res <- fit_grid(
    dat = dat_use,
    varNames = varNames,
    id_col = id_col,
    time_col = time_col,
    d_vec = d_vec,
    ng_vec = ng_vec,
    n_start = n_start,
    seed_base = seed_base,
    scaling = scaling,
    verbose = verbose
  )

  quality_list <- lapply(names(res$grid), function(k) {
    ms <- res$grid[[k]]
    q <- extract_quality2(ms$best_model, dat_use, id_col = id_col)
    q$key <- k
    q
  })
  names(quality_list) <- sapply(quality_list, `[[`, "key")

  quality_summary <- bind_rows(lapply(quality_list, function(q) {
    tibble(
      key = q$key,
      n_id = q$n_id,
      K = q$K,
      min_AvePP  = if (all(is.na(q$appa))) NA else min(q$appa, na.rm = TRUE),
      mean_AvePP = if (all(is.na(q$appa))) NA else mean(q$appa, na.rm = TRUE),
      min_OCC    = if (all(is.na(q$occ)))  NA else min(q$occ,  na.rm = TRUE),
      min_prop   = if (all(is.na(q$prop))) NA else min(as.numeric(q$prop), na.rm = TRUE),
      max_prop   = if (all(is.na(q$prop))) NA else max(as.numeric(q$prop), na.rm = TRUE),
      loglik  = q$loglik,
      entropy = q$entropy,
      AIC = q$AIC, BIC = q$BIC, SABIC = q$SABIC
    )
  }))

  model_choice <- res$summary %>%
    left_join(quality_summary, by = "key") %>%
    mutate(
      BIC_use = dplyr::coalesce(.data$BIC, .data$bic, .data$BIC),
      AIC_use = dplyr::coalesce(.data$AIC, .data$aic, .data$AIC)
    )

  cand <- model_choice %>%
    filter(
      is.na(min_AvePP) | min_AvePP >= thr_appa,
      is.na(min_OCC)   | min_OCC   >  thr_occ,
      is.na(min_prop)  | min_prop  >= thr_prop
    ) %>%
    arrange(BIC_use) %>%
    mutate(
      delta_BIC = BIC_use - min(BIC_use, na.rm = TRUE),
      delta_AIC = AIC_use - min(AIC_use, na.rm = TRUE)
    )

  if (nrow(cand) > 0) {
    best_key <- cand$key[which.min(cand$BIC_use)]
    pick_note <- "picked from candidates (constraints satisfied)"
  } else {
    best_key <- model_choice$key[which.min(model_choice$BIC_use)]
    pick_note <- "fallback: picked from all models by min BIC"
  }

  best_model <- res$grid[[best_key]]$best_model

  prefix <- paste0("GBMTM_", scenario_name)

  write.csv(model_choice, paste0(prefix, "_model_selection_summary.csv"), row.names = FALSE)
  write.csv(cand,        paste0(prefix, "_candidate_models.csv"), row.names = FALSE)

  to_long_vec <- function(v, key, type) {
    if (all(is.na(v))) return(tibble(key=key, type=type, class=NA_integer_, value=NA_real_))
    tibble(key = key, type = type, class = seq_along(v), value = as.numeric(v))
  }
  appa_long <- bind_rows(lapply(quality_list, function(q) to_long_vec(q$appa, q$key, "AvePP")))
  occ_long  <- bind_rows(lapply(quality_list, function(q) to_long_vec(q$occ,  q$key, "OCC")))
  prop_long <- bind_rows(lapply(quality_list, function(q) to_long_vec(q$prop, q$key, "Prop")))
  quality_long <- bind_rows(appa_long, occ_long, prop_long)
  write.csv(quality_long, paste0(prefix, "_quality_long.csv"), row.names = FALSE)

  appa_tbl <- pad_vectors(lapply(quality_list, `[[`, "appa"), row_names = names(quality_list))
  occ_tbl  <- pad_vectors(lapply(quality_list, `[[`, "occ"),  row_names = names(quality_list))
  prop_tbl <- pad_vectors(lapply(quality_list, `[[`, "prop"), row_names = names(quality_list))
  write.csv(appa_tbl, paste0(prefix, "_AvePP_wide.csv"))
  write.csv(occ_tbl,  paste0(prefix, "_OCC_wide.csv"))
  write.csv(prop_tbl, paste0(prefix, "_Prop_wide.csv"))

  plot_file <- paste0(prefix, "_BEST_", best_key, "_traj.png")
  png(plot_file, width = 1600, height = 950, res = 150)
  plot(best_model, bands = FALSE,
       xlab = "Day", ylab = "Value",
       titles = varNames)
  dev.off()

  out_g <- save_group_and_means_strict(
    dat_use = as_tibble(dat_use),
    best_model = best_model,
    varNames = varNames,
    out_prefix = paste0(prefix, "_BEST_", best_key),
    id_col = id_col,
    time_col = time_col
  )

  post_mat <- NULL
  if (exists("posterior", mode = "function")) {
    try({
      pm <- posterior(best_model)
      if (!is.null(pm)) post_mat <- pm
    }, silent = TRUE)
  }
  if (is.null(post_mat) && ("posterior" %in% names(best_model))) post_mat <- best_model$posterior
  if (is.null(post_mat) && ("pprob" %in% names(best_model)))     post_mat <- best_model$pprob
  if (is.null(post_mat) && ("post" %in% names(best_model)))      post_mat <- best_model$post

  if (!is.null(post_mat)) {
    post_df <- as.data.frame(post_mat)
    post_df <- post_df %>%
      mutate(!!id_col := out_g$group_df[[id_col]]) %>%
      relocate(!!id_col)
    write.csv(post_df, paste0(prefix, "_BEST_", best_key, "_posterior.csv"), row.names = FALSE)
  }

  cat("\n====================================================\n")
  cat("[Scenario] ", scenario_name, "\n", sep = "")
  cat("Vars: ", paste(varNames, collapse = ", "), "\n", sep = "")
  cat("Best key: ", best_key, " (", pick_note, ")\n", sep = "")
  cat("Saved: \n",
      " - ", paste0(prefix, "_model_selection_summary.csv"), "\n",
      " - ", paste0(prefix, "_candidate_models.csv"), "\n",
      " - ", plot_file, "\n",
      " - ", paste0(prefix, "_BEST_", best_key, "_group.csv"), "\n",
      " - ", paste0(prefix, "_BEST_", best_key, "_group_day_means.csv"), "\n",
      sep = "")
  cat("====================================================\n")

  invisible(list(
    scenario = scenario_name,
    varNames = varNames,
    res = res,
    model_choice = model_choice,
    cand = cand,
    best_key = best_key,
    best_model = best_model
  ))
}

scenarios <- tibble::tribble(
  ~scenario_name,                ~lactate_col,    ~sf_col,    ~pf_col,
  "LAC_MIN_SF_PF_MIN",           "lactate_min",   "sf_min",   "pf_min",
  "LAC_MEAN_SF_PF_MIN",          "lactate_mean",  "sf_min",   "pf_min",
  "LAC_MAX_SF_PF_MIN",           "lactate_max",   "sf_min",   "pf_min",
  "LAC_MIN_SF_PF_MEAN",          "lactate_min",   "sf_mean",  "pf_mean",
  "LAC_MEAN_SF_PF_MEAN",         "lactate_mean",  "sf_mean",  "pf_mean",
  "LAC_MAX_SF_PF_MEAN",          "lactate_max",   "sf_mean",  "pf_mean"
)

need_any <- unique(c("lactate_min","lactate_mean","lactate_max",
                     "sf_min","sf_mean","pf_min","pf_mean",
                     "ne_total_mcg_log","stay_id_code","day"))
miss2 <- setdiff(need_any, names(dat0))
if (length(miss2) > 0) {
  stop("imp_out/dat0 缺少这些关键列：", paste(miss2, collapse = ", "))
}

all_results <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  all_results[[i]] <- run_one_scenario(
    dat0 = dat0,
    scenario_name = sc$scenario_name,
    lactate_col = sc$lactate_col,
    sf_col = sc$sf_col,
    pf_col = sc$pf_col,
    include_ne = TRUE,
    ne_col = "ne_total_mcg_log",
    id_col = "stay_id_code",
    time_col = "day",
    d_vec = c(1,2,3),
    ng_vec = 1:6,
    n_start = 20,
    seed_base = 202601 + i*10000,
    scaling = 0,
    verbose = TRUE,
    thr_appa = 0.70,
    thr_occ  = 5,
    thr_prop = 0.05
  )
}

names(all_results) <- scenarios$scenario_name

best_summary <- bind_rows(lapply(all_results, function(x) {
  mc <- x$model_choice
  bk <- x$best_key
  row <- mc %>% filter(key == bk) %>% slice(1)
  tibble(
    scenario = x$scenario,
    best_key = bk,
    d = row$d,
    K = row$ng,
    BIC = row$BIC_use,
    AIC = row$AIC_use,
    SABIC = row$SABIC,
    min_AvePP = row$min_AvePP,
    min_OCC   = row$min_OCC,
    min_prop  = row$min_prop,
    entropy   = row$entropy
  )
})) %>% arrange(BIC)

write.csv(best_summary, "GBMTM_BEST_OF_6_summary.csv", row.names = FALSE)
print(best_summary)
cat("\n已输出：GBMTM_BEST_OF_6_summary.csv\n")
