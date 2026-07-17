options(warn = 1)

parse_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(FALSE), value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[[length(hit)]])
}

parse_int_vec <- function(name, default) {
  value <- gsub("\\s+", "", parse_arg(name, default))
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    x <- as.integer(strsplit(value, ":", fixed = FALSE)[[1]])
    return(seq(x[1], x[2]))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

append_csv <- function(row, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.table(
    row,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file),
    qmethod = "double"
  )
}

adjusted_rand_index <- function(x, y) {
  x <- as.integer(as.factor(x))
  y <- as.integer(as.factor(y))
  if (length(x) != length(y)) stop("ARI vectors have different lengths.", call. = FALSE)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  sum_ij <- sum(choose2(tab))
  sum_i <- sum(choose2(rowSums(tab)))
  sum_j <- sum(choose2(colSums(tab)))
  n <- sum(tab)
  total <- choose2(n)
  expected <- sum_i * sum_j / total
  max_index <- (sum_i + sum_j) / 2
  denom <- max_index - expected
  if (denom == 0) return(1)
  (sum_ij - expected) / denom
}

mode_int <- function(x) {
  tx <- sort(table(x), decreasing = TRUE)
  as.integer(names(tx)[1])
}

permutations <- function(v) {
  if (length(v) == 1) return(matrix(v, nrow = 1))
  out <- lapply(seq_along(v), function(i) {
    rest <- v[-i]
    cbind(v[i], permutations(rest))
  })
  do.call(rbind, out)
}

best_permutation <- function(cost) {
  k <- nrow(cost)
  perms <- permutations(seq_len(k))
  scores <- apply(perms, 1, function(p) sum(cost[cbind(seq_len(k), p)]))
  perms[which.min(scores), ]
}

select_best_model <- function(grid_df, thr_appa = 0.70, thr_occ = 5, thr_prop = 0.05, fallback = TRUE) {
  if (nrow(grid_df) == 0) return(NULL)
  grid_df <- grid_df[!is.na(grid_df$BIC), , drop = FALSE]
  if (nrow(grid_df) == 0) return(NULL)
  grid_df$valid <- with(
    grid_df,
    (is.na(min_AvePP) | min_AvePP >= thr_appa) &
      (is.na(min_OCC) | min_OCC > thr_occ) &
      (is.na(min_prop) | min_prop >= thr_prop)
  )
  cand <- grid_df[grid_df$valid %in% TRUE, , drop = FALSE]
  if (nrow(cand) == 0) {
    if (!fallback) return(NULL)
    cand <- grid_df
  }
  cand[which.min(cand$BIC), , drop = FALSE]
}

model_summary_row <- function(b, seed_base, scenario_name, key, d, k, ms, dat_use, n_start) {
  fit <- ms$best_model
  q <- extract_quality2(fit, dat_use, id_col = "boot_id")
  data.frame(
    B = b,
    seed_base = seed_base,
    scenario = scenario_name,
    key = key,
    d = d,
    K = k,
    best_start = ms$best_start,
    n_start = n_start,
    n_id = q$n_id,
    K_observed = q$K,
    AIC = q$AIC,
    BIC = q$BIC,
    SABIC = q$SABIC,
    loglik = q$loglik,
    min_AvePP = if (all(is.na(q$appa))) NA_real_ else min(q$appa, na.rm = TRUE),
    mean_AvePP = if (all(is.na(q$appa))) NA_real_ else mean(q$appa, na.rm = TRUE),
    min_OCC = if (all(is.na(q$occ))) NA_real_ else min(q$occ, na.rm = TRUE),
    min_prop = if (all(is.na(q$prop))) NA_real_ else min(as.numeric(q$prop), na.rm = TRUE),
    max_prop = if (all(is.na(q$prop))) NA_real_ else max(as.numeric(q$prop), na.rm = TRUE),
    entropy = q$entropy,
    error = FALSE,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

model_error_row <- function(b, seed_base, scenario_name, key, d, k, err, n_start) {
  data.frame(
    B = b,
    seed_base = seed_base,
    scenario = scenario_name,
    key = key,
    d = d,
    K = k,
    best_start = NA_integer_,
    n_start = n_start,
    n_id = NA_integer_,
    K_observed = NA_integer_,
    AIC = NA_real_,
    BIC = NA_real_,
    SABIC = NA_real_,
    loglik = NA_real_,
    min_AvePP = NA_real_,
    mean_AvePP = NA_real_,
    min_OCC = NA_real_,
    min_prop = NA_real_,
    max_prop = NA_real_,
    entropy = NA_real_,
    error = TRUE,
    error_message = conditionMessage(err),
    stringsAsFactors = FALSE
  )
}

script_path <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
  else normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_path)
revise_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
code_dir <- normalizePath(file.path(revise_dir, ".."), winslash = "/", mustWork = TRUE)
root_dir <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = TRUE)

scenario_name <- parse_arg("scenario", "LAC_MIN_SF_PF_MEAN")
d_vec <- parse_int_vec("d-vec", "1")
k_vec <- parse_int_vec("k-vec", "5,6,7")
B_total <- as.integer(parse_arg("B", "500"))
b_start <- as.integer(parse_arg("b-start", "1"))
b_end <- as.integer(parse_arg("b-end", as.character(B_total)))
n_start <- as.integer(parse_arg("n-start", "5"))
seed <- as.integer(parse_arg("seed", "20260603"))
ref_k <- as.integer(parse_arg("ref-k", "6"))
light_env <- parse_arg("light-env", "")
quiet_fit <- parse_arg("quiet-fit", "true")
quiet_fit <- tolower(quiet_fit) %in% c("true", "t", "1", "yes", "y")
thr_appa <- as.numeric(parse_arg("thr-appa", "0.70"))
thr_occ <- as.numeric(parse_arg("thr-occ", "5"))
thr_prop <- as.numeric(parse_arg("thr-prop", "0.05"))
out_root <- file.path(revise_dir, "outputs", "bootstrap_stability")

run_args <- list(
  scenario_name = scenario_name,
  d_vec = d_vec,
  k_vec = k_vec,
  B_total = B_total,
  b_start = b_start,
  b_end = b_end,
  n_start = n_start,
  seed = seed,
  ref_k = ref_k,
  light_env = light_env,
  quiet_fit = quiet_fit,
  thr_appa = thr_appa,
  thr_occ = thr_occ,
  thr_prop = thr_prop
)

run_tag <- paste0(
  scenario_name,
  "_d", paste(d_vec, collapse = "-"),
  "_K", paste(k_vec, collapse = "-"),
  "_B", b_start, "-", b_end,
  "_seed", seed,
  "_nstart", n_start
)
out_dir <- file.path(out_root, run_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "bootstrap_stability.log")
sink(log_file, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("GBMTM bootstrap stability\n")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Root: ", root_dir, "\n", sep = "")
cat("Scenario: ", scenario_name, "\n", sep = "")
cat("d_vec: ", paste(d_vec, collapse = ","), "\n", sep = "")
cat("k_vec: ", paste(k_vec, collapse = ","), "\n", sep = "")
cat("B range: ", b_start, "-", b_end, " of ", B_total, "\n", sep = "")
cat("n_start: ", n_start, "\n", sep = "")
cat("Seed: ", seed, "\n", sep = "")
cat("quiet_fit: ", quiet_fit, "\n", sep = "")
cat("Output dir: ", out_dir, "\n\n", sep = "")

expected_r_version <- "4.4.2"
cat("R version in manuscript Methods: ", expected_r_version, "\n", sep = "")
cat("R version used for this run: ", R.version.string, "\n", sep = "")
if (as.character(getRversion()) != expected_r_version) {
  stop("This script must be run with R ", expected_r_version,
       ". Current version is ", as.character(getRversion()), call. = FALSE)
}

if (nzchar(light_env)) {
  env_file <- normalizePath(light_env, winslash = "/", mustWork = TRUE)
  cat("Loading lightweight bootstrap RData environment: ", env_file, "\n", sep = "")
} else {
  rdata_candidates <- list.files(code_dir, pattern = "10.*[.]RData$", full.names = TRUE)
  if (length(rdata_candidates) == 0) stop("No '*10*.RData' file found.", call. = FALSE)
  env_file <- rdata_candidates[which.max(file.info(rdata_candidates)$mtime)]
  cat("Loading RData environment: ", env_file, "\n", sep = "")
}
load_time <- system.time(load(env_file))
cat("RData loaded in seconds: ", unname(load_time[["elapsed"]]), "\n", sep = "")

scenario_name <- run_args$scenario_name
d_vec <- run_args$d_vec
k_vec <- run_args$k_vec
B_total <- run_args$B_total
b_start <- run_args$b_start
b_end <- run_args$b_end
n_start <- run_args$n_start
seed <- run_args$seed
ref_k <- run_args$ref_k
light_env <- run_args$light_env
quiet_fit <- run_args$quiet_fit
thr_appa <- run_args$thr_appa
thr_occ <- run_args$thr_occ
thr_prop <- run_args$thr_prop

script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)
revise_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
code_dir <- normalizePath(file.path(revise_dir, ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(revise_dir, "outputs", "bootstrap_stability", run_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

r_files <- list.files(code_dir, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
has_gbmtm_marker <- vapply(r_files, function(f) {
  x <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  any(grepl("^run_one_scenario\\s*<-\\s*function\\s*\\(", x))
}, logical(1))
source_script <- r_files[has_gbmtm_marker][1]
if (is.na(source_script)) stop("Could not find GBMTM modeling script.", call. = FALSE)
source_lines <- readLines(source_script, warn = FALSE)
stop_at <- grep("^scenarios\\s*<-\\s*tibble::tribble\\(", source_lines)
if (length(stop_at) == 0) stop("Could not find scenario definition in GBMTM modeling script.", call. = FALSE)
eval(parse(text = source_lines[seq_len(stop_at[[1]] - 1)]), envir = .GlobalEnv)
cat("Original helper functions loaded from: ", source_script, "\n", sep = "")

scenario_map <- tibble::tribble(
  ~scenario,                  ~lactate_col,    ~sf_col,    ~pf_col,
  "LAC_MIN_SF_PF_MIN",        "lactate_min",   "sf_min",   "pf_min",
  "LAC_MEAN_SF_PF_MIN",       "lactate_mean",  "sf_min",   "pf_min",
  "LAC_MAX_SF_PF_MIN",        "lactate_max",   "sf_min",   "pf_min",
  "LAC_MIN_SF_PF_MEAN",       "lactate_min",   "sf_mean",  "pf_mean",
  "LAC_MEAN_SF_PF_MEAN",      "lactate_mean",  "sf_mean",  "pf_mean",
  "LAC_MAX_SF_PF_MEAN",       "lactate_max",   "sf_mean",  "pf_mean"
)
sc <- scenario_map[scenario_map$scenario == scenario_name, , drop = FALSE]
if (nrow(sc) != 1) stop("Unknown scenario: ", scenario_name, call. = FALSE)
var_names <- c(sc$lactate_col, sc$sf_col, sc$pf_col, "ne_total_mcg_log")
need_cols <- c("stay_id_code", "day", var_names)

if (exists("scenario_list") && is.list(scenario_list) && scenario_name %in% names(scenario_list)) {
  dat_source <- as.data.frame(scenario_list[[scenario_name]])
  cat("Using imputed scenario data frame from scenario_list: ", scenario_name, "\n", sep = "")
} else {
  dat_source <- dat0
  cat("Using dat0 compatibility data frame for scenario: ", scenario_name, "\n", sep = "")
}
missing_need_cols <- setdiff(need_cols, names(dat_source))
if (length(missing_need_cols) > 0) {
  stop("Bootstrap input is missing required columns for ", scenario_name, ": ",
       paste(missing_need_cols, collapse = ", "), call. = FALSE)
}

dat_use <- dat_source %>%
  dplyr::select(dplyr::all_of(need_cols)) %>%
  dplyr::arrange(.data$stay_id_code, .data$day) %>%
  as.data.frame()

ids <- sort(unique(dat_use$stay_id_code))
n_id <- length(ids)
dat_by_id <- split(dat_use, dat_use$stay_id_code)

ref_group_path <- file.path(code_dir, paste0("GBMTM_", scenario_name, "_BEST_d1_K6_group.csv"))
ref_traj_path <- file.path(code_dir, paste0("GBMTM_", scenario_name, "_BEST_d1_K6_group_day_means.csv"))
if (!file.exists(ref_group_path)) stop("Reference group file not found: ", ref_group_path, call. = FALSE)
if (!file.exists(ref_traj_path)) stop("Reference trajectory file not found: ", ref_traj_path, call. = FALSE)
ref_group <- read.csv(ref_group_path, stringsAsFactors = FALSE)
ref_traj <- read.csv(ref_traj_path, stringsAsFactors = FALSE)
ref_assign_map <- setNames(as.integer(ref_group$group), as.character(ref_group$stay_id_code))
ref_group_levels <- sort(unique(ref_group$group))
if (length(ref_group_levels) != ref_k) stop("Reference group count is not ref_k.", call. = FALSE)

make_boot_data <- function(sampled_ids) {
  pieces <- vector("list", length(sampled_ids))
  for (i in seq_along(sampled_ids)) {
    one <- dat_by_id[[as.character(sampled_ids[i])]]
    one$orig_stay_id_code <- sampled_ids[i]
    one$boot_id <- i
    pieces[[i]] <- one[, c("boot_id", "orig_stay_id_code", "day", var_names)]
  }
  out <- dplyr::bind_rows(pieces)
  out[order(out$boot_id, out$day), , drop = FALSE]
}

group_means_for_model <- function(boot_dat, fit, group_col = "group") {
  id_assign <- data.frame(
    boot_id = sort(unique(boot_dat$boot_id)),
    group = as.integer(fit$assign)
  )
  if (nrow(id_assign) != length(fit$assign)) stop("Assignment length mismatch.", call. = FALSE)
  dat_g <- dplyr::left_join(boot_dat, id_assign, by = "boot_id")
  dat_g %>%
    dplyr::group_by(.data[[group_col]], .data$day) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      dplyr::across(dplyr::all_of(var_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = .data[[group_col]]) %>%
    dplyr::arrange(.data$group, .data$day)
}

trajectory_cost <- function(ref_traj, boot_traj) {
  ref_days <- sort(unique(ref_traj$day))
  scales <- vapply(var_names, function(v) {
    s <- stats::sd(ref_traj[[v]], na.rm = TRUE)
    if (is.na(s) || s == 0) 1 else s
  }, numeric(1))
  cost <- matrix(0, nrow = ref_k, ncol = ref_k)
  for (rg in seq_len(ref_k)) {
    for (bg in seq_len(ref_k)) {
      rr <- ref_traj[ref_traj$group == rg & ref_traj$day %in% ref_days, c("day", var_names)]
      bb <- boot_traj[boot_traj$group == bg & boot_traj$day %in% ref_days, c("day", var_names)]
      mm <- merge(rr, bb, by = "day", suffixes = c("_ref", "_boot"))
      if (nrow(mm) == 0) {
        cost[rg, bg] <- Inf
      } else {
        parts <- vapply(var_names, function(v) {
          mean(((mm[[paste0(v, "_ref")]] - mm[[paste0(v, "_boot")]]) / scales[[v]])^2, na.rm = TRUE)
        }, numeric(1))
        cost[rg, bg] <- sum(parts, na.rm = TRUE)
      }
    }
  }
  cost
}

alignment_map_for_k6 <- function(boot_traj) {
  if (!all(seq_len(ref_k) %in% boot_traj$group)) return(NULL)
  cost <- trajectory_cost(ref_traj, boot_traj)
  if (any(!is.finite(cost))) return(NULL)
  perm <- best_permutation(cost)
  map <- rep(NA_integer_, ref_k)
  for (rg in seq_len(ref_k)) map[perm[rg]] <- rg
  setNames(map, as.character(seq_len(ref_k)))
}

assignments_for_model <- function(boot_dat, fit) {
  data.frame(
    boot_id = sort(unique(boot_dat$boot_id)),
    orig_stay_id_code = as.integer(boot_dat$orig_stay_id_code[match(sort(unique(boot_dat$boot_id)), boot_dat$boot_id)]),
    group = as.integer(fit$assign)
  )
}

apply_alignment <- function(assign_df, boot_dat, fit) {
  boot_traj <- group_means_for_model(boot_dat, fit)
  map <- alignment_map_for_k6(boot_traj)
  if (is.null(map)) {
    assign_df$aligned_group <- NA_integer_
    return(list(assign = assign_df, trajectory = NULL, map = NULL))
  }
  assign_df$aligned_group <- as.integer(map[as.character(assign_df$group)])
  dat_g <- dplyr::left_join(
    boot_dat,
    assign_df[, c("boot_id", "aligned_group")],
    by = "boot_id"
  )
  traj <- dat_g %>%
    dplyr::group_by(.data$aligned_group, .data$day) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      dplyr::across(dplyr::all_of(var_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = .data$aligned_group) %>%
    dplyr::arrange(.data$group, .data$day)
  list(assign = assign_df, trajectory = traj, map = map)
}

fit_boot_grid <- function(b, boot_dat) {
  rows <- list()
  models <- list()
  for (d in d_vec) {
    for (k in k_vec) {
      key <- paste0("d", d, "_K", k)
      seed_base <- seed + b * 100000L + d * 1000L + k * 100L
      cat("  B=", b, " fitting ", key, " seed_base=", seed_base, " ... ", sep = "")
      one <- tryCatch(
        {
          if (quiet_fit) {
            ans <- NULL
            invisible(capture.output({
              ans <- fit_one_multistart(
                dat = boot_dat,
                varNames = var_names,
                id_col = "boot_id",
                time_col = "day",
                d = d,
                ng = k,
                n_start = n_start,
                seed_base = seed_base,
                scaling = 0,
                verbose = FALSE
              )
            }, file = NULL))
            ans
          } else {
            fit_one_multistart(
              dat = boot_dat,
              varNames = var_names,
              id_col = "boot_id",
              time_col = "day",
              d = d,
              ng = k,
              n_start = n_start,
              seed_base = seed_base,
              scaling = 0,
              verbose = FALSE
            )
          }
        },
        error = function(e) e
      )
      if (inherits(one, "error")) {
        rows[[key]] <- model_error_row(b, seed_base, scenario_name, key, d, k, one, n_start)
        cat("ERROR: ", conditionMessage(one), "\n", sep = "")
      } else {
        rows[[key]] <- model_summary_row(b, seed_base, scenario_name, key, d, k, one, boot_dat, n_start)
        models[[key]] <- one$best_model
        cat("BIC=", round(rows[[key]]$BIC, 2), ", min_prop=", round(rows[[key]]$min_prop, 4), "\n", sep = "")
      }
    }
  }
  list(summary = dplyr::bind_rows(rows), models = models)
}

iter_file <- file.path(out_dir, "bootstrap_iteration_summary.csv")
best_per_k_file <- file.path(out_dir, "bootstrap_best_per_K.csv")
full_grid_file <- file.path(out_dir, "bootstrap_full_grid_summary.csv")
props_file <- file.path(out_dir, "bootstrap_group_props_K6_aligned.csv")
traj_file <- file.path(out_dir, "bootstrap_trajectories_K6_aligned_long.csv")

co_counts <- matrix(0, nrow = n_id, ncol = n_id, dimnames = list(as.character(ids), as.character(ids)))
co_denoms <- matrix(0, nrow = n_id, ncol = n_id, dimnames = list(as.character(ids), as.character(ids)))

for (b in seq.int(b_start, b_end)) {
  cat("\nBootstrap ", b, "/", B_total, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  set.seed(seed + b)
  sampled_ids <- sample(ids, size = n_id, replace = TRUE)
  boot_dat <- make_boot_data(sampled_ids)
  fitted <- fit_boot_grid(b, boot_dat)
  grid_df <- fitted$summary
  append_csv(grid_df, full_grid_file)

  selected <- select_best_model(grid_df, thr_appa = thr_appa, thr_occ = thr_occ, thr_prop = thr_prop, fallback = FALSE)
  selected_key <- NA_character_
  ari_selected <- NA_real_
  if (!is.null(selected)) {
    selected_key <- selected$key[1]
    selected_fit <- fitted$models[[selected_key]]
    selected_assign <- assignments_for_model(boot_dat, selected_fit)
    selected_ref <- as.integer(ref_assign_map[as.character(selected_assign$orig_stay_id_code)])
    ari_selected <- adjusted_rand_index(selected_assign$group, selected_ref)
  }

  best_per_k_rows <- list()
  for (k in k_vec) {
    sub <- grid_df[grid_df$K == k, , drop = FALSE]
    sub_best <- select_best_model(sub, thr_appa = thr_appa, thr_occ = thr_occ, thr_prop = thr_prop, fallback = TRUE)
    if (is.null(sub_best)) next
    key <- sub_best$key[1]
    fit <- fitted$models[[key]]
    assign <- assignments_for_model(boot_dat, fit)
    ref <- as.integer(ref_assign_map[as.character(assign$orig_stay_id_code)])
    sub_best$ARI_vs_reference_K6 <- adjusted_rand_index(assign$group, ref)
    sub_best$selected_key <- selected_key
    best_per_k_rows[[as.character(k)]] <- sub_best
  }
  if (length(best_per_k_rows) > 0) append_csv(dplyr::bind_rows(best_per_k_rows), best_per_k_file)

  k6_best <- select_best_model(grid_df[grid_df$K == ref_k, , drop = FALSE], thr_appa = thr_appa, thr_occ = thr_occ, thr_prop = thr_prop, fallback = TRUE)
  ari_k6 <- NA_real_
  k6_key <- NA_character_
  if (!is.null(k6_best)) {
    k6_key <- k6_best$key[1]
    k6_fit <- fitted$models[[k6_key]]
    k6_assign <- assignments_for_model(boot_dat, k6_fit)
    k6_ref <- as.integer(ref_assign_map[as.character(k6_assign$orig_stay_id_code)])
    ari_k6 <- adjusted_rand_index(k6_assign$group, k6_ref)
    aligned <- apply_alignment(k6_assign, boot_dat, k6_fit)

    if (!all(is.na(aligned$assign$aligned_group))) {
      prop_tab <- prop.table(table(factor(aligned$assign$aligned_group, levels = seq_len(ref_k))))
      prop_df <- data.frame(
        B = b,
        group = seq_len(ref_k),
        prop = as.numeric(prop_tab),
        k6_key = k6_key
      )
      append_csv(prop_df, props_file)

      traj <- aligned$trajectory
      if (!is.null(traj)) {
        traj$B <- b
        traj$k6_key <- k6_key
        append_csv(traj[, c("B", "k6_key", "group", "day", "n_obs", var_names)], traj_file)
      }

      modal <- stats::aggregate(
        aligned_group ~ orig_stay_id_code,
        data = aligned$assign,
        FUN = mode_int
      )
      present_idx <- match(as.character(modal$orig_stay_id_code), as.character(ids))
      co_denoms[present_idx, present_idx] <- co_denoms[present_idx, present_idx] + 1
      for (g in seq_len(ref_k)) {
        idx <- present_idx[modal$aligned_group == g]
        if (length(idx) > 0) co_counts[idx, idx] <- co_counts[idx, idx] + 1
      }
    }
  }

  append_csv(data.frame(
    B = b,
    selected_key = selected_key,
    selected_d = if (is.null(selected)) NA_integer_ else selected$d[1],
    selected_K = if (is.null(selected)) NA_integer_ else selected$K[1],
    selected_BIC = if (is.null(selected)) NA_real_ else selected$BIC[1],
    selected_valid = if (is.null(selected)) FALSE else selected$valid[1],
    ARI_selected_vs_reference = ari_selected,
    K6_key = k6_key,
    ARI_K6_vs_reference = ari_k6,
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  ), iter_file)

  saveRDS(
    list(co_counts = co_counts, co_denoms = co_denoms, ids = ids),
    file.path(out_dir, "bootstrap_co_class_checkpoint.rds")
  )
}

iter_summary <- if (file.exists(iter_file)) read.csv(iter_file, stringsAsFactors = FALSE) else data.frame()
if (nrow(iter_summary) > 0) {
  selected_label <- ifelse(is.na(iter_summary$selected_K), "none_acceptable", as.character(iter_summary$selected_K))
  k_freq <- as.data.frame(table(selected_label), stringsAsFactors = FALSE)
  names(k_freq) <- c("selected_K", "n")
  k_freq$frequency <- k_freq$n / sum(k_freq$n)
  write.csv(k_freq, file.path(out_dir, "bootstrap_K_selection_frequency.csv"), row.names = FALSE)

  ari_summary <- data.frame(
    metric = c("ARI_selected_vs_reference", "ARI_K6_vs_reference"),
    n = c(sum(!is.na(iter_summary$ARI_selected_vs_reference)), sum(!is.na(iter_summary$ARI_K6_vs_reference))),
    mean = c(mean(iter_summary$ARI_selected_vs_reference, na.rm = TRUE), mean(iter_summary$ARI_K6_vs_reference, na.rm = TRUE)),
    median = c(stats::median(iter_summary$ARI_selected_vs_reference, na.rm = TRUE), stats::median(iter_summary$ARI_K6_vs_reference, na.rm = TRUE)),
    q025 = c(stats::quantile(iter_summary$ARI_selected_vs_reference, 0.025, na.rm = TRUE), stats::quantile(iter_summary$ARI_K6_vs_reference, 0.025, na.rm = TRUE)),
    q975 = c(stats::quantile(iter_summary$ARI_selected_vs_reference, 0.975, na.rm = TRUE), stats::quantile(iter_summary$ARI_K6_vs_reference, 0.975, na.rm = TRUE))
  )
  write.csv(ari_summary, file.path(out_dir, "bootstrap_ARI_summary.csv"), row.names = FALSE)
}

if (file.exists(props_file)) {
  props <- read.csv(props_file, stringsAsFactors = FALSE)
  prop_intervals <- props %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data$prop, na.rm = TRUE),
      q025 = stats::quantile(.data$prop, 0.025, na.rm = TRUE),
      q50 = stats::quantile(.data$prop, 0.50, na.rm = TRUE),
      q975 = stats::quantile(.data$prop, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(prop_intervals, file.path(out_dir, "bootstrap_group_prop_intervals_K6.csv"), row.names = FALSE)
}

if (file.exists(traj_file)) {
  traj <- read.csv(traj_file, stringsAsFactors = FALSE)
  traj_long <- tidyr::pivot_longer(
    traj,
    cols = dplyr::all_of(var_names),
    names_to = "variable",
    values_to = "value"
  )
  traj_intervals <- traj_long %>%
    dplyr::group_by(.data$group, .data$day, .data$variable) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data$value, na.rm = TRUE),
      q025 = stats::quantile(.data$value, 0.025, na.rm = TRUE),
      q50 = stats::quantile(.data$value, 0.50, na.rm = TRUE),
      q975 = stats::quantile(.data$value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(traj_intervals, file.path(out_dir, "bootstrap_trajectory_intervals_K6.csv"), row.names = FALSE)
}

co_prob <- co_counts / co_denoms
co_prob[co_denoms == 0] <- NA_real_
write.csv(co_prob, file.path(out_dir, "bootstrap_co_class_probability_matrix.csv"))
saveRDS(
  list(
    scenario = scenario_name,
    d_vec = d_vec,
    k_vec = k_vec,
    B_total = B_total,
    b_start = b_start,
    b_end = b_end,
    n_start = n_start,
    seed = seed,
    ref_k = ref_k,
    var_names = var_names,
    co_counts = co_counts,
    co_denoms = co_denoms,
    co_prob = co_prob
  ),
  file.path(out_dir, "bootstrap_stability_summary.rds")
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  if (exists("k_freq")) {
    p <- ggplot2::ggplot(k_freq, ggplot2::aes(x = factor(.data$selected_K), y = .data$frequency)) +
      ggplot2::geom_col(fill = "#3B6FB6", width = 0.7) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = "Selected K", y = "Bootstrap frequency") +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_K_selection_frequency.png"), p, width = 6, height = 4, dpi = 300)
  }
  if (nrow(iter_summary) > 0 && "ARI_K6_vs_reference" %in% names(iter_summary)) {
    iter_summary$ARI_K6_vs_reference <- suppressWarnings(as.numeric(iter_summary$ARI_K6_vs_reference))
  }
  if (nrow(iter_summary) > 0 && "ARI_K6_vs_reference" %in% names(iter_summary) && any(is.finite(iter_summary$ARI_K6_vs_reference))) {
    p <- ggplot2::ggplot(iter_summary[is.finite(iter_summary$ARI_K6_vs_reference), , drop = FALSE], ggplot2::aes(x = .data$ARI_K6_vs_reference)) +
      ggplot2::geom_histogram(bins = 25, fill = "#4C956C", color = "white") +
      ggplot2::labs(x = "ARI vs reference K=6", y = "Bootstrap count") +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_ARI_K6_distribution.png"), p, width = 6, height = 4, dpi = 300)
  }
  if (exists("prop_intervals")) {
    p <- ggplot2::ggplot(prop_intervals, ggplot2::aes(x = factor(.data$group), y = .data$mean, ymin = .data$q025, ymax = .data$q975)) +
      ggplot2::geom_pointrange(color = "#7B2CBF") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = "Reference group", y = "Bootstrap class proportion") +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_group_prop_intervals_K6.png"), p, width = 6, height = 4, dpi = 300)
  }
  if (exists("traj_intervals")) {
    p <- ggplot2::ggplot(traj_intervals, ggplot2::aes(x = .data$day, y = .data$q50, ymin = .data$q025, ymax = .data$q975, group = .data$group)) +
      ggplot2::geom_ribbon(fill = "#9ECAE1", alpha = 0.35) +
      ggplot2::geom_line(color = "#08519C", linewidth = 0.4) +
      ggplot2::facet_grid(variable ~ group, scales = "free_y") +
      ggplot2::labs(x = "Day", y = "Bootstrap trajectory") +
      ggplot2::theme_minimal(base_size = 9)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_trajectory_intervals_K6.png"), p, width = 12, height = 8, dpi = 300)
  }
}

cat("\nCompleted: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
