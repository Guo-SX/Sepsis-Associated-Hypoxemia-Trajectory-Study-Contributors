
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
})

stopifnot(exists("all_results"))
stopifnot(is.list(all_results))
stopifnot(length(names(all_results)) > 0)

make_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

class_pct6 <- function(assign, K) {
  out <- rep(NA_real_, 6)
  if (!is.null(assign) && length(assign) > 0 && !is.na(K) && K >= 1) {
    tb <- prop.table(table(assign)) * 100
    for (g in seq_len(min(K, 6))) {
      if (as.character(g) %in% names(tb)) out[g] <- as.numeric(tb[[as.character(g)]])
    }
  }
  setNames(out, paste0("%class", 1:6))
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x
}

df_to_latex_tabular <- function(df,
                                caption = NULL,
                                label = NULL,
                                align = NULL,
                                digits = 2,
                                booktabs = TRUE,
                                sanitize = TRUE) {

  if (is.null(align)) align <- paste0(rep("l", ncol(df)), collapse = "")

  df2 <- df
  for (j in seq_len(ncol(df2))) {
    if (is.numeric(df2[[j]])) {
      df2[[j]] <- format(round(df2[[j]], digits), nsmall = digits, trim = TRUE)
    }
  }

  if (sanitize) {
    df2[] <- lapply(df2, latex_escape)
    names(df2) <- latex_escape(names(df2))
  }

  header <- paste(names(df2), collapse = " & ")
  body <- apply(df2, 1, function(r) paste(r, collapse = " & "))

  lines <- c()
  if (!is.null(caption) || !is.null(label)) lines <- c(lines, "\\begin{table}[htbp]", "\\centering")

  if (booktabs) {
    lines <- c(lines,
               sprintf("\\begin{tabular}{%s}", align),
               "\\toprule",
               paste0(header, " \\\\"),
               "\\midrule",
               paste0(body, " \\\\"),
               "\\bottomrule",
               "\\end{tabular}")
  } else {
    lines <- c(lines,
               sprintf("\\begin{tabular}{%s}", align),
               "\\hline",
               paste0(header, " \\\\ \\hline"),
               paste0(body, " \\\\"),
               "\\hline",
               "\\end{tabular}")
  }

  if (!is.null(caption)) lines <- c(lines, sprintf("\\caption{%s}", latex_escape(caption)))
  if (!is.null(label))   lines <- c(lines, sprintf("\\label{%s}", latex_escape(label)))
  if (!is.null(caption) || !is.null(label)) lines <- c(lines, "\\end{table}")

  paste(lines, collapse = "\n")
}

make_perf_table_full_grid <- function(res, scenario_name) {
  stopifnot(!is.null(res$model_choice))
  stopifnot(!is.null(res$res$grid))

  mc <- res$model_choice %>%
    mutate(
      Scenario = scenario_name,
      d = as.integer(d),
      K = as.integer(ng),
      BIC = as.numeric(BIC_use),
      AIC = as.numeric(AIC_use)
    )

  grid <- res$res$grid

  pct_df <- purrr::map_dfr(names(grid), function(key) {
    ms <- grid[[key]]
    K_here <- NA_integer_
    K_here <- suppressWarnings(as.integer(gsub(".*_K", "", key)))
    assign_vec <- ms$best_model$assign
    pct <- class_pct6(assign_vec, K_here)
    tibble(
      key = key,
      `%class1` = pct["%class1"],
      `%class2` = pct["%class2"],
      `%class3` = pct["%class3"],
      `%class4` = pct["%class4"],
      `%class5` = pct["%class5"],
      `%class6` = pct["%class6"]
    )
  })

  out <- mc %>%
    left_join(pct_df, by = "key") %>%
    select(
      Scenario, key, d, K,
      BIC, AIC, SABIC,
      min_AvePP, mean_AvePP, min_OCC, min_prop, max_prop, entropy,
      `%class1`, `%class2`, `%class3`, `%class4`, `%class5`, `%class6`,
      best_start
    ) %>%
    arrange(d, K)

  out
}

all_perf_tables <- list()

for (scenario_name in names(all_results)) {

  res <- all_results[[scenario_name]]
  if (is.null(res$model_choice) || is.null(res$res$grid)) {
    warning("Skip (missing model_choice or res$grid): ", scenario_name)
    next
  }

  sc_safe <- make_safe_name(scenario_name)

  perf_df <- make_perf_table_full_grid(res, scenario_name)

  all_perf_tables[[scenario_name]] <- perf_df
  assign(paste0("perf_all_tbl_", sc_safe), perf_df, envir = .GlobalEnv)

  perf_tex <- df_to_latex_tabular(
    perf_df,
    caption = paste0("GBMTM performance across all candidate models (d=1..3, K=1..6): ", scenario_name),
    label   = paste0("tab:perf_all_", sc_safe),
    align   = paste0("l", paste(rep("c", ncol(perf_df)-1), collapse = "")),
    digits  = 2,
    booktabs = TRUE
  )

  assign(paste0("perf_all_tex_", sc_safe), perf_tex, envir = .GlobalEnv)

  fn <- paste0("TEX_PERF_ALL_", sc_safe, ".tex")
  writeLines(perf_tex, con = fn)

  cat("Wrote: ", fn, " (rows=", nrow(perf_df), ")\n", sep = "")
}

perf_all_combined <- bind_rows(all_perf_tables)

assign("perf_all_tbl_ALL", perf_all_combined, envir = .GlobalEnv)

perf_all_combined_tex <- df_to_latex_tabular(
  perf_all_combined,
  caption = "GBMTM performance across all scenarios and all candidate models (d=1..3, K=1..6).",
  label   = "tab:perf_all_scenarios",
  align   = paste0("l", paste(rep("c", ncol(perf_all_combined)-1), collapse = "")),
  digits  = 2,
  booktabs = TRUE
)

assign("perf_all_tex_ALL", perf_all_combined_tex, envir = .GlobalEnv)
writeLines(perf_all_combined_tex, con = "TEX_PERF_ALL_SCENARIOS.tex")

cat("Wrote combined: TEX_PERF_ALL_SCENARIOS.tex (rows=", nrow(perf_all_combined), ")\n", sep = "")










for_merge_baseline <- baseline3 %>% select(-stay_id,-subject_id,-hadm_id,-intime,-outtime,-pf_ratio_min,-sf_ratio_min,-pao2_n)


for_baseline_table <- merge(dat_with_group_LAC_MIN_SF_PF_MEAN,for_merge_baseline,merge = "stay_id_code",all.x = TRUE)

colnames(for_baseline_table)


suppressPackageStartupMessages({
  library(dplyr)
  library(gtsummary)
  library(gt)
  library(tibble)
})

stopifnot(exists("for_baseline_table"))

df0 <- for_baseline_table %>%
  filter(!is.na(stay_id_code))

cat("[baseline] after drop NA stay_id_code:", nrow(df0), "rows\n")

baseline_vars <- c(
  "stay_id_code", "group",
  "gender", "age_at_admit",
  "copd", "chf", "ckd", "liver_disease", "diabetes", "malignancy",
  "icu_los_hours", "hosp_los_days",
  "death_28d", "death_90d", "icu_mortality", "hosp_mortality"
)

miss_b <- setdiff(baseline_vars, names(df0))
if (length(miss_b) > 0) stop("for_baseline_table 缺少基线表需要列：", paste(miss_b, collapse = ", "))

binary_vars <- c(
  "copd", "chf", "ckd", "liver_disease", "diabetes", "malignancy",
  "death_28d", "death_90d", "icu_mortality", "hosp_mortality"
)

table1_df <- df0 %>%
  select(all_of(baseline_vars)) %>%
  group_by(stay_id_code) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    group  = factor(group),
    gender = factor(gender),
    across(all_of(binary_vars), ~ dplyr::case_when(
      . %in% c(0, "0") ~ "NO",
      . %in% c(1, "1") ~ "YES",
      TRUE ~ NA_character_
    )),
    across(all_of(binary_vars), ~ factor(., levels = c("NO", "YES")))
  )

fisher_or_chisq <- function(data, variable, by, ...) {
  x <- data[[variable]]
  g <- data[[by]]

  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]

  if (length(x) == 0 || length(unique(x)) < 2 || length(unique(g)) < 2) {
    return(tibble(p.value = NA_real_))
  }

  tab <- table(x, g)

  p <- tryCatch({
    fisher.test(tab)$p.value
  }, error = function(e) {
    tryCatch({
      suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 20000)$p.value)
    }, error = function(e2) NA_real_)
  })

  tibble(p.value = p)
}

kruskal_safe <- function(data, variable, by, ...) {
  x <- data[[variable]]
  g <- data[[by]]

  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]

  if (length(x) == 0 || length(unique(g)) < 2 || length(unique(x)) < 2) {
    return(tibble(p.value = NA_real_))
  }

  p <- tryCatch({
    kruskal.test(x ~ g)$p.value
  }, error = function(e) NA_real_)

  tibble(p.value = p)
}

tbl1 <- table1_df %>%
  select(-stay_id_code) %>%
  tbl_summary(
    by = group,
    type = list(all_of(binary_vars) ~ "categorical"),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous()  ~ kruskal_safe,
      all_categorical() ~ fisher_or_chisq
    )
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels() %>%
  modify_fmt_fun(
    p.value ~ function(x) ifelse(is.na(x), "NA", x)
  )

gt_tbl  <- as_gt(tbl1)
tex_code <- gt::as_latex(gt_tbl)

out_tex <- "Table1_baseline_by_group.tex"
writeLines(tex_code, out_tex)

cat("Saved LaTeX to: ", normalizePath(out_tex), "\n", sep = "")
cat("\n--- LaTeX (copy to Overleaf) ---\n")
cat(paste(tex_code, collapse = "\n"))
cat("\n--- end ---\n")
