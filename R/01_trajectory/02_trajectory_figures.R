
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(patchwork)
})

if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}
suppressPackageStartupMessages(library(ggridges))

stopifnot(exists("dat0"))
stopifnot(exists("all_results"))
stopifnot(is.list(all_results))
stopifnot(length(names(all_results)) > 0)

make_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

theme_pub <- function() {
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text  = element_text(size = 8),
    plot.title   = element_text(size = 8),
    plot.subtitle = element_text(size = 8),
    plot.caption  = element_text(size = 8)
  )
}

stack_AB <- function(plotA, plotB) {
  (plotA + theme_pub() +
     labs(tag = "A") +
     theme(plot.tag = element_text(size = 10, face = "bold"),
           plot.tag.position = c(0, 1))) /
    (plotB + theme_pub() +
       labs(tag = "B") +
       theme(plot.tag = element_text(size = 10, face = "bold"),
             plot.tag.position = c(0, 1))) +
    plot_layout(heights = c(1.1, 1))
}

var_labels_plain <- c(
  lactate_min      = "Lactate (Min)",
  lactate_mean     = "Lactate (Mean)",
  lactate_max      = "Lactate (Max)",
  ne_total_mcg_log = "Norepinephrine (log1p mcg)",
  pf_min           = "PaO2/FiO2 (Min)",
  pf_mean          = "PaO2/FiO2 (Mean)",
  sf_min           = "SpO2/FiO2 (Min)",
  sf_mean          = "SpO2/FiO2 (Mean)"
)

facet_labeller <- function(mapping) {
  function(x) {
    out <- mapping[x]
    out[is.na(out)] <- x[is.na(out)]
    unname(out)
  }
}

plot_obs_mean_se_from_data <- function(data_with_group,
                                       vars,
                                       group_col = "group",
                                       time_col = "day",
                                       x_label = "Day",
                                       y_label = "Value",
                                       facet_titles_map = NULL,
                                       conf_z = 1.96,
                                       show_points = TRUE,
                                       smooth = FALSE,
                                       loess_span = 0.75,
                                       ribbon_alpha = 0.20) {

  stopifnot(all(c(group_col, time_col) %in% names(data_with_group)))
  stopifnot(all(vars %in% names(data_with_group)))

  long_df <- data_with_group %>%
    mutate(
      Group = factor(.data[[group_col]]),
      TimePt = as.numeric(.data[[time_col]])
    ) %>%
    pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value")

  sum_df <- long_df %>%
    group_by(Group, TimePt, Variable) %>%
    summarise(
      N    = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      Sd   = sd(Value, na.rm = TRUE),
      Se   = Sd / sqrt(pmax(N, 1)),
      Lwr  = Mean - conf_z * Se,
      Upr  = Mean + conf_z * Se,
      .groups = "drop"
    )

  lab_fun <- if (!is.null(facet_titles_map)) {
    as_labeller(facet_labeller(facet_titles_map))
  } else {
    label_value
  }

  p <- ggplot(sum_df, aes(x = TimePt, y = Mean, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = Lwr, ymax = Upr), alpha = ribbon_alpha, color = NA, show.legend = FALSE) +
    {if (smooth) geom_smooth(method = "loess", span = loess_span, se = FALSE, linewidth = 1.05)
      else       geom_line(linewidth = 1.05)} +
    {if (show_points) geom_point(size = 2)} +
    facet_wrap(~ Variable, nrow = 1, scales = "free_y", labeller = lab_fun) +
    scale_x_continuous(breaks = sort(unique(sum_df$TimePt))) +
    labs(x = x_label, y = y_label, color = "Group", fill = "Group") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  list(plot = p, data = sum_df)
}

plot_ridge_simple <- function(data_with_group,
                              vars,
                              group_col = "group",
                              facet_titles_map = NULL,
                              ridge_alpha = 0.75,
                              ridge_scale = 1.0) {

  stopifnot(all(c(group_col, vars) %in% names(data_with_group)))

  ridge_long <- data_with_group %>%
    mutate(Group = factor(.data[[group_col]])) %>%
    pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>%
    filter(is.finite(Value))

  lab_fun <- if (!is.null(facet_titles_map)) {
    as_labeller(facet_labeller(facet_titles_map))
  } else {
    label_value
  }

  ggplot(ridge_long, aes(x = Value, y = Group, fill = Group)) +
    ggridges::geom_density_ridges(alpha = ridge_alpha, scale = ridge_scale, color = NA) +
    facet_wrap(~ Variable, nrow = 1, scales = "free_x", labeller = lab_fun) +
    labs(x = "Value", y = "Group", fill = "Group") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

for (scenario_name in names(all_results)) {

  res <- all_results[[scenario_name]]

  if (is.null(res$best_key) || is.null(res$varNames)) {
    warning("Skip scenario (missing best_key/varNames): ", scenario_name)
    next
  }

  scenario_safe <- make_safe_name(scenario_name)
  best_key <- res$best_key
  vars_use <- as.character(res$varNames)

  group_file <- paste0("GBMTM_", scenario_name, "_BEST_", best_key, "_group.csv")
  if (!file.exists(group_file)) {
    warning("Missing group file: ", group_file, " (skip scenario: ", scenario_name, ")")
    next
  }

  group_df <- read.csv(group_file) %>%
    mutate(
      stay_id_code = as.integer(stay_id_code),
      group = as.integer(group)
    )

  need_cols <- unique(c("stay_id_code", "day", vars_use))
  miss <- setdiff(need_cols, names(dat0))
  if (length(miss) > 0) stop("dat0 missing columns: ", paste(miss, collapse = ", "))

  dat_with_group <- dat0 %>%
    select(all_of(need_cols)) %>%
    mutate(
      stay_id_code = as.integer(stay_id_code),
      day = as.integer(day)
    ) %>%
    left_join(group_df, by = "stay_id_code")

  if (any(is.na(dat_with_group$group))) {
    stop("Join produced NA group for scenario: ", scenario_name,
         ". Check stay_id_code alignment between dat0 and group file.")
  }

  outA <- plot_obs_mean_se_from_data(
    data_with_group = dat_with_group,
    vars = vars_use,
    facet_titles_map = var_labels_plain,
    show_points = TRUE,
    smooth = FALSE,
    ribbon_alpha = 0.20
  )

  pB <- plot_ridge_simple(
    data_with_group = dat_with_group,
    vars = vars_use,
    facet_titles_map = var_labels_plain,
    ridge_alpha = 0.80,
    ridge_scale = 1.0
  )

  figAB <- stack_AB(outA$plot, pB)

  assign(paste0("dat_with_group_", scenario_safe), dat_with_group, envir = .GlobalEnv)

  assign(paste0("dfObsSe_", scenario_safe), outA$data, envir = .GlobalEnv)
  assign(paste0("pA_obsse_", scenario_safe), outA$plot, envir = .GlobalEnv)

  assign(paste0("pB_ridge_", scenario_safe), pB, envir = .GlobalEnv)

  assign(paste0("figAB_", scenario_safe), figAB, envir = .GlobalEnv)

  file_out <- paste0("FIG_", scenario_safe, "_AB.png")
  ggsave(file_out, figAB, width = 11, height = 9, dpi = 300)
  assign(paste0("fileFigAB_", scenario_safe), file_out, envir = .GlobalEnv)

  file_A <- paste0("FIG_", scenario_safe, "_A_OBS_SE.png")
  file_B <- paste0("FIG_", scenario_safe, "_B_RIDGE.png")
  ggsave(file_A, outA$plot, width = 11, height = 4.2, dpi = 300)
  ggsave(file_B, pB,        width = 11, height = 4.2, dpi = 300)
  assign(paste0("fileA_", scenario_safe), file_A, envir = .GlobalEnv)
  assign(paste0("fileB_", scenario_safe), file_B, envir = .GlobalEnv)

  cat("OK: ", scenario_name, "\n",
      "  Output object: dat_with_group_", scenario_safe, "\n",
      "           pA_obsse_", scenario_safe, "\n",
      "           pB_ridge_", scenario_safe, "\n",
      "           figAB_", scenario_safe, "\n",
      "  Saved:   ", file_out, "\n", sep = "")
}













class(all_results$LAC_MIN_SF_PF_MEAN$best_model)


