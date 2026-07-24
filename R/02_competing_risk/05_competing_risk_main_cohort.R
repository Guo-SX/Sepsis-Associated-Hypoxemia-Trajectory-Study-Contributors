
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(cmprsk)
  library(patchwork)
  library(forestploter)
  library(grid)
  library(scales)
})

cfg <- list(
  id_stay = "stay_id",
  id_code = "stay_id_code",
  group_col = "group",
  day_col   = "day",

  death_flag = "hosp_mortality",
  # Must be the terminal hospital-event time on the same scale as
  # hosp_mortality (death time for hospital deaths; live-discharge time
  # for survivors), expressed in hours from the chosen analysis origin.
  time_hours_col = "hospital_event_time_hours",

  age_col = "age_at_admit",
  cbc_cols = c("wbc","rbc","hb","hct","plt"),
  comorb_cols = c("copd","chf","liver_disease","diabetes"),

  ref_group = "1",

  resp_cols = c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"),
  ribbon_alpha = 0.18,

  forest_xlim_default = c(0, 4),

  scenario_label_fun = function(s) trimws(gsub("_+", " ", s)),

  tag_size = 12
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

stop_has_cols <- function(df, cols, dfname = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(dfname, " 缺少列: ", paste(miss, collapse = ", "))
}

fmt_p <- function(p, is_ref = FALSE) {
  if (is_ref) return("")
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_hrci <- function(hr, l, u, is_ref = FALSE, digits = 3) {
  if (is_ref) return("Reference")
  if (any(!is.finite(c(hr, l, u))) || any(is.na(c(hr, l, u)))) return("")
  sprintf(paste0("%.",digits,"f (%.",digits,"f–%.",digits,"f)"), hr, l, u)
}

make_group_df_from_long <- function(dat_long,
                                    id_code = cfg$id_code,
                                    group_col = cfg$group_col) {
  stop_has_cols(dat_long, c(id_code, group_col), "dat_long")

  dat_long %>%
    select(.data[[id_code]], .data[[group_col]]) %>%
    filter(!is.na(.data[[id_code]]), !is.na(.data[[group_col]])) %>%
    distinct() %>%
    group_by(.data[[id_code]]) %>%
    summarise(group = first(as.character(.data[[group_col]])), .groups = "drop")
}

attach_stay_id <- function(group_df_code,
                           imp_back = NULL,
                           id_code = cfg$id_code,
                           id_stay = cfg$id_stay) {
  if (id_stay %in% names(group_df_code)) return(group_df_code)

  if (is.null(imp_back)) {
    stop("group_df_code 没有 stay_id，同时你没提供 imp_back 用于 stay_id_code -> stay_id 对应。")
  }
  stop_has_cols(imp_back, c(id_code, id_stay), "imp_back")

  group_df_code %>%
    left_join(
      imp_back %>% select(.data[[id_code]], .data[[id_stay]]) %>% distinct(),
      by = setNames(id_code, id_code)
    )
}

build_fg_data <- function(baseline,
                          group_df_stay,
                          cfg = cfg) {
  if (is.null(cfg$time_hours_col) ||
      identical(cfg$time_hours_col, "icu_los_hours")) {
    stop(
      "Fine-Gray analysis requires a verified hospital terminal-event time. ",
      "Do not pair hosp_mortality with icu_los_hours."
    )
  }
  stop_has_cols(
    baseline,
    c(cfg$id_stay, cfg$death_flag, cfg$age_col, cfg$time_hours_col),
    "baseline"
  )
  stop_has_cols(group_df_stay, c(cfg$id_stay, "group"), "group_df_stay")

  dat <- baseline %>%
    inner_join(group_df_stay %>% select(.data[[cfg$id_stay]], group),
               by = setNames(cfg$id_stay, cfg$id_stay)) %>%
    mutate(
      group = factor(as.character(group)),
      death = as.integer(.data[[cfg$death_flag]])
    )

  dat <- dat %>% mutate(ftime = as.numeric(.data[[cfg$time_hours_col]]))

  dat %>%
    mutate(
      fstatus = case_when(
        death == 1 ~ 1L,
        death == 0 ~ 2L,
        TRUE ~ 0L
      )
    ) %>%
    filter(is.finite(ftime), ftime > 0, fstatus %in% c(0,1,2)) %>%
    droplevels()
}

fit_fg_models <- function(reg_df, cfg = cfg) {
  reg_df <- reg_df %>%
    mutate(
      group = factor(group),
      age = as.numeric(.data[[cfg$age_col]])
    )

  cbc_ok    <- intersect(cfg$cbc_cols, names(reg_df))
  comorb_ok <- intersect(cfg$comorb_cols, names(reg_df))

  for (v in cbc_ok)    reg_df[[v]] <- as.numeric(reg_df[[v]])
  for (v in comorb_ok) reg_df[[v]] <- as.numeric(reg_df[[v]])

  mkX <- function(df, rhs_terms) {
    ff <- as.formula(paste0("~ ", paste(rhs_terms, collapse = " + ")))
    mm <- model.matrix(ff, data = df)
    if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    mm
  }

  clean_X <- function(X) {
    if (is.null(X) || ncol(X) == 0) return(X)
    X <- X[, colSums(!is.finite(X)) == 0, drop = FALSE]
    if (ncol(X) == 0) return(X)
    v <- apply(X, 2, var)
    X <- X[, is.finite(v) & v > 0, drop = FALSE]
    if (ncol(X) == 0) return(X)
    dup <- duplicated(as.data.frame(t(X)))
    if (any(dup)) X <- X[, !dup, drop = FALSE]
    X
  }

  safe_crr <- function(ftime, fstatus, X, failcode = 1, cencode = 0) {
    X <- clean_X(X)
    if (is.null(X) || ncol(X) == 0) {
      X <- matrix(0, nrow = length(ftime), ncol = 1)
      colnames(X) <- "dummy0"
    }

    fit_try <- function(X_use) cmprsk::crr(ftime, fstatus, cov1 = X_use, failcode = failcode, cencode = cencode)

    out <- try(fit_try(X), silent = TRUE)
    if (!inherits(out, "try-error")) return(list(fit = out, X = X))

    X2 <- X
    while (ncol(X2) > 1) {
      X2 <- X2[, -ncol(X2), drop = FALSE]
      out2 <- try(fit_try(X2), silent = TRUE)
      if (!inherits(out2, "try-error")) return(list(fit = out2, X = X2))
    }

    stop("crr 拟合失败（可能分离/某组无事件）。原始错误：\n", out)
  }

  vars_m1 <- c("ftime","fstatus","group")
  vars_m2 <- unique(c(vars_m1, "age", cbc_ok))
  vars_m3 <- unique(c(vars_m2, comorb_ok))

  d1 <- reg_df %>% filter(complete.cases(across(all_of(vars_m1))))
  d2 <- reg_df %>% filter(complete.cases(across(all_of(vars_m2))))
  d3 <- reg_df %>% filter(complete.cases(across(all_of(vars_m3))))

  if (cfg$ref_group %in% levels(d1$group)) d1$group <- relevel(d1$group, ref = cfg$ref_group)
  if (cfg$ref_group %in% levels(d2$group)) d2$group <- relevel(d2$group, ref = cfg$ref_group)
  if (cfg$ref_group %in% levels(d3$group)) d3$group <- relevel(d3$group, ref = cfg$ref_group)

  X1 <- mkX(d1, c("group"))
  X2 <- mkX(d2, c("group", "age", cbc_ok))
  X3 <- mkX(d3, c("group", "age", cbc_ok, comorb_ok))

  f1 <- safe_crr(d1$ftime, d1$fstatus, X1)
  f2 <- safe_crr(d2$ftime, d2$fstatus, X2)
  f3 <- safe_crr(d3$ftime, d3$fstatus, X3)

  list(
    Model1 = f1$fit, data_m1 = d1,
    Model2 = f2$fit, data_m2 = d2,
    Model3 = f3$fit, data_m3 = d3,
    cbc_used = cbc_ok, comorb_used = comorb_ok
  )
}

extract_group_from_crr <- function(fit, ref_group, group_levels) {
  sm <- summary(fit)
  coef_tab <- sm$coef
  if (is.null(coef_tab) || nrow(coef_tab) == 0) {
    return(tibble(group = factor(group_levels, levels = group_levels),
                  HR = NA_real_, L = NA_real_, U = NA_real_, P = NA_real_))
  }

  rn <- rownames(coef_tab)
  b  <- as.numeric(coef_tab[, 1])
  se <- as.numeric(coef_tab[, 3])
  p  <- suppressWarnings(as.numeric(coef_tab[, ncol(coef_tab)]))
  names(b) <- names(se) <- names(p) <- rn

  out <- tibble(group = factor(group_levels, levels = group_levels),
                HR = NA_real_, L = NA_real_, U = NA_real_, P = NA_real_)

  out[out$group == ref_group, c("HR","L","U")] <- 1
  out[out$group == ref_group, "P"] <- NA_real_

  for (g in group_levels) {
    if (g == ref_group) next
    key <- paste0("group", g)
    if (!key %in% rn) next
    bi  <- b[key]; sei <- se[key]
    out[out$group == g, "HR"] <- exp(bi)
    out[out$group == g, "L"]  <- exp(bi - 1.96 * sei)
    out[out$group == g, "U"]  <- exp(bi + 1.96 * sei)
    out[out$group == g, "P"]  <- p[key]
  }
  out
}

get_gray_p <- function(ci_obj) {
  if (is.null(ci_obj$Tests)) return(NA_real_)
  tt <- ci_obj$Tests
  pcol <- intersect(c("pv","p","pvalue","P","Pvalue"), colnames(tt))
  if (length(pcol) >= 1) return(as.numeric(tt[1, pcol[1]]))
  suppressWarnings(as.numeric(tt[1, ncol(tt)]))
}

tidy_cuminc_ci <- function(ci_obj, event_code = 1) {
  nms <- names(ci_obj)
  keep <- grepl(paste0(" ", event_code, "$"), nms)
  nms  <- nms[keep]
  if (length(nms) == 0) return(tibble(group=factor(), time=numeric(), est=numeric(), lwr=numeric(), upr=numeric()))

  out <- lapply(nms, function(nm) {
    x <- ci_obj[[nm]]
    grp <- sub(paste0(" ", event_code, "$"), "", nm)
    est <- as.numeric(x$est)
    tt  <- as.numeric(x$time)
    vv  <- suppressWarnings(as.numeric(x$var))
    se  <- if (!is.null(vv) && length(vv) == length(est)) sqrt(pmax(vv, 0)) else rep(NA_real_, length(est))
    lwr <- pmax(est - 1.96 * se, 0)
    upr <- pmin(est + 1.96 * se, 1)
    tibble(group = grp, time = tt, est = est, lwr = lwr, upr = upr)
  })

  df <- bind_rows(out)
  df$group <- factor(df$group, levels = unique(df$group))
  df
}

make_plot_cif_death <- function(reg_df, tag_short, cfg = cfg) {
  k <- nlevels(reg_df$group)
  pal <- cfg$resp_cols[seq_len(min(k, length(cfg$resp_cols)))]

  ci_obj <- cmprsk::cuminc(reg_df$ftime, reg_df$fstatus, reg_df$group, cencode = 0)
  cif_df <- tidy_cuminc_ci(ci_obj, event_code = 1)
  gray_p <- get_gray_p(ci_obj)
  subtxt <- if (!is.na(gray_p)) sprintf("Gray p %.3g", gray_p) else NULL

  p <- ggplot(cif_df, aes(time, est, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = cfg$ribbon_alpha, linewidth = 0) +
    geom_step(linewidth = 1.1) +
    scale_color_manual(values = pal, labels = paste0("Group ", levels(reg_df$group))) +
    scale_fill_manual(values = pal, labels = paste0("Group ", levels(reg_df$group))) +
    labs(title = paste0(tag_short, "CIF"),
         subtitle = subtxt,
         x = "Hours", y = "CIF (95% CI)") +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4), limits = c(0, 1)) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = alpha("white", 0.75), color = "grey80")
    )

  list(ci_obj = ci_obj, p = p)
}

make_plot_event_comp2 <- function(reg_df, tag_short) {
  df <- reg_df %>%
    mutate(event_lab = factor(as.character(fstatus),
                              levels = c("1","2","0"),
                              labels = c("Death","Discharge alive","Censored"))) %>%
    count(group, event_lab) %>%
    group_by(group) %>%
    mutate(p = n / sum(n)) %>%
    ungroup()

  ggplot(df, aes(x = group, y = p, fill = event_lab)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100*x),
                       breaks = c(0, 0.5, 1)) +
    labs(title = paste0("Event types"),
         x = NULL, y = "Percent", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "top",
      legend.text = element_text(size = 8)
    )
}

make_plot_death_time_ecdf <- function(reg_df, tag_short, cfg = cfg) {
  d <- reg_df %>% filter(fstatus == 1)
  if (nrow(d) == 0) return(NULL)

  k <- nlevels(reg_df$group)
  pal <- cfg$resp_cols[seq_len(min(k, length(cfg$resp_cols)))]

  ggplot(d, aes(ftime, color = group)) +
    stat_ecdf(linewidth = 1.1) +
    scale_color_manual(values = pal, labels = paste0("Group ", levels(reg_df$group))) +
    labs(title = paste0(tag_short, "(ECDF)"),
         x = "Hours", y = "ECDF") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = alpha("white", 0.75), color = "grey80")
    )
}

make_table3_fg <- function(fit_pack, ref_group, digits = 3) {
  lv <- levels(fit_pack$data_m1$group)
  if (!(ref_group %in% lv)) ref_group <- lv[1]

  r1 <- extract_group_from_crr(fit_pack$Model1, ref_group, lv)
  r2 <- extract_group_from_crr(fit_pack$Model2, ref_group, lv)
  r3 <- extract_group_from_crr(fit_pack$Model3, ref_group, lv)

  tab_np <- fit_pack$data_m1 %>%
    group_by(group) %>%
    summarise(Participants = n(),
              Deaths = sum(fstatus == 1, na.rm = TRUE),
              .groups = "drop") %>%
    right_join(tibble(group = factor(lv, levels = lv)), by = "group") %>%
    mutate(Participants = ifelse(is.na(Participants), 0L, Participants),
           Deaths = ifelse(is.na(Deaths), 0L, Deaths)) %>%
    arrange(group)

  tibble(group = factor(lv, levels = lv)) %>%
    left_join(tab_np, by = "group") %>%
    mutate(
      M1_txt = mapply(fmt_hrci, r1$HR, r1$L, r1$U, group == ref_group, MoreArgs = list(digits = digits)),
      M1_p   = mapply(fmt_p,   r1$P, group == ref_group),
      M2_txt = mapply(fmt_hrci, r2$HR, r2$L, r2$U, group == ref_group, MoreArgs = list(digits = digits)),
      M2_p   = mapply(fmt_p,   r2$P, group == ref_group),
      M3_txt = mapply(fmt_hrci, r3$HR, r3$L, r3$U, group == ref_group, MoreArgs = list(digits = digits)),
      M3_p   = mapply(fmt_p,   r3$P, group == ref_group)
    ) %>%
    transmute(
      Group = paste0(as.character(group)),
      Deaths, Participants,
      `Model 1 SHR (95% CI)` = M1_txt, `Model 1 P` = M1_p,
      `Model 2 SHR (95% CI)` = M2_txt, `Model 2 P` = M2_p,
      `Model 3 SHR (95% CI)` = M3_txt, `Model 3 P` = M3_p,
      m1 = r1$HR, m1l = r1$L, m1u = r1$U,
      m2 = r2$HR, m2l = r2$L, m2u = r2$U,
      m3 = r3$HR, m3l = r3$L, m3u = r3$U
    )
}

plot_parallel_forest_fg_grob <- function(table3_df, title_short, xlim = cfg$forest_xlim_default) {

  show_df <- table3_df %>%
    mutate(
      Forest1 = strrep("\u00A0", 22),
      Forest2 = strrep("\u00A0", 22),
      Forest3 = strrep("\u00A0", 22)
    ) %>%
    transmute(
      Group = Group,

      `Model 1` = `Model 1 SHR (95% CI)`,
      `P` = `Model 1 P`,
      Forest1 = Forest1,

      `Model 2` = `Model 2 SHR (95% CI)`,
      `P ` = `Model 2 P`,
      Forest2 = Forest2,

      `Model 3` = `Model 3 SHR (95% CI)`,
      `P  ` = `Model 3 P`,
      Forest3 = Forest3
    )

  bg <- rep(c("#FFFFFF", "#F5F5F5"), length.out = nrow(show_df) + 1)

  forest_grob <- grid::grid.grabExpr({
    gt <- forestploter::forest(
      show_df,
      est   = list(table3_df$m1, table3_df$m2, table3_df$m3),
      lower = list(table3_df$m1l, table3_df$m2l, table3_df$m3l),
      upper = list(table3_df$m1u, table3_df$m2u, table3_df$m3u),

      ci_column = c(4, 7, 10),

      ref_line = 1,
      xlim = xlim,
      ticks_at = pretty(xlim, 5),
      title = "",
      bg_col = bg,
      ci_pch = 15,
      ci_lwd = 1.6,
      ref_line_lty = 2
    )
    grid::grid.newpage()
    grid::grid.draw(gt)

    grid::grid.text(
      label = paste0(title_short, " — Fine-Gray"),
      x = 0.5, y = 0.98,
      just = "center",
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  })

  forest_grob
}

combine_ABCD <- function(pA, pB, pC, forest_grob,
                         heights = c(1.05, 2.55),
                         tag_size = cfg$tag_size) {
  if (is.null(pB)) pB <- ggplot() + theme_void() + ggtitle("ECDF not available")

  d <- patchwork::wrap_elements(full = forest_grob)

  ((pA | pB | pC) / d) +
    plot_layout(heights = heights) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = tag_size))
}

run_one <- function(dat_with_group_long,
                    scenario_name,
                    baseline,
                    imp_back = NULL,
                    cfg = cfg,
                    forest_xlim = NULL) {

  cat("\n=============================\n")
  cat("Running: ", scenario_name, "\n\n", sep = "")

  tag_short <- cfg$scenario_label_fun(scenario_name)

  g_code <- make_group_df_from_long(dat_with_group_long, cfg$id_code, cfg$group_col)
  g_stay <- attach_stay_id(g_code, imp_back, cfg$id_code, cfg$id_stay)
  if (!cfg$id_stay %in% names(g_stay)) stop("无法得到 stay_id：请检查 imp_back 映射。")

  reg_df <- build_fg_data(baseline, g_stay, cfg)

  lv <- levels(reg_df$group)
  ref <- cfg$ref_group
  if (!(ref %in% lv)) ref <- lv[1]
  reg_df <- reg_df %>% mutate(group = relevel(group, ref = ref))

  fit_pack <- fit_fg_models(reg_df, cfg)
  table3 <- make_table3_fg(fit_pack, ref, digits = 3)

  pA <- make_plot_cif_death(reg_df, tag_short, cfg)$p
  pB <- make_plot_death_time_ecdf(reg_df, tag_short, cfg)
  pC <- make_plot_event_comp2(reg_df, tag_short)

  xlim_use <- forest_xlim %||% cfg$forest_xlim_default
  forest_grob <- plot_parallel_forest_fg_grob(table3, tag_short, xlim = xlim_use)

  fig <- combine_ABCD(pA, pB, pC, forest_grob,
                      heights = c(1.05, 2.55),
                      tag_size = cfg$tag_size)

  list(
    scenario = scenario_name,
    scenario_label = tag_short,
    reg_df = reg_df,
    fit = fit_pack,
    table3 = table3,
    fig = fig
  )
}

imp_back <- if (exists("imp_back")) imp_back else NULL
stopifnot(exists("baseline"))

need_objs <- c(
  "dat_with_group_LAC_MIN_SF_PF_MIN",
  "dat_with_group_LAC_MEAN_SF_PF_MIN",
  "dat_with_group_LAC_MAX_SF_PF_MIN",
  "dat_with_group_LAC_MIN_SF_PF_MEAN",
  "dat_with_group_LAC_MEAN_SF_PF_MEAN",
  "dat_with_group_LAC_MAX_SF_PF_MEAN"
)
miss_obj <- need_objs[!vapply(need_objs, exists, logical(1))]
if (length(miss_obj) > 0) stop("缺少这些对象：", paste(miss_obj, collapse = ", "))

results <- list()

results[["LAC_MIN_SF_PF_MIN"]]  <- run_one(dat_with_group_LAC_MIN_SF_PF_MIN,  "LAC_MIN_SF_PF_MIN",  baseline, imp_back, cfg, forest_xlim = c(0,4))
results[["LAC_MEAN_SF_PF_MIN"]] <- run_one(dat_with_group_LAC_MEAN_SF_PF_MIN, "LAC_MEAN_SF_PF_MIN", baseline, imp_back, cfg, forest_xlim = c(0,4))
results[["LAC_MAX_SF_PF_MIN"]]  <- run_one(dat_with_group_LAC_MAX_SF_PF_MIN,  "LAC_MAX_SF_PF_MIN",  baseline, imp_back, cfg, forest_xlim = c(0,4))

results[["LAC_MIN_SF_PF_MEAN"]]  <- run_one(dat_with_group_LAC_MIN_SF_PF_MEAN,  "LAC_MIN_SF_PF_MEAN",  baseline, imp_back, cfg, forest_xlim = c(0,4))
results[["LAC_MEAN_SF_PF_MEAN"]] <- run_one(dat_with_group_LAC_MEAN_SF_PF_MEAN, "LAC_MEAN_SF_PF_MEAN", baseline, imp_back, cfg, forest_xlim = c(0,8))
results[["LAC_MAX_SF_PF_MEAN"]]  <- run_one(dat_with_group_LAC_MAX_SF_PF_MEAN,  "LAC_MAX_SF_PF_MEAN",  baseline, imp_back, cfg, forest_xlim = c(0,4))

out_dir <- "FG_Figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (nm in names(results)) {
  ggsave(
    filename = file.path(out_dir, paste0("Figure_", nm, "_ABCD.png")),
    plot = results[[nm]]$fig,
    width = 30, height = 20, units = "cm",
    dpi = 600, limitsize = FALSE
  )
}

cat("Saved to: ", normalizePath(out_dir), "\n", sep = "")










dat_with_group_LAC_MIN_SF_PF_MEAN



imp_out


stay_code <- imp_back %>% select(stay_id,stay_id_code)
baseline3 <- merge(data_baseline,stay_code , by = "stay_id",all.x = TRUE)
baseline3 <- unique(baseline3)
dat_with_group_LAC_MIN_SF_PF_MEAN2 <- merge(dat_with_group_LAC_MIN_SF_PF_MEAN,baseline3,by = "stay_id_code",all.x = TRUE)




colnames(traj_data_imp)
colnames(data_baseline)



library(dplyr)
library(tibble)
library(stringr)
library(purrr)

vars_cont <- c(
  "age_at_admit",
  "icu_los_hours",
  "hosp_los_days"
)

vars_cat <- c(
  "gender",
  "copd","chf","ckd","liver_disease","diabetes","malignancy",
  "death_28d","death_90d",
  "icu_mortality","hosp_mortality"
)

var_labels <- c(
  age_at_admit   = "Age (years)",
  gender         = "Sex (Male)",
  copd           = "COPD",
  chf            = "CHF",
  ckd            = "CKD",
  liver_disease  = "Liver disease",
  diabetes       = "Diabetes",
  malignancy     = "Malignancy",
  icu_los_hours  = "ICU length of stay (hours)",
  hosp_los_days  = "Hospital length of stay (days)",
  death_28d      = "28-day mortality",
  death_90d      = "90-day mortality",
  icu_mortality  = "ICU mortality",
  hosp_mortality = "Hospital mortality"
)

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

fmt_med_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 2)
  paste0(fmt_num(q[2], digits), " (", fmt_num(q[1], digits), ", ", fmt_num(q[3], digits), ")")
}

fmt_n_pct_yes <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  if (!is.factor(x)) x <- factor(x)
  lev <- levels(x)

  yes_lev <- dplyr::case_when(
    "Yes"  %in% lev ~ "Yes",
    "1"    %in% lev ~ "1",
    "TRUE" %in% lev ~ "TRUE",
    "Male" %in% lev ~ "Male",
    TRUE ~ lev[length(lev)]
  )

  n_yes <- sum(x == yes_lev)
  n_all <- length(x)
  pct <- if (n_all == 0) NA_real_ else 100 * n_yes / n_all
  paste0(n_yes, " (", fmt_num(pct, 1), "\\%)")
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  out <- sprintf("%.3f", p)
  sub("^0", "", out)
}

p_kw <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  if (length(x) < 2 || nlevels(g) < 2) return(NA_real_)
  if (length(unique(x)) < 2) return(NA_real_)
  suppressWarnings(tryCatch(kruskal.test(x ~ g)$p.value, error = function(e) NA_real_))
}

p_cat_auto <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- droplevels(as.factor(x[ok]))
  g <- droplevels(as.factor(g[ok]))
  if (nlevels(x) < 2 || nlevels(g) < 2) return(NA_real_)

  tab <- table(x, g)
  chi <- suppressWarnings(tryCatch(chisq.test(tab, correct = FALSE), error = function(e) NULL))
  if (is.null(chi)) return(NA_real_)
  exp <- chi$expected

  suppressWarnings(tryCatch({
    if (any(exp < 5) || any(tab == 0)) {
      sim <- (nrow(tab) > 2 || ncol(tab) > 2)
      fisher.test(tab, simulate.p.value = sim, B = if (sim) 20000 else 0)$p.value
    } else {
      chi$p.value
    }
  }, error = function(e) NA_real_))
}

collapse_to_stay <- function(dat) {
  if (!("stay_id" %in% names(dat))) stop("dat must contain stay_id")
  if (!("group" %in% names(dat))) stop("dat must contain group")

  dat %>%
    arrange(stay_id, day) %>%
    group_by(stay_id) %>%
    summarise(
      group = dplyr::first(group),

      across(all_of(intersect(vars_cont, names(.))), ~ dplyr::first(.x)),

      across(all_of(intersect(vars_cat, names(.))), ~ dplyr::first(.x)),

      .groups = "drop"
    )
}

prep_baseline_data <- function(dat_stay) {
  out <- dat_stay

  glev <- out %>% filter(!is.na(group)) %>% distinct(group) %>% pull(group)
  if (length(glev) == 0) stop("No non-NA group in data.")
  if (all(grepl("^\\d+$", glev))) glev <- as.character(sort(as.integer(glev))) else glev <- sort(glev)
  out <- out %>% mutate(group = factor(as.character(group), levels = glev))

  if ("gender" %in% names(out)) {
    out <- out %>%
      mutate(gender = case_when(
        gender %in% c("M","Male","male",1,"1",TRUE) ~ "Male",
        gender %in% c("F","Female","female",0,"0",FALSE) ~ "Female",
        TRUE ~ as.character(gender)
      )) %>%
      mutate(gender = factor(gender, levels = c("Female","Male")))
  }

  cat_vars_01 <- setdiff(intersect(vars_cat, names(out)), "gender")
  out <- out %>%
    mutate(across(all_of(cat_vars_01), ~ {
      if (is.numeric(.x) || is.integer(.x) || is.logical(.x) || is.character(.x)) {
        x0 <- suppressWarnings(as.numeric(as.character(.x)))
        if (all(na.omit(x0) %in% c(0,1))) {
          factor(x0, levels = c(0,1), labels = c("No","Yes"))
        } else {
          factor(.x)
        }
      } else factor(.x)
    }))

  out
}

build_baseline_obj <- function(dat_stay_prepped) {
  dat <- dat_stay_prepped
  groups <- levels(dat$group)

  n_by_group <- setNames(sapply(groups, function(gg) sum(dat$group == gg, na.rm = TRUE)), groups)
  n_overall  <- sum(!is.na(dat$group))

  make_row <- function(var_name, overall, by_group_named, pval) {
    row_list <- c(
      list(Variable = as.character(var_name),
           Overall  = as.character(overall)),
      as.list(by_group_named),
      list(`p.value` = as.character(pval))
    )
    tibble::as_tibble_row(row_list, .name_repair = "minimal")
  }

  rows <- list()

  for (v in intersect(vars_cont, names(dat))) {
    disp <- if (!is.na(var_labels[v])) var_labels[v] else v
    by_grp <- setNames(sapply(groups, function(gg) fmt_med_iqr(dat[[v]][dat$group == gg])), groups)

    rows[[length(rows) + 1]] <- make_row(
      var_name = disp,
      overall  = fmt_med_iqr(dat[[v]]),
      by_group_named = by_grp,
      pval = fmt_p(p_kw(dat[[v]], dat$group))
    )
  }

  for (v in intersect(vars_cat, names(dat))) {
    disp <- if (!is.na(var_labels[v])) var_labels[v] else v
    by_grp <- setNames(sapply(groups, function(gg) fmt_n_pct_yes(dat[[v]][dat$group == gg])), groups)

    rows[[length(rows) + 1]] <- make_row(
      var_name = disp,
      overall  = fmt_n_pct_yes(dat[[v]]),
      by_group_named = by_grp,
      pval = fmt_p(p_cat_auto(dat[[v]], dat$group))
    )
  }

  df <- bind_rows(rows) %>%
    select(Variable, Overall, all_of(groups), `p.value`)

  list(df = df, groups = groups, n_by_group = n_by_group, n_overall = n_overall)
}

baseline_to_latex <- function(obj, caption, label) {
  df <- obj$df
  groups <- obj$groups

  align <- paste0("l", paste(rep("r", ncol(df) - 1), collapse = ""))

  esc_var <- function(x) {
    x %>%
      str_replace_all("&", "\\\\&") %>%
      str_replace_all("%", "\\\\%") %>%
      str_replace_all("_", "\\\\_")
  }

  header1 <- c("Variable", "Overall", groups, "p-value")
  header2 <- c("", paste0("N = ", obj$n_overall),
               paste0("N = ", unname(obj$n_by_group[groups])), "")

  lines <- c(
    "\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{\\textbf{", caption, "}}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{\\linewidth}{!}{%",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0("\\textbf{", header1, "}", collapse = " & "), " \\\\",
    paste0(header2, collapse = " & "), " \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(df))) {
    r <- df[i, , drop = FALSE]
    r$Variable <- paste0("\\textbf{", esc_var(r$Variable), "}")
    lines <- c(lines, paste(as.character(r), collapse = " & "), " \\\\")
  }

  lines <- c(
    lines,
    "\\bottomrule",
    paste0("\\multicolumn{", ncol(df),
           "}{l}{\\footnotesize Data are presented as median (Q1, Q3) or n (\\%). Continuous variables: Kruskal--Wallis; categorical variables: chi-square or Fisher's exact test as appropriate.}\\\\"),
    "\\end{tabular}}%",
    "\\end{table}"
  )

  paste(lines, collapse = "\n")
}

dat_stay <- collapse_to_stay(dat_with_group_LAC_MIN_SF_PF_MEAN2)
dat_stay <- prep_baseline_data(dat_stay)

obj <- build_baseline_obj(dat_stay)

latex_code <- baseline_to_latex(
  obj,
  caption = "Baseline characteristics by trajectory group",
  label   = "tab:baseline_traj"
)

cat(latex_code)



data_baseline$gender <- ifelse(data_baseline$gender == "M", 1, 0)
