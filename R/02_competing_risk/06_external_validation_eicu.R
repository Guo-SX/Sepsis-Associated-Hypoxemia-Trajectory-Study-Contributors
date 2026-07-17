eicu_data <-  read.table("data/eicu/eicu_sepsis_hypoxemia_ml_features.txt",sep="\t",header=TRUE)

eicu_traj_data <-  read.table("data/eicu/eicu_sepsis_hypoxemia_traj_day_final_v2.txt",sep="\t",header=TRUE)


colnames(eicu_traj_data)
colnames(eicu_data)



methods("posterior")
methods("predict")

args(posterior)
args(predict)





res_target  <- all_results[["LAC_MIN_SF_PF_MEAN"]]$res
best_key    <- res_target$best_key
model_train <- res_target$best_model
varNames    <- res_target$varNames





res$ne_total_mcg_log







getwd()






suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(cmprsk)
  library(forestploter)
  library(grid)
  library(scales)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")
suppressPackageStartupMessages(library(ggridges))

if (!requireNamespace("mice", quietly = TRUE)) install.packages("mice")
suppressPackageStartupMessages(library(mice))

if (!requireNamespace("gtsummary", quietly = TRUE)) install.packages("gtsummary")
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
suppressPackageStartupMessages({
  library(gtsummary)
  library(gt)
})

stopifnot(exists("eicu_traj_data"))
stopifnot(exists("eicu_data"))
stopifnot(exists("all_results"))
stopifnot("LAC_MIN_SF_PF_MEAN" %in% names(all_results))

best_model <- all_results$LAC_MIN_SF_PF_MEAN$best_model
stopifnot(inherits(best_model, "gbmt"))

pred_train <- predict(best_model)
G <- length(pred_train)
if (G < 1) stop("predict(best_model) is empty.")

train_vars_pred <- names(pred_train[[1]])
cat("[Train vars from predict]: ", paste(train_vars_pred, collapse=", "), "\n")

need_3 <- c("lactate_min", "sf_mean", "pf_mean")
if (!all(need_3 %in% train_vars_pred)) {
  stop("这个 best_model 的 predict 输出不包含 lactate_min/sf_mean/pf_mean。\n",
       "实际 predict vars = ", paste(train_vars_pred, collapse=", "), "\n",
       "=> 你现在的 best_model 不是你要的那个（乳酸min+sf/pf均值）模型。")
}

v4_candidates <- setdiff(train_vars_pred, need_3)
train_v4 <- NULL
if (length(v4_candidates) >= 1) {
  dose_like <- v4_candidates[grepl("mcg|dose|vaso|norepi|epi|ne_", v4_candidates, ignore.case = TRUE)]
  train_v4 <- (dose_like %||% v4_candidates)[1]
}
if (!is.null(train_v4)) cat("[Train 4th var detected]: ", train_v4, "\n") else cat("[Train 4th var]: none\n")

need_log_v4 <- if (!is.null(train_v4)) grepl("log", train_v4, ignore.case = TRUE) else FALSE

need_cols_traj <- c("patientunitstayid","day","lactate_min","sf_mean","pf_mean","epi_total_mcg")
miss_traj <- setdiff(need_cols_traj, names(eicu_traj_data))
if (length(miss_traj) > 0) stop("eicu_traj_data 缺少列：", paste(miss_traj, collapse=", "))

eicu_long <- eicu_traj_data %>%
  transmute(
    patientunitstayid = as.integer(patientunitstayid),
    day = as.integer(day),
    lactate_min = as.numeric(lactate_min),
    sf_mean = as.numeric(sf_mean),
    pf_mean = as.numeric(pf_mean),
    epi_total_mcg = as.numeric(epi_total_mcg)
  ) %>%
  filter(day %in% 1:3) %>%
  mutate(
    epi_total_mcg = ifelse(is.na(epi_total_mcg), 0, epi_total_mcg)
  ) %>%
  arrange(patientunitstayid, day)

keep_ids <- eicu_long %>%
  group_by(patientunitstayid) %>%
  summarise(
    n_day_lac = sum(!is.na(lactate_min)),
    n_day_sf  = sum(!is.na(sf_mean)),
    n_day_pf  = sum(!is.na(pf_mean)),
    .groups = "drop"
  ) %>%
  filter(n_day_lac >= 2, n_day_sf >= 2, n_day_pf >= 2) %>%
  pull(patientunitstayid)

eicu_long <- eicu_long %>% filter(patientunitstayid %in% keep_ids)
cat("[filter] kept patients:", dplyr::n_distinct(eicu_long$patientunitstayid), "\n")

eicu_aligned <- eicu_long %>%
  transmute(
    patientunitstayid,
    day,
    lactate_min,
    sf_mean,
    pf_mean
  )

if (!is.null(train_v4)) {
  v4_src <- eicu_long$epi_total_mcg
  eicu_aligned[[train_v4]] <- if (need_log_v4) log1p(pmax(v4_src, 0)) else v4_src
}

vars_used <- c(need_3, train_v4) %>% unique() %>% .[!is.na(.)]

cat("[eICU aligned vars]: ", paste(vars_used, collapse=", "), "\n")

set.seed(20260207)

imp_vars <- intersect(vars_used, names(eicu_aligned))
if (length(imp_vars) == 0) stop("没有可插补的变量列。")

imp_input <- eicu_aligned %>%
  select(all_of(imp_vars)) %>%
  mutate(across(where(is.character), as.factor))

mice_methods <- rep("pmm", ncol(imp_input))
names(mice_methods) <- names(imp_input)
mice_fit <- mice::mice(
  data = imp_input,
  m = 5,
  maxit = 10,
  method = mice_methods,
  printFlag = TRUE,
  seed = 20260207
)

eicu_aligned_imp <- eicu_aligned %>%
  select(patientunitstayid, day) %>%
  bind_cols(as_tibble(mice::complete(mice_fit, action = 1)))

eicu_aligned <- eicu_aligned_imp

cat("[mice] done. Method: predictive mean matching; m=5; maxit=10.\n")

days_use <- 1:3

vec_to_df <- function(vec, group_id, var_out) {
  day <- suppressWarnings(as.integer(names(vec)))
  if (all(is.na(day))) day <- seq_along(vec)
  tibble(group = group_id, day = day, !!var_out := as.numeric(vec))
}

train_mu <- bind_rows(lapply(seq_len(G), function(k) {
  pk <- pred_train[[k]]
  missk <- setdiff(vars_used, names(pk))
  if (length(missk) > 0) {
    stop("pred_train[[", k, "]] 缺少变量：", paste(missk, collapse=", "),
         "\n实际有：", paste(names(pk), collapse=", "))
  }

  dfs <- list(
    vec_to_df(pk[["lactate_min"]], k, "lactate_min"),
    vec_to_df(pk[["sf_mean"]],     k, "sf_mean"),
    vec_to_df(pk[["pf_mean"]],     k, "pf_mean")
  )
  if (!is.null(train_v4)) dfs <- c(dfs, list(vec_to_df(pk[[train_v4]], k, train_v4)))

  out <- Reduce(function(a,b) inner_join(a,b, by=c("group","day")), dfs) %>%
    filter(day %in% days_use) %>%
    arrange(group, day)

  out
}))

cat("[train_mu head]\n")
print(head(train_mu))

assign_by_SEE <- function(long_df, train_mu, id_col="patientunitstayid", day_col="day", vars, days=1:3) {
  stopifnot(all(c(id_col, day_col, vars) %in% names(long_df)))
  stopifnot(all(c("group","day", vars) %in% names(train_mu)))

  ids <- sort(unique(long_df[[id_col]]))
  groups <- sort(unique(train_mu$group))

  mu_list <- lapply(groups, function(g) {
    train_mu %>% filter(group == g, day %in% days) %>% arrange(day)
  })
  names(mu_list) <- as.character(groups)

  out <- vector("list", length(ids))

  for (i in seq_along(ids)) {
    pid <- ids[i]
    dfp <- long_df %>% filter(.data[[id_col]] == pid, .data[[day_col]] %in% days) %>%
      arrange(.data[[day_col]])

    mat <- matrix(NA_real_, nrow = length(days), ncol = length(vars),
                  dimnames = list(day = days, var = vars))
    for (d in days) {
      row <- dfp %>% filter(.data[[day_col]] == d) %>% slice(1)
      if (nrow(row) == 1) mat[as.character(d), ] <- as.numeric(row[1, vars])
    }

    mse <- rep(NA_real_, length(groups))
    nn  <- rep(0L, length(groups))

    for (gi in seq_along(groups)) {
      g <- groups[gi]
      mug <- mu_list[[as.character(g)]]
      ss <- 0; ncell <- 0
      for (td in seq_along(days)) {
        for (vj in seq_along(vars)) {
          y <- mat[td, vj]
          if (is.finite(y)) {
            muval <- mug[[vars[vj]]][td]
            ss <- ss + (y - muval)^2
            ncell <- ncell + 1
          }
        }
      }
      mse[gi] <- if (ncell > 0) ss / ncell else NA_real_
      nn[gi]  <- ncell
    }

    gbest_idx <- if (all(is.na(mse))) NA_integer_ else which.min(mse)

    out[[i]] <- tibble(
      patientunitstayid = pid,
      group = if (is.na(gbest_idx)) NA_integer_ else groups[gbest_idx],
      mse = if (is.na(gbest_idx)) NA_real_ else mse[gbest_idx],
      n_cells = if (is.na(gbest_idx)) 0L else nn[gbest_idx]
    )
  }

  bind_rows(out)
}

eicu_assign <- assign_by_SEE(
  long_df = eicu_aligned,
  train_mu = train_mu,
  vars = vars_used,
  days = days_use
)

write.csv(eicu_assign, "eicu_SEE_assignment.csv", row.names = FALSE)
cat("Saved: eicu_SEE_assignment.csv\n")

eicu_with_group <- eicu_aligned %>%
  left_join(eicu_assign %>% select(patientunitstayid, group), by = "patientunitstayid") %>%
  mutate(group = factor(group))

var_labels_plain <- setNames(vars_used, vars_used)
var_labels_plain["lactate_min"] <- "Lactate (min)"
var_labels_plain["sf_mean"] <- "SpO2/FiO2 (mean)"
var_labels_plain["pf_mean"] <- "PaO2/FiO2 (mean)"
if (!is.null(train_v4)) {
  var_labels_plain[train_v4] <- if (need_log_v4) "Dose (log1p mcg)" else "Dose (mcg)"
}

facet_labeller <- function(mapping) {
  function(x) {
    out <- mapping[x]
    out[is.na(out)] <- x[is.na(out)]
    unname(out)
  }
}

sum_df <- eicu_with_group %>%
  pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
  mutate(day = as.integer(day)) %>%
  group_by(group, day, Variable) %>%
  summarise(
    N = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    Sd = sd(Value, na.rm = TRUE),
    Se = Sd / sqrt(pmax(N, 1)),
    Lwr = Mean - 1.96 * Se,
    Upr = Mean + 1.96 * Se,
    .groups = "drop"
  )

pA <- ggplot(sum_df, aes(x = day, y = Mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), alpha = 0.20, linewidth = 0, show.legend = FALSE) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2) +
  facet_wrap(~ Variable, nrow = 1, scales = "free_y",
             labeller = as_labeller(facet_labeller(var_labels_plain))) +
  scale_x_continuous(breaks = sort(unique(sum_df$day))) +
  labs(x = "Day", y = "Mean (95% CI)", color = "Group") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ridge_long <- eicu_with_group %>%
  pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value))

pB <- ggplot(ridge_long, aes(x = Value, y = group, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.80, scale = 1.0, color = NA) +
  facet_wrap(~ Variable, nrow = 1, scales = "free_x",
             labeller = as_labeller(facet_labeller(var_labels_plain))) +
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

figAB <- (pA + labs(tag = "A") +
            theme(plot.tag = element_text(size = 12, face = "bold"),
                  plot.tag.position = c(0, 1))) /
  (pB + labs(tag = "B") +
     theme(plot.tag = element_text(size = 12, face = "bold"),
           plot.tag.position = c(0, 1))) +
  plot_layout(heights = c(1.1, 1))

ggsave("FIG_eICU_AB.png", figAB, width = 11, height = 9, dpi = 300)
ggsave("FIG_eICU_A_TRAJ_CI.png", pA, width = 11, height = 4.2, dpi = 300)
ggsave("FIG_eICU_B_RIDGE.png", pB, width = 11, height = 4.2, dpi = 300)
cat("Saved: FIG_eICU_AB.png, FIG_eICU_A_TRAJ_CI.png, FIG_eICU_B_RIDGE.png\n")

assign("eicu_with_group", eicu_with_group, envir = .GlobalEnv)
assign("pA_traj_ci", pA, envir = .GlobalEnv)
assign("pB_ridge", pB, envir = .GlobalEnv)
assign("figAB_trajCI_ridge", figAB, envir = .GlobalEnv)

need_cols_base <- c("patientunitstayid","gender","ethnicity","age_years",
                    "copd","chf","ckd","liver_disease","diabetes","malignancy",
                    "icu_los_hours","hosp_los_days","icu_mortality","hosp_mortality")
miss_base <- setdiff(need_cols_base, names(eicu_data))
if (length(miss_base) > 0) {
  warning("eicu_data 缺少部分基线列，Table1 将跳过缺失列：", paste(miss_base, collapse=", "))
}

base_cols_ok <- intersect(need_cols_base, names(eicu_data))

table1_df <- eicu_data %>%
  select(all_of(base_cols_ok)) %>%
  mutate(patientunitstayid = as.integer(patientunitstayid)) %>%
  inner_join(eicu_assign %>% mutate(patientunitstayid = as.integer(patientunitstayid)),
             by = "patientunitstayid") %>%
  mutate(group = factor(group))

tbl1 <- table1_df %>%
  select(-patientunitstayid) %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_overall() %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test")) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

tex_code <- gtsummary::as_latex(tbl1)
writeLines(tex_code, "Table1_eICU_by_group.tex")
cat("Saved: Table1_eICU_by_group.tex (Overleaf LaTeX)\n")

cfg <- list(
  id = "patientunitstayid",
  group = "group",
  time_hours = "icu_los_hours",
  death_flag = "hosp_mortality",
  age = "age_years",
  cbc_cols = c(),
  comorb_cols = c("copd","chf","ckd","liver_disease","diabetes","malignancy"),
  ref_group = "1",
  ribbon_alpha = 0.18,
  forest_xlim_default = c(0, 4),
  tag_size = 12
)

stop_has_cols <- function(df, cols, dfname="data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(dfname, " 缺少列: ", paste(miss, collapse=", "))
}

stop_has_cols(eicu_data, c(cfg$id, cfg$death_flag, cfg$age, cfg$time_hours), "eicu_data")

reg_df <- eicu_data %>%
  mutate(!!cfg$id := as.integer(.data[[cfg$id]])) %>%
  inner_join(eicu_assign %>% mutate(!!cfg$id := as.integer(.data[[cfg$id]])),
             by = cfg$id) %>%
  mutate(
    group = factor(as.character(.data[[cfg$group]])),
    death = as.integer(.data[[cfg$death_flag]]),
    ftime = as.numeric(.data[[cfg$time_hours]]),
    fstatus = case_when(
      death == 1 ~ 1L,
      death == 0 ~ 2L,
      TRUE ~ 0L
    )
  ) %>%
  filter(is.finite(ftime), ftime > 0, fstatus %in% c(0,1,2)) %>%
  droplevels()

lv <- levels(reg_df$group)
ref <- cfg$ref_group
if (!(ref %in% lv)) ref <- lv[1]
reg_df <- reg_df %>% mutate(group = relevel(group, ref = ref))

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

safe_crr <- function(ftime, fstatus, X, failcode=1, cencode=0) {
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

cbc_ok    <- intersect(cfg$cbc_cols, names(reg_df))
comorb_ok <- intersect(cfg$comorb_cols, names(reg_df))

reg_df <- reg_df %>% mutate(age = as.numeric(.data[[cfg$age]]))
for (v in cbc_ok)    reg_df[[v]] <- as.numeric(reg_df[[v]])
for (v in comorb_ok) reg_df[[v]] <- as.numeric(reg_df[[v]])

vars_m1 <- c("ftime","fstatus","group")
vars_m2 <- unique(c(vars_m1, "age", cbc_ok))
vars_m3 <- unique(c(vars_m2, comorb_ok))

d1 <- reg_df %>% filter(complete.cases(across(all_of(vars_m1))))
d2 <- reg_df %>% filter(complete.cases(across(all_of(vars_m2))))
d3 <- reg_df %>% filter(complete.cases(across(all_of(vars_m3))))

d1$group <- relevel(d1$group, ref = ref)
d2$group <- relevel(d2$group, ref = ref)
d3$group <- relevel(d3$group, ref = ref)

X1 <- mkX(d1, c("group"))
X2 <- mkX(d2, c("group","age", cbc_ok))
X3 <- mkX(d3, c("group","age", cbc_ok, comorb_ok))

f1 <- safe_crr(d1$ftime, d1$fstatus, X1)
f2 <- safe_crr(d2$ftime, d2$fstatus, X2)
f3 <- safe_crr(d3$ftime, d3$fstatus, X3)

fit_pack <- list(Model1=f1$fit, data_m1=d1, Model2=f2$fit, data_m2=d2, Model3=f3$fit, data_m3=d3)

fmt_p2 <- function(p, is_ref=FALSE) {
  if (is_ref) return("")
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_hrci <- function(hr,l,u,is_ref=FALSE,digits=3){
  if (is_ref) return("Reference")
  if (any(!is.finite(c(hr,l,u))) || any(is.na(c(hr,l,u)))) return("")
  sprintf(paste0("%.",digits,"f (%.",digits,"f–%.",digits,"f)"), hr,l,u)
}

extract_group_from_crr <- function(fit, ref_group, group_levels) {
  sm <- summary(fit)
  coef_tab <- sm$coef
  if (is.null(coef_tab) || nrow(coef_tab) == 0) {
    return(tibble(group=factor(group_levels, levels=group_levels),
                  HR=NA_real_, L=NA_real_, U=NA_real_, P=NA_real_))
  }
  rn <- rownames(coef_tab)
  b  <- as.numeric(coef_tab[,1])
  se <- as.numeric(coef_tab[,3])
  p  <- suppressWarnings(as.numeric(coef_tab[, ncol(coef_tab)]))
  names(b)<-names(se)<-names(p)<-rn

  out <- tibble(group=factor(group_levels, levels=group_levels),
                HR=NA_real_, L=NA_real_, U=NA_real_, P=NA_real_)
  out[out$group==ref_group, c("HR","L","U")] <- 1
  out[out$group==ref_group, "P"] <- NA_real_

  for (g in group_levels) {
    if (g==ref_group) next
    key <- paste0("group", g)
    if (!key %in% rn) next
    bi <- b[key]; sei <- se[key]
    out[out$group==g, "HR"] <- exp(bi)
    out[out$group==g, "L"]  <- exp(bi - 1.96*sei)
    out[out$group==g, "U"]  <- exp(bi + 1.96*sei)
    out[out$group==g, "P"]  <- p[key]
  }
  out
}

make_table3_fg <- function(fit_pack, ref_group, digits=3) {
  lv <- levels(fit_pack$data_m1$group)
  if (!(ref_group %in% lv)) ref_group <- lv[1]

  r1 <- extract_group_from_crr(fit_pack$Model1, ref_group, lv)
  r2 <- extract_group_from_crr(fit_pack$Model2, ref_group, lv)
  r3 <- extract_group_from_crr(fit_pack$Model3, ref_group, lv)

  tab_np <- fit_pack$data_m1 %>%
    group_by(group) %>%
    summarise(Participants=n(), Deaths=sum(fstatus==1, na.rm=TRUE), .groups="drop") %>%
    right_join(tibble(group=factor(lv, levels=lv)), by="group") %>%
    mutate(Participants=ifelse(is.na(Participants),0L,Participants),
           Deaths=ifelse(is.na(Deaths),0L,Deaths)) %>%
    arrange(group)

  tibble(group=factor(lv, levels=lv)) %>%
    left_join(tab_np, by="group") %>%
    mutate(
      M1_txt = mapply(fmt_hrci, r1$HR, r1$L, r1$U, group==ref_group, MoreArgs=list(digits=digits)),
      M1_p   = mapply(fmt_p2,   r1$P,  group==ref_group),
      M2_txt = mapply(fmt_hrci, r2$HR, r2$L, r2$U, group==ref_group, MoreArgs=list(digits=digits)),
      M2_p   = mapply(fmt_p2,   r2$P,  group==ref_group),
      M3_txt = mapply(fmt_hrci, r3$HR, r3$L, r3$U, group==ref_group, MoreArgs=list(digits=digits)),
      M3_p   = mapply(fmt_p2,   r3$P,  group==ref_group)
    ) %>%
    transmute(
      Group = as.character(group),
      Deaths, Participants,
      `Model 1 SHR (95% CI)` = M1_txt, `Model 1 P` = M1_p,
      `Model 2 SHR (95% CI)` = M2_txt, `Model 2 P` = M2_p,
      `Model 3 SHR (95% CI)` = M3_txt, `Model 3 P` = M3_p,
      m1=r1$HR,m1l=r1$L,m1u=r1$U,
      m2=r2$HR,m2l=r2$L,m2u=r2$U,
      m3=r3$HR,m3l=r3$L,m3u=r3$U
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
      `Model 1` = `Model 1 SHR (95% CI)`, `P` = `Model 1 P`, Forest1 = Forest1,
      `Model 2` = `Model 2 SHR (95% CI)`, `P `= `Model 2 P`, Forest2 = Forest2,
      `Model 3` = `Model 3 SHR (95% CI)`, `P  `= `Model 3 P`, Forest3 = Forest3
    )

  bg <- rep(c("#FFFFFF", "#F5F5F5"), length.out = nrow(show_df) + 1)

  grid::grid.grabExpr({
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
  nms <- nms[keep]
  if (length(nms) == 0) return(tibble(group=factor(), time=numeric(), est=numeric(), lwr=numeric(), upr=numeric()))

  out <- lapply(nms, function(nm) {
    x <- ci_obj[[nm]]
    grp <- sub(paste0(" ", event_code, "$"), "", nm)
    est <- as.numeric(x$est)
    tt  <- as.numeric(x$time)
    vv  <- suppressWarnings(as.numeric(x$var))
    se  <- if (!is.null(vv) && length(vv) == length(est)) sqrt(pmax(vv, 0)) else rep(NA_real_, length(est))
    lwr <- pmax(est - 1.96*se, 0)
    upr <- pmin(est + 1.96*se, 1)
    tibble(group=grp, time=tt, est=est, lwr=lwr, upr=upr)
  })
  df <- bind_rows(out)
  df$group <- factor(df$group, levels = unique(df$group))
  df
}

make_plot_cif_death <- function(reg_df) {
  ci_obj <- cmprsk::cuminc(reg_df$ftime, reg_df$fstatus, reg_df$group, cencode = 0)
  cif_df <- tidy_cuminc_ci(ci_obj, event_code = 1)
  gray_p <- get_gray_p(ci_obj)
  subtxt <- if (!is.na(gray_p)) sprintf("Gray p %.3g", gray_p) else NULL

  ggplot(cif_df, aes(time, est, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = cfg$ribbon_alpha, linewidth = 0) +
    geom_step(linewidth = 1.1) +
    labs(title = "CIF", subtitle = subtxt, x = "Hours", y = "CIF (95% CI)") +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4), limits = c(0,1)) +
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
}

make_plot_event_comp2 <- function(reg_df) {
  df <- reg_df %>%
    mutate(event_lab = factor(as.character(fstatus),
                              levels=c("1","2","0"),
                              labels=c("Death","Discharge alive","Censored"))) %>%
    count(group, event_lab) %>%
    group_by(group) %>%
    mutate(p = n / sum(n)) %>%
    ungroup()

  ggplot(df, aes(x=group, y=p, fill=event_lab)) +
    geom_col(width=0.72) +
    coord_flip() +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100*x),
                       breaks = c(0,0.5,1)) +
    labs(title="Event types", x=NULL, y="Percent", fill=NULL) +
    theme_bw(base_size=10) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=10),
      legend.position = "top",
      legend.text = element_text(size=8)
    )
}

make_plot_death_time_ecdf <- function(reg_df) {
  d <- reg_df %>% filter(fstatus == 1)
  if (nrow(d) == 0) return(ggplot() + theme_void() + ggtitle("ECDF not available"))
  ggplot(d, aes(ftime, color=group)) +
    stat_ecdf(linewidth=1.1) +
    labs(title="ECDF (Deaths only)", x="Hours", y="ECDF") +
    theme_bw(base_size=10) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=10),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title = element_blank(),
      legend.text = element_text(size=8),
      legend.background = element_rect(fill = alpha("white", 0.75), color = "grey80")
    )
}

combine_ABCD <- function(pA, pB, pC, forest_grob, tag_size = cfg$tag_size) {
  d <- patchwork::wrap_elements(full = forest_grob)
  ((pA | pB | pC) / d) +
    plot_layout(heights = c(1.05, 2.55)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face="bold", size=tag_size))
}

table3 <- make_table3_fg(fit_pack, ref_group = ref, digits = 3)
forest_grob <- plot_parallel_forest_fg_grob(table3, title_short = "eICU", xlim = cfg$forest_xlim_default)

pCIF  <- make_plot_cif_death(reg_df)
pECDF <- make_plot_death_time_ecdf(reg_df)
pEvent<- make_plot_event_comp2(reg_df)

figABCD <- combine_ABCD(pCIF, pECDF, pEvent, forest_grob, tag_size = cfg$tag_size)

out_dir <- "FG_Figures_eICU"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Figure_eICU_ABCD.png"), figABCD,
       width = 30, height = 20, units = "cm", dpi = 600, limitsize = FALSE)

cat("Saved: ", normalizePath(out_dir), "/Figure_eICU_ABCD.png\n", sep="")

cat("\n================= DONE =================\n")
cat("Model: all_results$LAC_MIN_SF_PF_MEAN$best_model (gbmt)  [NO refit]\n")
cat("Filter: keep >=2 non-missing days for lactate_min/sf_mean/pf_mean (day1:3)\n")
cat("Imputation: mice predictive mean matching on vars_used (exclude id/day); m=5\n")
cat("Matching: SEE (min MSE to class mean trajectories)\n")
cat("vars_used: ", paste(vars_used, collapse=", "), "\n", sep="")
cat("Outputs:\n")
cat(" - eicu_SEE_assignment.csv\n")
cat(" - FIG_eICU_AB.png (A: mean±95%CI trajectories; B: ridges)\n")
cat(" - FIG_eICU_A_TRAJ_CI.png, FIG_eICU_B_RIDGE.png\n")
cat(" - FG_Figures_eICU/Figure_eICU_ABCD.png\n")
cat(" - Table1_eICU_by_group.tex (Overleaf)\n")
cat("Objects in GlobalEnv:\n")
cat(" - eicu_with_group, pA_traj_ci, pB_ridge, figAB_trajCI_ridge\n")
cat("========================================\n")
















base_vars <- c(
  "gender","ethnicity","age_years",
  "copd","chf","ckd","liver_disease","diabetes","malignancy",
  "icu_los_hours","hosp_los_days",
  "icu_mortality","hosp_mortality"
)
base_vars <- intersect(base_vars, names(eicu_data))

table1_df <- eicu_data %>%
  mutate(patientunitstayid = as.integer(patientunitstayid)) %>%
  inner_join(eicu_assign %>% mutate(patientunitstayid = as.integer(patientunitstayid)),
             by = "patientunitstayid") %>%
  transmute(
    group = factor(group),
    across(all_of(base_vars))
  )

is_single_level <- function(x) {
  x2 <- x[!is.na(x)]
  if (length(x2) == 0) return(TRUE)
  length(unique(x2)) < 2
}

drop_vars <- names(table1_df)[sapply(table1_df, is_single_level)]
drop_vars <- setdiff(drop_vars, "group")
if (length(drop_vars) > 0) {
  message("[Table1] Drop single-level vars: ", paste(drop_vars, collapse=", "))
  table1_df <- table1_df %>% select(-all_of(drop_vars))
}

tbl1 <- table1_df %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    )
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

gt_tbl <- tbl1 %>% as_gt()

gt::gtsave(gt_tbl, "Table1_eICU_by_group.tex")
cat("Saved: Table1_eICU_by_group.tex\n")

assign("tbl1", tbl1, envir = .GlobalEnv)
assign("gt_tbl1", gt_tbl, envir = .GlobalEnv)







