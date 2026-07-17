nwicu_data <-  read.table("data/mimic3/mimic3_sepsis_hypoxemia_ml_features.txt",sep="\t",header=TRUE)

nwicu_traj_data <-  read.table("data/mimic3/mimic3_sepsis_hypoxemia_traj_day_final.txt",sep="\t",header=TRUE)




colnames(nwicu_traj_data)
colnames(nwicu_data)















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

if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")
if (!requireNamespace("mice", quietly = TRUE)) install.packages("mice")
if (!requireNamespace("gtsummary", quietly = TRUE)) install.packages("gtsummary")
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")

suppressPackageStartupMessages({
  library(ggridges)
  library(mice)
  library(gtsummary)
  library(gt)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

stopifnot(exists("nwicu_traj_data"))
stopifnot(exists("nwicu_data"))
stopifnot(exists("all_results"))
stopifnot("LAC_MIN_SF_PF_MEAN" %in% names(all_results))

best_model <- all_results$LAC_MIN_SF_PF_MEAN$best_model
stopifnot(inherits(best_model, "gbmt"))

pred_train <- predict(best_model)
stopifnot(length(pred_train) >= 1)
train_vars <- names(pred_train[[1]])
cat("[Train vars from predict]: ", paste(train_vars, collapse = ", "), "\n")

lac_name <- intersect(c("lactate_min","lactate_mean","lactate_max"), train_vars)
sf_name  <- intersect(c("sf_mean","sf_min"), train_vars)
pf_name  <- intersect(c("pf_mean","pf_min"), train_vars)
if (length(lac_name) == 0) stop("Training predict缺少 lactate_*，实际：", paste(train_vars, collapse=", "))
if (length(sf_name)  == 0) stop("Training predict缺少 sf_*，实际：", paste(train_vars, collapse=", "))
if (length(pf_name)  == 0) stop("Training predict缺少 pf_*，实际：", paste(train_vars, collapse=", "))

lac_name <- lac_name[1]
sf_name  <- sf_name[1]
pf_name  <- pf_name[1]

core3 <- c(lac_name, sf_name, pf_name)
v4_candidates <- setdiff(train_vars, core3)
train_v4 <- if (length(v4_candidates) >= 1) v4_candidates[1] else NULL
need_log_v4 <- if (!is.null(train_v4)) grepl("log", train_v4, ignore.case = TRUE) else FALSE

cat("[Train mapping] lactate=", lac_name,
    " sf=", sf_name, " pf=", pf_name,
    " v4=", train_v4 %||% "(none)",
    " log1p?=", need_log_v4, "\n", sep="")

need_traj_cols <- c("icustay_id","day","lactate_min","sf_mean","pf_mean","ne_total_mcg")
miss_traj <- setdiff(need_traj_cols, names(nwicu_traj_data))
if (length(miss_traj) > 0) stop("nwicu_traj_data 缺少列：", paste(miss_traj, collapse=", "))

nw_long <- nwicu_traj_data %>%
  mutate(
    icustay_id = as.integer(icustay_id),
    day = as.integer(day)
  ) %>%
  filter(day %in% 1:3) %>%
  mutate(
    ne_total_mcg = ifelse(is.na(ne_total_mcg), 0, as.numeric(ne_total_mcg)),
    ne_total_mcg_log = log1p(pmax(ne_total_mcg, 0))
  )

key_vars_for_filter <- c("lactate_min","sf_mean","pf_mean")
keep_ids <- nw_long %>%
  group_by(icustay_id) %>%
  summarise(
    n_day_lac = sum(!is.na(lactate_min)),
    n_day_sf  = sum(!is.na(sf_mean)),
    n_day_pf  = sum(!is.na(pf_mean)),
    .groups = "drop"
  ) %>%
  filter(n_day_lac >= 2, n_day_sf >= 2, n_day_pf >= 2) %>%
  pull(icustay_id)

nw_long <- nw_long %>% filter(icustay_id %in% keep_ids)
cat("[filter] kept icustay_id:", length(unique(nw_long$icustay_id)), "\n")

nw_aligned <- nw_long %>%
  transmute(
    icustay_id,
    day,
    lactate_min = as.numeric(lactate_min),
    sf_mean     = as.numeric(sf_mean),
    pf_mean     = as.numeric(pf_mean)
  )

if (lac_name != "lactate_min") {
  nw_aligned[[lac_name]] <- nw_aligned$lactate_min
  nw_aligned$lactate_min <- NULL
}

if (sf_name != "sf_mean") {
  nw_aligned[[sf_name]] <- nw_aligned$sf_mean
  nw_aligned$sf_mean <- NULL
}
if (pf_name != "pf_mean") {
  nw_aligned[[pf_name]] <- nw_aligned$pf_mean
  nw_aligned$pf_mean <- NULL
}

if (!is.null(train_v4)) {
  nw_aligned[[train_v4]] <- if (need_log_v4) nw_long$ne_total_mcg_log else nw_long$ne_total_mcg
}

need_cols <- c("icustay_id","day", train_vars)
miss_need <- setdiff(need_cols, names(nw_aligned))
if (length(miss_need) > 0) {
  stop("对齐后 nwicu_aligned 缺少训练模型列：", paste(miss_need, collapse=", "),
       "\n当前列：", paste(names(nw_aligned), collapse=", "))
}

nw_aligned <- nw_aligned %>%
  select(all_of(need_cols)) %>%
  arrange(icustay_id, day)

cat("[nw_aligned cols]: ", paste(names(nw_aligned), collapse=", "), "\n")

imp_vars <- train_vars
imp_df <- nw_aligned %>%
  select(icustay_id, day, all_of(imp_vars)) %>%
  mutate(across(where(is.character), as.factor))

set.seed(20260207)
mice_input <- imp_df %>% select(all_of(imp_vars))
mice_methods <- rep("pmm", ncol(mice_input))
names(mice_methods) <- names(mice_input)
mice_fit <- mice::mice(
  data = mice_input,
  m = 5,
  maxit = 10,
  method = mice_methods,
  printFlag = TRUE,
  seed = 20260207
)

nw_aligned_imp <- imp_df %>%
  select(icustay_id, day) %>%
  bind_cols(as_tibble(mice::complete(mice_fit, action = 1))) %>%
  arrange(icustay_id, day)

cat("[mice] done. Method: predictive mean matching; m=5; maxit=10.\n")

nw_use <- nw_aligned_imp

extract_mu_from_gbmt <- function(model, vars, days = 1:3) {
  pr <- predict(model)
  K <- length(pr)
  if (K < 1) stop("predict(model) empty.")
  vars_pred <- names(pr[[1]])
  vars_use <- intersect(vars, vars_pred)
  if (length(vars_use) == 0) stop("predict输出找不到 vars：", paste(vars_pred, collapse=", "))

  v0 <- pr[[1]][[vars_use[1]]]
  Tlen <- length(v0)
  idx <- match(days, 1:Tlen)
  if (any(is.na(idx))) stop("训练预测长度=", Tlen, " 但你要 days=", paste(days, collapse=", "))

  mu <- array(NA_real_, dim = c(K, length(days), length(vars_use)),
              dimnames = list(class=paste0("G",1:K), day=paste0("D",days), var=vars_use))
  for (g in 1:K) {
    for (j in seq_along(vars_use)) {
      mu[g, , j] <- as.numeric(pr[[g]][[vars_use[j]]][idx])
    }
  }
  list(mu=mu, vars=vars_use, days=days, K=K)
}

mu_pack <- extract_mu_from_gbmt(best_model, vars=train_vars, days=1:3)
mu <- mu_pack$mu
K  <- mu_pack$K
vars_used <- mu_pack$vars
days_used <- mu_pack$days

cat("[SEE] vars_used: ", paste(vars_used, collapse=", "), "\n")
cat("[SEE] K=", K, " days=", paste(days_used, collapse=", "), "\n", sep="")

assign_by_SEE <- function(long_df, id_col, day_col, vars, mu, days = 1:3) {
  stopifnot(all(c(id_col, day_col, vars) %in% names(long_df)))
  ids <- sort(unique(long_df[[id_col]]))
  out_list <- vector("list", length(ids))

  for (ii in seq_along(ids)) {
    pid <- ids[ii]
    dfp <- long_df %>%
      filter(.data[[id_col]] == pid, .data[[day_col]] %in% days) %>%
      arrange(.data[[day_col]])

    mat <- matrix(NA_real_, nrow = length(days), ncol = length(vars),
                  dimnames = list(day = days, var = vars))

    for (d in days) {
      row <- dfp %>% filter(.data[[day_col]] == d) %>% slice(1)
      if (nrow(row) == 1) mat[as.character(d), ] <- as.numeric(row[1, vars])
    }

    mse <- rep(NA_real_, dim(mu)[1])
    ncell <- rep(0L, dim(mu)[1])
    for (g in seq_len(dim(mu)[1])) {
      ss <- 0; nn <- 0
      for (td in seq_along(days)) {
        for (j in seq_along(vars)) {
          y <- mat[td, j]
          if (is.finite(y)) {
            ss <- ss + (y - mu[g, td, j])^2
            nn <- nn + 1
          }
        }
      }
      mse[g] <- if (nn > 0) ss / nn else NA_real_
      ncell[g] <- nn
    }

    gbest <- if (all(is.na(mse))) NA_integer_ else which.min(mse)
    out_list[[ii]] <- tibble(
      !!id_col := pid,
      group = as.integer(gbest),
      mse = if (is.na(gbest)) NA_real_ else mse[gbest],
      n_cells = if (is.na(gbest)) 0L else ncell[gbest]
    )
  }
  bind_rows(out_list)
}

nw_assign <- assign_by_SEE(
  long_df = nw_use,
  id_col = "icustay_id",
  day_col = "day",
  vars = vars_used,
  mu = mu,
  days = days_used
)

write.csv(nw_assign, "nwicu_SEE_assignment.csv", row.names = FALSE)
cat("Saved: nwicu_SEE_assignment.csv\n")

nw_with_group <- nw_use %>%
  left_join(nw_assign %>% select(icustay_id, group), by = "icustay_id") %>%
  mutate(group = factor(group))

traj_long <- nw_with_group %>%
  pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
  mutate(day = as.integer(day)) %>%
  filter(!is.na(group))

sum_df <- traj_long %>%
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
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), alpha = 0.20, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2) +
  facet_wrap(~ Variable, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(sum_df$day))) +
  labs(x = "Day", y = "Observed mean (95% CI)", color = "Group") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ridge_long <- nw_with_group %>%
  pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value), !is.na(group))

pB <- ggplot(ridge_long, aes(x = Value, y = group, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.80, scale = 1.0, color = NA) +
  facet_wrap(~ Variable, nrow = 1, scales = "free_x") +
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

ggsave("FIG_nwICU_AB.png", figAB, width = 11, height = 9, dpi = 300)
ggsave("FIG_nwICU_A_TRAJ_CI.png", pA, width = 11, height = 4.2, dpi = 300)
ggsave("FIG_nwICU_B_RIDGE.png", pB, width = 11, height = 4.2, dpi = 300)
cat("Saved: FIG_nwICU_AB.png, FIG_nwICU_A_TRAJ_CI.png, FIG_nwICU_B_RIDGE.png\n")

cfg <- list(
  id = "icustay_id",
  group = "group",
  time_hours = "icu_los_hours",
  death_flag = "hosp_mortality",
  age = "age_at_admit",
  comorb_cols = c("copd","chf","ckd","liver_disease","diabetes","malignancy"),
  ref_group = "1",
  resp_cols = c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"),
  ribbon_alpha = 0.18,
  forest_xlim_default = c(0, 4),
  tag_size = 12
)

stop_has_cols <- function(df, cols, dfname="data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(dfname, " 缺少列: ", paste(miss, collapse=", "))
}

stop_has_cols(nwicu_data, c(cfg$id, cfg$death_flag, cfg$age, cfg$time_hours), "nwicu_data")

reg_df <- nwicu_data %>%
  mutate(!!cfg$id := as.integer(.data[[cfg$id]])) %>%
  inner_join(nw_assign %>% mutate(!!cfg$id := as.integer(.data[[cfg$id]])), by = cfg$id) %>%
  mutate(
    group = factor(as.character(.data[[cfg$group]])),
    death = as.integer(.data[[cfg$death_flag]]),
    ftime = as.numeric(.data[[cfg$time_hours]]),
    fstatus = case_when(
      death == 1 ~ 1L,
      death == 0 ~ 2L,
      TRUE ~ 0L
    ),
    age = as.numeric(.data[[cfg$age]])
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

comorb_ok <- intersect(cfg$comorb_cols, names(reg_df))
for (v in comorb_ok) reg_df[[v]] <- as.numeric(reg_df[[v]])

vars_m1 <- c("ftime","fstatus","group")
vars_m2 <- unique(c(vars_m1, "age"))
vars_m3 <- unique(c(vars_m2, comorb_ok))

d1 <- reg_df %>% filter(complete.cases(across(all_of(vars_m1))))
d2 <- reg_df %>% filter(complete.cases(across(all_of(vars_m2))))
d3 <- reg_df %>% filter(complete.cases(across(all_of(vars_m3))))

d1$group <- relevel(d1$group, ref = ref)
d2$group <- relevel(d2$group, ref = ref)
d3$group <- relevel(d3$group, ref = ref)

X1 <- mkX(d1, c("group"))
X2 <- mkX(d2, c("group","age"))
X3 <- mkX(d3, c("group","age", comorb_ok))

f1 <- safe_crr(d1$ftime, d1$fstatus, X1)
f2 <- safe_crr(d2$ftime, d2$fstatus, X2)
f3 <- safe_crr(d3$ftime, d3$fstatus, X3)

fit_pack <- list(Model1=f1$fit, data_m1=d1, Model2=f2$fit, data_m2=d2, Model3=f3$fit, data_m3=d3)

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
make_plot_cif_death <- function(reg_df, cfg) {
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
  list(ci_obj=ci_obj, p=p)
}

make_plot_death_time_ecdf <- function(reg_df, cfg) {
  d <- reg_df %>% filter(fstatus == 1)
  if (nrow(d) == 0) return(ggplot() + theme_void() + ggtitle("ECDF not available"))
  k <- nlevels(reg_df$group)
  pal <- cfg$resp_cols[seq_len(min(k, length(cfg$resp_cols)))]
  ggplot(d, aes(ftime, color=group)) +
    stat_ecdf(linewidth=1.1) +
    scale_color_manual(values = pal, labels = paste0("Group ", levels(reg_df$group))) +
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
      label = paste0("MIMIC subset — Fine-Gray"),
      x = 0.5, y = 0.98,
      just = "center",
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  })
}

combine_ABCD <- function(pA, pB, pC, forest_grob, tag_size = cfg$tag_size) {
  d <- patchwork::wrap_elements(full = forest_grob)
  ((pA | pB | pC) / d) +
    plot_layout(heights = c(1.05, 2.55)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face="bold", size=tag_size))
}

cif_pack <- make_plot_cif_death(reg_df, cfg)
pA_cif   <- cif_pack$p
pB_ecdf  <- make_plot_death_time_ecdf(reg_df, cfg)
pC_event <- make_plot_event_comp2(reg_df)

table3_fg <- make_table3_fg(fit_pack, ref_group = ref, digits = 3)
forest_grob <- plot_parallel_forest_fg_grob(table3_fg, title_short = "nwICU", xlim = cfg$forest_xlim_default)

figABCD <- combine_ABCD(pA_cif, pB_ecdf, pC_event, forest_grob, tag_size = cfg$tag_size)

out_dir <- "FG_Figures_nwICU"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Figure_nwICU_ABCD.png"), figABCD,
       width = 30, height = 20, units = "cm", dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_nwICU_A_CIF.png"),   pA_cif,   width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure_nwICU_B_ECDF.png"),  pB_ecdf,  width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure_nwICU_C_Event.png"), pC_event, width = 10, height = 7, dpi = 300)

png(file.path(out_dir, "Figure_nwICU_D_Forest.png"), width = 1800, height = 1200, res = 200)
grid::grid.newpage()
grid::grid.draw(forest_grob)
dev.off()

cat("Saved competing risk panels to: ", normalizePath(out_dir), "\n", sep = "")
cat(" - Figure_nwICU_ABCD.png\n")
cat(" - Figure_nwICU_A_CIF.png\n - Figure_nwICU_B_ECDF.png\n - Figure_nwICU_C_Event.png\n - Figure_nwICU_D_Forest.png\n")

need_base_cols <- c("icustay_id","gender","age_at_admit","ethnicity","copd","chf","ckd","liver_disease","diabetes","malignancy",
                    "icu_los_hours","hosp_los_days","icu_mortality","hosp_mortality")
miss_base <- setdiff(need_base_cols, names(nwicu_data))
if (length(miss_base) > 0) {
  cat("WARN: nwicu_data缺少以下Table1列（将自动跳过）： ", paste(miss_base, collapse=", "), "\n", sep="")
}
base_cols_ok <- intersect(need_base_cols, names(nwicu_data))

table1_df <- nwicu_data %>%
  select(all_of(base_cols_ok)) %>%
  mutate(icustay_id = as.integer(icustay_id)) %>%
  left_join(nw_assign %>% select(icustay_id, group), by="icustay_id") %>%
  mutate(
    group = factor(group),
    gender = if ("gender" %in% names(.)) factor(gender) else NULL,
    ethnicity = if ("ethnicity" %in% names(.)) factor(ethnicity) else NULL
  )

drop_low_level_cols <- function(df, by_col) {
  keep <- names(df)
  keep <- setdiff(keep, by_col)
  ok <- sapply(keep, function(v) {
    x <- df[[v]]
    if (all(is.na(x))) return(FALSE)
    if (is.numeric(x)) return(TRUE)
    length(unique(na.omit(x))) >= 2
  })
  keep2 <- c(by_col, keep[ok])
  df[, keep2, drop=FALSE]
}
table1_df2 <- drop_low_level_cols(table1_df, "group")

tbl1 <- table1_df2 %>%
  select(-icustay_id) %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

gt_tbl <- gtsummary::as_gt(tbl1)
latex_code <- gt::as_latex(gt_tbl)

writeLines(latex_code, "Table1_nwICU_byGroup.tex")
cat("Saved: Table1_nwICU_byGroup.tex\n")

assign("nw_use_long", nw_use, envir = .GlobalEnv)
assign("nw_assign", nw_assign, envir = .GlobalEnv)
assign("nw_with_group", nw_with_group, envir = .GlobalEnv)
assign("pA_nw_traj_ci", pA, envir = .GlobalEnv)
assign("pB_nw_ridge", pB, envir = .GlobalEnv)
assign("figAB_nw", figAB, envir = .GlobalEnv)
assign("reg_df_nw", reg_df, envir = .GlobalEnv)
assign("table3_fg_nw", table3_fg, envir = .GlobalEnv)
assign("forest_grob_nw", forest_grob, envir = .GlobalEnv)
assign("figABCD_nw", figABCD, envir = .GlobalEnv)
assign("tbl1_nw", tbl1, envir = .GlobalEnv)

cat("\n================= DONE (nwICU) =================\n")
cat("Model: all_results$LAC_MIN_SF_PF_MEAN$best_model (gbmt)  [NO refit]\n")
cat("Matching: SEE (mean squared error to class mean mu)\n")
cat("Filter: day 1:3 + >=2days non-missing (lactate_min/sf_mean/pf_mean)\n")
cat("Imputation: mice predictive mean matching on model vars (train_vars); m=5\n")
cat("Outputs:\n")
cat(" - nwicu_SEE_assignment.csv\n")
cat(" - FIG_nwICU_AB.png (A: mean+95%CI trajectories; B: ridges)\n")
cat(" - FG_Figures_nwICU/Figure_nwICU_ABCD.png (+ A/B/C/D separate)\n")
cat(" - Table1_nwICU_byGroup.tex\n")
cat("===============================================\n")















traj_long_noCI <- nw_with_group %>%
  pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
  mutate(day = as.integer(day)) %>%
  filter(!is.na(group))

sum_df_noCI <- traj_long_noCI %>%
  group_by(group, day, Variable) %>%
  summarise(
    N = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

pA_noCI <- ggplot(sum_df_noCI, aes(x = day, y = Mean, color = group)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2) +
  facet_wrap(~ Variable, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(sum_df_noCI$day))) +
  labs(x = "Day", y = "Observed mean", color = "Group") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

if (!exists("pB")) {
  ridge_long2 <- nw_with_group %>%
    pivot_longer(cols = all_of(vars_used), names_to = "Variable", values_to = "Value") %>%
    filter(is.finite(Value), !is.na(group))

  pB <- ggplot(ridge_long2, aes(x = Value, y = group, fill = group)) +
    ggridges::geom_density_ridges(alpha = 0.80, scale = 1.0, color = NA) +
    facet_wrap(~ Variable, nrow = 1, scales = "free_x") +
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

figAB_noCI <- (pA_noCI + labs(tag = "A") +
                 theme(plot.tag = element_text(size = 12, face = "bold"),
                       plot.tag.position = c(0, 1))) /
  (pB + labs(tag = "B") +
     theme(plot.tag = element_text(size = 12, face = "bold"),
           plot.tag.position = c(0, 1))) +
  plot_layout(heights = c(1.1, 1))

ggsave("FIG_nwICU_AB_noCI.png", figAB_noCI, width = 11, height = 9, dpi = 300)
ggsave("FIG_nwICU_A_TRAJ_noCI.png", pA_noCI, width = 11, height = 4.2, dpi = 300)

assign("pA_nw_traj_noCI", pA_noCI, envir = .GlobalEnv)
assign("figAB_nw_noCI", figAB_noCI, envir = .GlobalEnv)

cat("Saved: FIG_nwICU_AB_noCI.png, FIG_nwICU_A_TRAJ_noCI.png\n")






















dat_with_group_LAC_MIN_SF_PF_MEAN

