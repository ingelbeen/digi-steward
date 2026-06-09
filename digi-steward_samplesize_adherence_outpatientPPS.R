#######################################################################
# digi-steward                                                        #
# stepped wedge power estimation for primary outcome: proportion with #
# correct antibiotic intake among patients w/ appropriately dispensed #
# antibiotic treatment                                                #
#######################################################################

# define parameters
T_ <- 12 # number of time periods (months); crossover steps = T_ - 1, each batch crosses once.
p0       <- 0.50    # baseline correct intake
p1       <- 0.65    # intervention correct intake (target)
alpha    <- 0.05    # two-sided
icc      <- 0.05    # intraclass correlation coefficient -> based on clustering of AMU in Ingelbeen, Valia et al Lancet Infect Dis 2026
power_target <- 0.80
prop_antibiotic   <- 0.50   # proportion of surveyed patients receiving antibiotic
prop_appropriate  <- 0.40   # proportion of those with appropriate dispensing
eligibility_rate  <- prop_antibiotic * prop_appropriate   # = 0.20
elig_rate <- 0.5 * 0.4

# two survey scenarios - 20 surveys per provider per month (interviewer spends one day at each provider, each month) or 40 per provider (two survey days per month)
surveys_per_month <- c(20, 40)
m_values          <- surveys_per_month * eligibility_rate   # eligible per provider per month
# => m = 4 (for 20 surveys) and m = 8 (for 40 surveys)

sw_power <- function(K, T_, m, p0, p1, icc, alpha = 0.05) {

  delta    <- p1 - p0
  p_bar    <- (p0 + p1) / 2
  sigma2   <- p_bar * (1 - p_bar)          # total variance (average p)
  sigma_w2 <- sigma2 * (1 - icc)            # within-cluster variance
  sigma_b2 <- sigma2 * icc                  # between-cluster variance

  # Cluster-period mean variance / covariance
  var_cp   <- sigma_b2 + sigma_w2 / m       # Var(Y-bar_kt)
  cov_cp   <- sigma_b2                      # Cov(Y-bar_kt, Y-bar_kt'), t != t'
# treatment assignment matrix (K x T_): cluster k crosses over after step ceil(k / clusters_per_step).
  steps              <- T_ - 1
  clusters_per_step  <- K / steps

  X_assign <- matrix(0, nrow = K, ncol = T_)
  for (k in seq_len(K)) {
    step_k <- ceiling(k / clusters_per_step)   # 1-indexed crossover step
    # Cluster k is in intervention from period (step_k + 1) onwards (1-indexed)
    # i.e., column indices (step_k + 1):T_ in 1-indexed = step_k:T_ 0-indexed
    if (step_k < T_) {
      X_assign[k, (step_k + 1):T_] <- 1
    } }
  # cluster covariance matrix (T_ x T_) 
  Sigma_k              <- matrix(cov_cp, nrow = T_, ncol = T_)
  diag(Sigma_k)        <- var_cp
  Sigma_k_inv          <- solve(Sigma_k)
  # Design matrix for cluster k: T_ rows, (T_ + 1) columns
  # Columns 1:T_ = time period indicators (fixed effects); column T_+1 = treatment
  # Sum information matrices over all K clusters (clusters are independent)
  n_params <- T_ + 1L
  Info     <- matrix(0, nrow = n_params, ncol = n_params)
  for (k in seq_len(K)) {
    D_k           <- matrix(0, nrow = T_, ncol = n_params)
    diag(D_k)     <- 1                        # time fixed effects
    D_k[, n_params] <- X_assign[k, ]         # treatment column
    Info          <- Info + t(D_k) %*% Sigma_k_inv %*% D_k
  }
  # Guard: Info is singular when K < T_ (fewer clusters than time periods means
  # not all time fixed effects are estimable). Return NA gracefully.
  Info_inv <- tryCatch(
    solve(Info),
    error = function(e) NULL)
  if (is.null(Info_inv)) {
    return(list(power = NA_real_, se_delta = NA_real_,
                var_delta = NA_real_, design_matrix = X_assign))
  }
  var_delta <- Info_inv[n_params, n_params]
  se_delta  <- sqrt(var_delta)
  z_alpha   <- qnorm(1 - alpha / 2)
  power     <- pnorm(abs(delta) / se_delta - z_alpha)
  list(
    power          = power,
    se_delta       = se_delta,
    var_delta      = var_delta,
    design_matrix  = X_assign)
}

# minimum number of clusters achieving a target power
ind_min_K <- function(T_, m, p0, p1, icc,
                       alpha = 0.05, target_power = 0.80, K_max = 200) {
  # The information matrix is singular for K < T_ (cannot estimate T_ time
  # fixed effects with fewer than T_ clusters), so start the search at K = T_.
  for (K in T_:K_max) {
    res <- sw_power(K, T_, m, p0, p1, icc, alpha)
    if (!is.na(res$power) && res$power >= target_power) return(K)
  }
  NA_integer_
}

# estimate power for each of both scenarios (20 vs 40 surveys)
for (i in seq_along(surveys_per_month)) {

  surveys <- surveys_per_month[i]
  m       <- m_values[i]

  cat(sprintf(">>> SCENARIO: %d surveys / provider / month  (m = %.0f eligible)\n",
              surveys, m))
  cat(sprintf("    (%.0f x %.0f%% antibiotic x %.0f%% appropriate = %.0f eligible/provider/month)\n\n",
              surveys, prop_antibiotic * 100, prop_appropriate * 100, m))

  # Minimum K (any integer)
  min_K_any <- find_min_K(T_, m, p0, p1, icc, alpha, power_target)
  res_min   <- sw_power(min_K_any, T_, m, p0, p1, icc, alpha)

  # Recommended K: smallest multiple of (T_-1) = 11 at or above min_K_any
  # (ensures equal batch sizes at each crossover step)
  step_size  <- T_ - 1L   # = 11
  min_K_rec  <- ceiling(min_K_any / step_size) * step_size
  res_rec    <- sw_power(min_K_rec, T_, m, p0, p1, icc, alpha)

  cat(sprintf("  Minimum K (any integer)          : K = %d  (power = %.1f%%)\n",
              min_K_any, res_min$power * 100))
  cat(sprintf("  Recommended K (multiple of %d)   : K = %d  (power = %.1f%%)\n",
              step_size, min_K_rec, res_rec$power * 100))
  cat(sprintf("  Batch size at each step           : %d provider(s) per month\n",
              min_K_rec / step_size))
  cat(sprintf("  SE(treatment effect)              : %.4f\n", res_rec$se_delta))
  cat(sprintf("  Total eligible patient-months     : %d\n\n",
              as.integer(min_K_rec * T_ * m)))
}

# explore power for different number of clusters/providers (ranging from 5 to 30)
K_seq <- 5:30

for (i in seq_along(surveys_per_month)) {
  surveys <- surveys_per_month[i]
  m       <- m_values[i]

  cat(sprintf("\n--- %d surveys/provider/month  (m = %.0f eligible) ---\n\n",
              surveys, m))
  cat(sprintf("  %-6s  %-10s  %-12s  %-10s\n",
              "K", "Power (%)", "SE(delta)", "Note"))
  cat(sprintf("  %-6s  %-10s  %-12s  %-10s\n",
              "------", "----------", "------------", "----------"))

  for (K in K_seq) {
    res  <- sw_power(K, T_, m, p0, p1, icc, alpha)

    # K < T_ yields a singular information matrix; flag and skip
    if (is.na(res$power)) {
      cat(sprintf("  %-6d  %-10s  %-12s  %-30s\n",
                  K, "n/a", "n/a",
                  paste0("singular (need K >= ", T_, ")")))
      next
    }

    note <- ""
    if (K %% (T_ - 1) == 0) note <- paste0(K / (T_ - 1), " per step")
    prev_power <- sw_power(K - 1, T_, m, p0, p1, icc, alpha)$power
    if (res$power >= power_target && (is.na(prev_power) || prev_power < power_target)) {
      note <- paste(note, "[MIN]")
    }
    cat(sprintf("  %-6d  %-10s  %-12s  %-10s\n",
                K,
                sprintf("%.1f%%", res$power * 100),
                sprintf("%.4f",   res$se_delta),
                trimws(note)))
  }
}


# sensitivity analyses, variation in ICC and in effect size
m_base <- 8   # base case: 40 surveys, m = 8

# ICC
cat("\n--- Varying ICC ---\n\n")
cat(sprintf("  %-8s  %-10s  %-12s  %-12s\n",
            "ICC", "Min K", "K=11 power", "K=22 power"))
cat(sprintf("  %-8s  %-10s  %-12s  %-12s\n",
            "--------", "----------", "------------", "------------"))
for (icc_val in c(0.01, 0.02, 0.05, 0.10, 0.15, 0.20)) {
  mk   <- find_min_K(T_, m_base, p0, p1, icc_val, alpha, power_target)
  p11  <- sw_power(11, T_, m_base, p0, p1, icc_val, alpha)$power
  p22  <- sw_power(22, T_, m_base, p0, p1, icc_val, alpha)$power
  flag <- if (icc_val == 0.05) " <- base case" else ""
  cat(sprintf("  %-8s  %-10d  %-12s  %-12s  %s\n",
              sprintf("%.2f", icc_val), mk,
              sprintf("%.1f%%", p11 * 100),
              sprintf("%.1f%%", p22 * 100),
              flag))
}

# effect size
cat("\n--- Varying assumed effect (p0 = 0.50 fixed) ---\n\n")
cat(sprintf("  %-6s  %-8s  %-10s  %-12s  %-12s\n",
            "p1", "delta", "Min K", "K=11 power", "K=22 power"))
cat(sprintf("  %-6s  %-8s  %-10s  %-12s  %-12s\n",
            "------", "--------", "----------", "------------", "------------"))
for (p1_val in c(0.55, 0.60, 0.65, 0.70, 0.75)) {
  mk  <- find_min_K(T_, m_base, p0, p1_val, icc, alpha, power_target)
  p11 <- sw_power(11, T_, m_base, p0, p1_val, icc, alpha)$power
  p22 <- sw_power(22, T_, m_base, p0, p1_val, icc, alpha)$power
  flag <- if (p1_val == 0.65) " <- base case" else ""
  cat(sprintf("  %-6s  %-8s  %-10s  %-12s  %-12s  %s\n",
              sprintf("%.2f", p1_val),
              sprintf("%.2f", p1_val - p0),
              ifelse(is.na(mk), ">200", as.character(mk)),
              sprintf("%.1f%%", p11 * 100),
              sprintf("%.1f%%", p22 * 100),
              flag))
}

# number of surveys per month (varying m)
cat("\n--- 6c. Varying surveys per provider per month ---\n\n")
cat(sprintf("  %-10s  %-10s  %-10s  %-12s  %-12s\n",
            "Surveys", "m (elig.)", "Min K", "K=11 power", "K=22 power"))
cat(sprintf("  %-10s  %-10s  %-10s  %-12s  %-12s\n",
            "----------", "----------", "----------", "------------", "------------"))
for (surv in c(10, 20, 30, 40, 50, 80)) {
  m_val <- surv * eligibility_rate
  mk    <- find_min_K(T_, m_val, p0, p1, icc, alpha, power_target)
  p11   <- sw_power(11, T_, m_val, p0, p1, icc, alpha)$power
  p22   <- sw_power(22, T_, m_val, p0, p1, icc, alpha)$power
  flag  <- if (surv == 40) " <- base case" else if (surv == 20) " <- scenario 2" else ""
  cat(sprintf("  %-10d  %-10s  %-10s  %-12s  %-12s  %s\n",
              surv,
              sprintf("%.0f", m_val),
              ifelse(is.na(mk), ">200", as.character(mk)),
              sprintf("%.1f%%", p11 * 100),
              sprintf("%.1f%%", p22 * 100),
              flag))
}


# make a plot with the stepped wedge design in a matrix
# install.packages(c("ggplot2", "tidyr", "dplyr", "scales"))
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

build_sw_data <- function(K, T_, surveys) {

  res <- sw_power(K, T_, m = surveys * eligibility_rate, p0, p1, icc, alpha)
  X   <- res$design_matrix   # K x T_

  # Long-format provider x month grid
  df <- expand.grid(provider = seq_len(K), month = seq_len(T_)) |>
    arrange(provider, month)

  df$arm         <- ifelse(X[cbind(df$provider, df$month)] == 0L,
                           "Control", "Intervention")
  df$n_surveyed  <- surveys
  df$n_eligible  <- surveys * eligibility_rate
  df
}

#' including months where one arm has zero providers (M1 intervention, M12 control).
build_monthly_counts <- function(K, T_, surveys) {

  res <- sw_power(K, T_, m = surveys * eligibility_rate, p0, p1, icc, alpha)
  X   <- res$design_matrix   # K x T_

  n_ctrl <- colSums(X == 0L)   # length T_
  n_int  <- colSums(X == 1L)   # length T_

  data.frame(
    month       = rep(seq_len(T_), 2L),
    arm         = rep(c("Control", "Intervention"), each = T_),
    n_providers = c(n_ctrl, n_int),
    n_surveyed  = c(n_ctrl, n_int) * surveys,
    n_eligible  = c(n_ctrl, n_int) * surveys * eligibility_rate,
    stringsAsFactors = FALSE
  )
}

plot_sw_diagram <- function(K, T_, surveys,
                            scenario  = "",
                            save_pdf  = TRUE) {

  m       <- surveys * eligibility_rate
  res     <- sw_power(K, T_, m, p0, p1, icc, alpha)
  df      <- build_sw_data(K, T_, surveys)
  monthly <- build_monthly_counts(K, T_, surveys)
  col_ctrl <- "#D0CFC4"    # light grey-beige  -> control
  col_int  <- "#2E75B6"    # ITM-ish blue      -> intervention
  col_txt  <- "#1A1A1A"
  p_main <- ggplot(df, aes(x = factor(month), y = factor(provider, levels = rev(seq_len(K))), fill = arm)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(arm == "Control", "C", "I")),
              size = 2.8, colour = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Control" = col_ctrl, "Intervention" = col_int), name   = NULL) +
    scale_x_discrete(
      labels = paste0("M", seq_len(T_)), position = "top") +
    scale_y_discrete(
      labels = function(x) paste0("P", formatC(as.integer(x), width = nchar(K), flag = "0"))) +
    labs(x = NULL, y = "Provider") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x      = element_text(size = 9, colour = col_txt),
      axis.text.y      = element_text(size = 9, colour = col_txt),
      axis.title.y     = element_text(size = 10, colour = col_txt, margin = margin(r = 6)),
      panel.grid       = element_blank(),
      legend.position  = "top",
      legend.text      = element_text(size = 10),
      legend.key.size  = unit(0.55, "cm"),
      plot.margin      = margin(4, 8, 0, 8))

strip_rows <- list(
    list(label = "Surveyed (C)",   arm = "Control",      col = "n_surveyed"),
    list(label = "Surveyed (I)",   arm = "Intervention", col = "n_surveyed"),
    list(label = "Eligible (C)",   arm = "Control",      col = "n_eligible"),
    list(label = "Eligible (I)",   arm = "Intervention", col = "n_eligible"))

  surv_ctrl <- monthly$n_surveyed[monthly$arm == "Control"][order(monthly$month[monthly$arm == "Control"])]
  surv_int  <- monthly$n_surveyed[monthly$arm == "Intervention"][order(monthly$month[monthly$arm == "Intervention"])]
  elig_ctrl <- monthly$n_eligible[monthly$arm == "Control"][order(monthly$month[monthly$arm == "Control"])]
  elig_int  <- monthly$n_eligible[monthly$arm == "Intervention"][order(monthly$month[monthly$arm == "Intervention"])]

  make_strip_row <- function(label, arm, vals) {
    data.frame(
      row_label = rep(label, T_ + 1L),
      arm       = rep(arm,   T_ + 1L),
      month     = c(seq_len(T_), T_ + 1L),
      count     = c(vals, sum(vals)),
      is_total  = c(rep(FALSE, T_), TRUE),
      stringsAsFactors = FALSE)  }
  strip_df <- rbind(
    make_strip_row("Surveyed (C)",  "Control",      surv_ctrl),
    make_strip_row("Surveyed (I)",  "Intervention", surv_int),
    make_strip_row("Eligible (C)",  "Control",      elig_ctrl),
    make_strip_row("Eligible (I)",  "Intervention", elig_int) )
  
  # add factor order: top row = Surveyed (C), bottom = Eligible (I)
  strip_df$row_label <- factor(strip_df$row_label,
                                levels = rev(c("Surveyed (C)", "Surveyed (I)",
                                               "Eligible (C)", "Eligible (I)")))
  # totals column gets a darker shade to distinguish it
  col_ctrl_tot <- "#9B9A90"
  col_int_tot  <- "#1A4F7A"

  strip_df <- strip_df |>
    mutate(
      fill_col = case_when(
        arm == "Control"      & !is_total ~ col_ctrl,
        arm == "Control"      &  is_total ~ col_ctrl_tot,
        arm == "Intervention" & !is_total ~ col_int,
        arm == "Intervention" &  is_total ~ col_int_tot))

  # x-axis labels: M1..M12 + "Total"
  x_labels <- c(paste0("M", seq_len(T_)), "Total")
  p_counts <- ggplot(strip_df,
                     aes(x     = factor(month, levels = seq_len(T_ + 1L)),
                         y     = row_label,
                         fill  = fill_col,
                         label = count)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(size = 2.5, colour = "white", fontface = "bold") +
    # vertical separator before totals column
    geom_vline(xintercept = T_ + 0.5, colour = col_txt,
               linewidth = 0.6, linetype = "solid") +
    scale_fill_identity() +
    scale_x_discrete(labels = x_labels) +
    scale_y_discrete() +
    labs(x = "Month", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(size = 8.5, colour = col_txt),
      axis.text.y     = element_text(size = 9,   colour = col_txt, hjust = 1),
      axis.title.x    = element_text(size = 10,  colour = col_txt, margin = margin(t = 4)),
      panel.grid      = element_blank(),
      plot.margin     = margin(0, 8, 4, 8))

  # monthly totals bar chart (surveyed only, stacked C + I)
  monthly_bar <- monthly |>
    select(month, arm, n_surveyed) |>
    mutate(arm = factor(arm, levels = c("Control", "Intervention")))

  p_bar <- ggplot(monthly_bar,
                  aes(x = factor(month), y = n_surveyed, fill = arm)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = surveys * K, linetype = "dashed",
               colour = col_txt, linewidth = 0.4) +
    scale_fill_manual(
      values = c("Control" = col_ctrl, "Intervention" = col_int),
      name   = NULL ) +
    scale_y_continuous(
      name   = "Patients surveyed",
      breaks = pretty_breaks(4),
      expand = expansion(mult = c(0, 0.10)) ) +
    scale_x_discrete(labels = paste0("M", seq_len(T_))) +
    labs(x = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x        = element_blank(),
      axis.text.y        = element_text(size = 9, colour = col_txt),
      axis.title.y       = element_text(size = 10, colour = col_txt, margin = margin(r = 4)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top",
      legend.text        = element_text(size = 10),
      legend.key.size    = unit(0.45, "cm"),
      plot.margin        = margin(4, 8, 0, 8) )

  # assemble all three plots in a single figure
  # install.packages("patchwork") if needed
  library(patchwork)

  # remove legend from p_main since p_bar already carries it
  p_main <- p_main + theme(legend.position = "none")

  subtitle_str <- sprintf(
    paste0("K = %d providers  |  T = %d months  |  %d surveys/provider/month  |  ",
           "m = %.0f eligible/provider/month  |  Power = %.1f%%  |  SE(\u03b4) = %.4f"),
    K, T_, surveys, surveys * eligibility_rate,
    res$power * 100, res$se_delta )
  # Heights: bar ~ fixed, main scales with K, strip ~ fixed
  combined <- p_bar / p_main / p_counts +
    plot_layout(heights = c(1.6, K * 0.30, 2.2)) +
    plot_annotation(
      title    = sprintf("Stepped wedge design — %s", scenario),
      subtitle = subtitle_str,
      caption  = paste0(
        "C = control  |  I = intervention  |  ",
        "Eligible = appropriately dispensed antibiotic (",
        round(eligibility_rate * 100), "% of surveyed)  |  ",
        "Dashed line = total monthly surveys  |  ",
        "Darker shading in ‘Total’ column"
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold",   colour = col_txt),
        plot.subtitle = element_text(size = 9.5,                 colour = "#555555"),
        plot.caption  = element_text(size = 8.5,  colour = "#777777", hjust = 0)))

  if (save_pdf) {
    fname <- sprintf("sw_design_%s.pdf",
                     gsub("[^A-Za-z0-9]", "_", tolower(scenario)))
    h_in  <- 2.5 + K * 0.28 + 1.8   # adaptive height
    ggsave(fname, plot = combined,
           width = 13, height = h_in, units = "in", device = "pdf")
    message("Saved: ", fname)
  }

  combined
}

# make a plot with two scenarios (22 providers and 11 providers per stratum)
steppedwedge_plot_20surveyspermonth_at22providers <- plot_sw_diagram(
  K        = 22,
  T_       = T_,
  surveys  = 20,
  scenario = "20 surveys per provider per month (K = 22)")
steppedwedge_plot_20surveyspermonth_at22providers

steppedwedge_plot_40surveyspermonth_at11providers <- plot_sw_diagram(
  K        = 11,
  T_       = T_,
  surveys  = 40,
  scenario = "40 surveys per provider per month (K = 11)")
steppedwedge_plot_40surveyspermonth_at11providers
