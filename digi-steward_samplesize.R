# Stepped wedge Sample size estimation
# digi-steward
# Primary outcome: proportion with correct antibiotic intake among patients with appropriately dispensed antibiotic treatment

# Parameters
T_ <- 12
p0 <- 0.50; p1 <- 0.65
icc <- 0.05; alpha <- 0.05
target_power <- 0.80
elig_rate <- 0.5 * 0.4

# Stepped-wedge power (Hussey & Hughes)
sw_power <- function(K, T_, m, p0, p1, icc, alpha = 0.05) {
  
  delta <- p1 - p0
  s2 <- ((p0 + p1) / 2) * (1 - (p0 + p1) / 2)
  vb <- s2 * icc
  vw <- s2 * (1 - icc)
  
  Sigma <- matrix(vb, T_, T_)
  diag(Sigma) <- vb + vw / m
  Si <- solve(Sigma)
  
  step <- K / (T_ - 1)
  X <- matrix(0, K, T_)
  for (k in 1:K) {
    sw <- ceiling(k / step)
    if (sw < T_) X[k, (sw + 1):T_] <- 1
  }
  
  Info <- matrix(0, T_ + 1, T_ + 1)
  for (k in 1:K) {
    D <- cbind(diag(T_), X[k, ])
    Info <- Info + t(D) %*% Si %*% D
  }
  
  V <- tryCatch(solve(Info), error = function(e) NULL)
  if (is.null(V)) return(NA)
  
  se <- sqrt(V[T_ + 1, T_ + 1])
  pnorm(abs(delta) / se - qnorm(1 - alpha / 2))
}

find_min_K <- function(T_, m, p0, p1, icc,
                       alpha = 0.05, target = 0.80, Kmax = 200) {
  for (K in T_:Kmax)
    if (sw_power(K, T_, m, p0, p1, icc, alpha) >= target) return(K)
  NA
}

# Results
surveys <- c(20, 40)

for (s in surveys) {
  m <- s * elig_rate
  Kmin <- find_min_K(T_, m, p0, p1, icc, alpha, target_power)
  Krec <- ceiling(Kmin / (T_ - 1)) * (T_ - 1)
  
  cat("\n", s, "surveys/provider/month\n")
  cat("Eligible per month:", m, "\n")
  cat("Minimum K:", Kmin, "\n")
  cat("Recommended K:", Krec, "\n")
  cat("Providers per step:", Krec / (T_ - 1), "\n")
}
                
# draw plot with design and number of surveys
# install.packages(c("ggplot2", "tidyr", "dplyr", "scales"))
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Build a tidy data frame describing the SW design + patient counts
#' @param K          Number of providers
#' @param T_         Number of time periods
#' @param surveys    Total patients surveyed per provider per month
#' @return A data.frame with one row per provider x month combination
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

#' Build monthly aggregate counts (for annotation strip at bottom of plot)
build_monthly_counts <- function(K, T_, surveys) {

  df <- build_sw_data(K, T_, surveys)
  df |>
    group_by(month, arm) |>
    summarise(
      n_providers = n(),
      n_surveyed  = sum(n_surveyed),
      n_eligible  = sum(n_eligible),
      .groups     = "drop")
}

#' Draw the full stepped wedge diagram with ggplot2
plot_sw_diagram <- function(K, T_, surveys,
                            scenario  = "",
                            save_pdf  = TRUE) {

  m       <- surveys * eligibility_rate
  res     <- sw_power(K, T_, m, p0, p1, icc, alpha)
  df      <- build_sw_data(K, T_, surveys)
  monthly <- build_monthly_counts(K, T_, surveys)

  # ---- Colour palette ----
  col_ctrl <- "#D0CFC4"    # light grey-beige  -> control
  col_int  <- "#2E75B6"    # ITM-ish blue      -> intervention
  col_txt  <- "#1A1A1A"

  # ---- 1. Main design matrix (provider x month) ----
  p_main <- ggplot(df, aes(x = factor(month), y = factor(provider, levels = rev(seq_len(K))),
                           fill = arm)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(arm == "Control", "C", "I")),
              size = 2.8, colour = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Control" = col_ctrl, "Intervention" = col_int),
      name   = NULL
    ) +
    scale_x_discrete(
      labels = paste0("M", seq_len(T_)),
      position = "top"
    ) +
    scale_y_discrete(
      labels = function(x) paste0("P", formatC(as.integer(x),
                                               width = nchar(K), flag = "0"))
    ) +
    labs(x = NULL, y = "Provider") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x      = element_text(size = 9, colour = col_txt),
      axis.text.y      = element_text(size = 9, colour = col_txt),
      axis.title.y     = element_text(size = 10, colour = col_txt, margin = margin(r = 6)),
      panel.grid       = element_blank(),
      legend.position  = "top",
      legend.text      = element_text(size = 10),
      legend.key.size  = unit(0.55, "cm"),
      plot.margin      = margin(4, 8, 0, 8)
    )

  # ---- 2. Count strip: 4 rows (Surveyed C/I, Eligible C/I) x T_ months ----
  # Build one row per count type x arm, with a "Total" pseudo-month appended.
  # month = T_ + 1 is used for the totals column.

  strip_rows <- list(
    list(label = "Surveyed (C)",   arm = "Control",      col = "n_surveyed"),
    list(label = "Surveyed (I)",   arm = "Intervention", col = "n_surveyed"),
    list(label = "Eligible (C)",   arm = "Control",      col = "n_eligible"),
    list(label = "Eligible (I)",   arm = "Intervention", col = "n_eligible")
  )

  strip_df <- lapply(strip_rows, function(r) {
    vals <- monthly |>
      filter(arm == r$arm) |>
      arrange(month) |>
      pull(r$col)
    data.frame(
      row_label = r$label,
      arm       = r$arm,
      month     = c(seq_len(T_), T_ + 1L),        # months 1..T_ + totals col
      count     = c(vals, sum(vals)),
      is_total  = c(rep(FALSE, T_), TRUE),
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()

  # Factor order: top row = Surveyed (C), bottom = Eligible (I)
  strip_df$row_label <- factor(strip_df$row_label,
                                levels = rev(c("Surveyed (C)", "Surveyed (I)",
                                               "Eligible (C)", "Eligible (I)")))

  # Totals column gets a darker shade to distinguish it
  col_ctrl_tot <- "#9B9A90"
  col_int_tot  <- "#1A4F7A"

  strip_df <- strip_df |>
    mutate(
      fill_col = case_when(
        arm == "Control"      & !is_total ~ col_ctrl,
        arm == "Control"      &  is_total ~ col_ctrl_tot,
        arm == "Intervention" & !is_total ~ col_int,
        arm == "Intervention" &  is_total ~ col_int_tot
      )
    )

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
      plot.margin     = margin(0, 8, 4, 8)
    )

  # ---- 3. Monthly totals bar chart (surveyed only, stacked C + I) ----
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
      name   = NULL
    ) +
    scale_y_continuous(
      name   = "Patients surveyed",
      breaks = pretty_breaks(4),
      expand = expansion(mult = c(0, 0.10))
    ) +
    scale_x_discrete(labels = paste0("M", seq_len(T_))) +
    labs(x = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x        = element_blank(),
      axis.text.y        = element_text(size = 9, colour = col_txt),
      axis.title.y       = element_text(size = 10, colour = col_txt, margin = margin(r = 4)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top",
      legend.text        = element_text(size = 10),
      legend.key.size    = unit(0.45, "cm"),
      plot.margin        = margin(4, 8, 0, 8)
    )

  # ---- Assemble with patchwork ----
  # install.packages("patchwork") if needed
  library(patchwork)

  # Remove legend from p_main since p_bar already carries it
  p_main <- p_main + theme(legend.position = "none")

  subtitle_str <- sprintf(
    paste0("K = %d providers  |  T = %d months  |  %d surveys/provider/month  |  ",
           "m = %.0f eligible/provider/month  |  Power = %.1f%%  |  SE(\u03b4) = %.4f"),
    K, T_, surveys, surveys * eligibility_rate,
    res$power * 100, res$se_delta
  )

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
        plot.caption  = element_text(size = 8.5,  colour = "#777777", hjust = 0)
      )
    )

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

# ---- Draw and save both scenarios ----
plot_sw_diagram(
  K        = 22,
  T_       = T_,
  surveys  = 20,
  scenario = "20 surveys per provider per month (K = 22)"
)

plot_sw_diagram(
  K        = 11,
  T_       = T_,
  surveys  = 40,
  scenario = "40 surveys per provider per month (K = 11)"
)
