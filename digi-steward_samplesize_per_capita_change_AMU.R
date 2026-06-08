####################################################################
# digi-steward                                                     #
# Minimum Detectable Effect for per-capita antibiotic use, SRS     #
####################################################################

# assumptions
#   - Antibiotic use follows a Poisson process (courses/person/year)
#   - Two visit streams: public health centres (PHC) and pharmacies/other vendors
#   - PHC share: 73%; other vendors: 27%
#   - Antibiotic probability per PHC visit: 50%; per vendor visit: 25%
#   - One year of follow-up; no loss to follow-up
#   - Simple random sample (no clustering, no design effect)
#   - No overdispersion
#   - alpha = 0.05 (two-sided), power = 80%

# parameters
alpha  <- 0.05
power  <- 0.80
z_a    <- qnorm(1 - alpha / 2)   # 1.960
z_b    <- qnorm(power)            # 0.842

phc_share    <- 0.73;  vendor_share <- 0.27
phc_rx_prob  <- 0.50;  vendor_rx_prob <- 0.25

scenarios <- data.frame(
  country = c("Burkina Faso", "South Africa"),
  visits   = c(2.00, 3.58),
  n        = c(500, 500))

# baseline lambda: weighted antibiotic rate per visit stream
scenarios$lambda0 <- with(scenarios,
  visits * phc_share * phc_rx_prob +  visits * vendor_share * vendor_rx_prob)

# minimum detectable effect numerically (exact Poisson power equation)
# power condition: sqrt(n) * (lambda0 - lambda1) = z_a*sqrt(lambda0) + z_b*sqrt(lambda1)
# tearranged as root of: f(lambda1) = 0

mde_solve <- function(lambda0, n) {
  f <- function(l1) sqrt(n) * (lambda0 - l1) - z_a * sqrt(lambda0) - z_b * sqrt(l1)
  uniroot(f, interval = c(1e-6, lambda0 - 1e-6))$root}

scenarios$lambda1     <- mapply(mde_solve, scenarios$lambda0, scenarios$n)
scenarios$abs_mde     <- scenarios$lambda0 - scenarios$lambda1
scenarios$rel_mde_pct <- round(scenarios$abs_mde / scenarios$lambda0 * 100, 1)

# output
cat("\n=== minimum detectable effect for per-capita antibiotic use (alpha=0.05, power=80%, n=500) ===\n\n")
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  cat(sprintf(
    "%s (%.2f visits/capita/year)\n  lambda0 = %.3f  |  detectable lambda1 = %.3f\n  MDE: %.3f courses/person/year  =>  %.1f%% relative change\n\n",
    s$country, s$visits, s$lambda0, s$lambda1, s$abs_mde, s$rel_mde_pct))}
