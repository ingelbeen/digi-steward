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