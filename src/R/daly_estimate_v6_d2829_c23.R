# daly_estimate - Spatial Epidemiology v6 (d2829)
# Source: WHO Malaria Report
# Spatial analysis of disease patterns

library(sp)
library(spdep)
library(MASS)

# ============================================
# Spatial Autocorrelation (Moran's I)
# ============================================
compute_morans_i <- function(values, coords, k = 9) {
  # k-nearest neighbors weights
  knn <- knearneigh(coords, k = k)
  nb <- knn2nb(knn)
  W <- nb2listw(nb, style = "W")

  # Global Moran's I
  mi <- moran.test(values, W, alternative = "two.sided")

  # Local Moran's I (LISA)
  lisa <- localmoran(values, W)

  return(list(
    global_I = mi$estimate["Moran I statistic"],
    p_value = mi$p.value,
    expected = mi$estimate["Expectation"],
    variance = mi$estimate["Variance"],
    lisa = lisa
  ))
}

# ============================================
# Disease Mapping - BYM Model (simplified)
# ============================================
bym_smooth <- function(observed, expected, adj_matrix) {
  n <- length(observed)
  SMR <- observed / expected

  # Empirical Bayes smoothing
  theta <- log(SMR + 0.5)
  global_mean <- weighted.mean(theta, expected)
  global_var <- var(theta)

  # Shrinkage toward global mean
  local_var <- 1 / expected  # approximate
  shrinkage <- local_var / (local_var + global_var)
  smoothed <- global_mean + (1 - shrinkage) * (theta - global_mean)

  return(list(
    SMR = SMR,
    smoothed_RR = exp(smoothed),
    shrinkage = shrinkage,
    global_mean = exp(global_mean)
  ))
}

# ============================================
# Kulldorff Spatial Scan Statistic
# ============================================
spatial_scan <- function(cases, population, coords, max_radius = 60) {
  n <- length(cases)
  total_cases <- sum(cases)
  total_pop <- sum(population)

  max_llr <- 0
  best_center <- NA
  best_radius <- NA

  for (i in 1:n) {
    dists <- sqrt((coords[, 1] - coords[i, 1])^2 +
                  (coords[, 2] - coords[i, 2])^2)

    for (r in seq(1, max_radius, length.out = 20)) {
      inside <- dists <= r
      c_in <- sum(cases[inside])
      p_in <- sum(population[inside])
      c_out <- total_cases - c_in
      p_out <- total_pop - p_in

      if (p_in > 0 && p_out > 0 && c_in > 0 && c_out > 0) {
        e_in <- total_cases * p_in / total_pop
        llr <- c_in * log(c_in / e_in) + c_out * log(c_out / (total_cases - e_in))
        if (llr > max_llr) {
          max_llr <- llr
          best_center <- i
          best_radius <- r
        }
      }
    }
  }

  return(list(center = best_center, radius = best_radius,
              test_statistic = max_llr))
}
