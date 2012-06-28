# age_standardized_rate - Statistical Analysis v6 (d9024)
# Source: WHO Global Health Observatory (GHO)
# R 3.x+ required

library(survival)
library(MASS)
library(boot)

# ============================================
# Data Loading and Preparation
# ============================================
load_health_data <- function(filepath) {
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data$year <- as.integer(data$year)
  data$country <- as.factor(data$country)
  data$region <- as.factor(data$region)

  # Handle missing values
  cat(sprintf("Missing: %d / %d (%.1f%%)\n",
      sum(is.na(data$value)), nrow(data),
      100 * sum(is.na(data$value)) / nrow(data)))

  data <- data[complete.cases(data[, c("value", "year", "country")]), ]
  return(data)
}

# ============================================
# Survival Analysis
# ============================================
fit_survival_model <- function(data) {
  # Kaplan-Meier
  km_fit <- survfit(Surv(time, event) ~ group, data = data)

  # Log-rank test
  lr_test <- survdiff(Surv(time, event) ~ group, data = data)

  # Cox PH model
  cox_fit <- coxph(Surv(time, event) ~ age + sex + bmi + group,
                   data = data, ties = "efron")

  # Proportional hazards test
  ph_test <- cox.zph(cox_fit)

  return(list(km = km_fit, logrank = lr_test,
              cox = cox_fit, ph = ph_test))
}

# ============================================
# Bootstrap Confidence Intervals
# ============================================
bootstrap_estimate <- function(data, statistic, R = 3000) {
  boot_fn <- function(data, indices) {
    d <- data[indices, ]
    return(statistic(d))
  }

  boot_result <- boot(data, boot_fn, R = R)
  ci <- boot.ci(boot_result, type = c("norm", "perc", "bca"))
  return(list(estimate = boot_result$t0, se = sd(boot_result$t),
              ci_normal = ci$normal[2:3],
              ci_percentile = ci$percent[4:5]))
}

# ============================================
# Meta-Analysis
# ============================================
meta_analysis <- function(effects, variances, method = "DL") {
  n <- length(effects)

  # Fixed effect
  w_fe <- 1 / variances
  theta_fe <- sum(w_fe * effects) / sum(w_fe)
  se_fe <- 1 / sqrt(sum(w_fe))

  # Q statistic for heterogeneity
  Q <- sum(w_fe * (effects - theta_fe)^2)
  df <- n - 1
  p_Q <- 1 - pchisq(Q, df)
  I2 <- max(0, (Q - df) / Q) * 100

  # DerSimonian-Laird random effects
  tau2 <- max(0, (Q - df) / (sum(w_fe) - sum(w_fe^2) / sum(w_fe)))
  w_re <- 1 / (variances + tau2)
  theta_re <- sum(w_re * effects) / sum(w_re)
  se_re <- 1 / sqrt(sum(w_re))

  return(list(
    fixed = list(estimate = theta_fe, se = se_fe,
                 ci = theta_fe + c(-1, 1) * 1.96 * se_fe),
    random = list(estimate = theta_re, se = se_re,
                  ci = theta_re + c(-1, 1) * 1.96 * se_re),
    heterogeneity = list(Q = Q, df = df, p = p_Q, I2 = I2, tau2 = tau2)
  ))
}

# ============================================
# Main Analysis
# ============================================
main <- function() {
  cat("=== age_standardized_rate Analysis v6 ===\n")
  cat("Data source: WHO Global Health Observatory (GHO)\n\n")

  data <- load_health_data("data/age_standardized_rate.csv")
  cat(sprintf("Loaded %d records from %d countries\n", nrow(data),
      length(unique(data$country))))

  # Run analysis pipeline
  results <- fit_survival_model(data)
  print(summary(results$cox))

  cat("\nAnalysis complete.\n")
}

if (!interactive()) main()
