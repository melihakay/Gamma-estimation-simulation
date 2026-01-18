# ============================================================
# STAT 303 Term Project - Figure Generation Script
# Point Estimation of Gamma Distribution Parameters
# Authors: Rubar Akyıldız - Melih Akay
# ============================================================

# Setting seed for reproducibility
set.seed(303)

# Create figures directory if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# ============================================================
# Figure 1: Histogram of Simulated Data (Section 4.1)
# ============================================================

# True parameters for descriptive statistics
true_k <- 3
true_theta <- 2
n <- 100

# Generate data
data_sample <- rgamma(n, shape = true_k, scale = true_theta)

# Calculate sample statistics
sample_mean <- mean(data_sample)
sample_var <- var(data_sample)

cat("=== Descriptive Statistics ===\n")
cat("Sample Mean (X_bar):", round(sample_mean, 4), "\n")
cat("Sample Variance (S^2):", round(sample_var, 4), "\n")

# Save histogram
pdf("figures/histogram_descriptive.pdf", width = 7, height = 5)
hist(data_sample, breaks = 15, probability = TRUE,
     main = "Histogram of Simulated Gamma Data",
     xlab = "Data Values", col = "lightblue", border = "white")
lines(density(data_sample), col = "darkblue", lwd = 2)
dev.off()

cat("Figure 1 saved: figures/histogram_descriptive.pdf\n")

# ============================================================
# Figure 2: Fitted Gamma Densities (Section 4.3)
# ============================================================

# Method of Moments Estimation
k_mom <- (mean(data_sample)^2) / var(data_sample)
theta_mom <- var(data_sample) / mean(data_sample)

# Maximum Likelihood Estimation
target_val <- log(mean(data_sample)) - mean(log(data_sample))

mle_eqn <- function(k) {
  log(k) - digamma(k) - target_val
}

mle_root <- uniroot(mle_eqn, interval = c(0.1, 20), extendInt = "yes")
k_mle <- mle_root$root
theta_mle <- mean(data_sample) / k_mle

cat("\n=== Parameter Estimates ===\n")
cat("True values: k =", true_k, ", theta =", true_theta, "\n")
cat("MoM: k =", round(k_mom, 4), ", theta =", round(theta_mom, 4), "\n")
cat("MLE: k =", round(k_mle, 4), ", theta =", round(theta_mle, 4), "\n")

# Save fitted densities plot
pdf("figures/fitted_densities.pdf", width = 7, height = 5)
hist(data_sample, probability = TRUE, breaks = 15,
     main = "Fitted Gamma Densities over Histogram",
     xlab = "x", col = "lightgray", border = "white",
     ylim = c(0, max(density(data_sample)$y) * 1.2))

x_vals <- seq(min(data_sample), max(data_sample), length.out = 100)

lines(x_vals, dgamma(x_vals, shape = k_mom, scale = theta_mom),
      col = "blue", lwd = 2, lty = 2)
lines(x_vals, dgamma(x_vals, shape = k_mle, scale = theta_mle),
      col = "red", lwd = 2)

legend("topright", legend = c("MoM Fit", "MLE Fit"),
       col = c("blue", "red"), lty = c(2, 1), lwd = 2)
dev.off()

cat("Figure 2 saved: figures/fitted_densities.pdf\n")

# ============================================================
# Simulation Study (Section 5)
# ============================================================

cat("\n=== Running Simulation Study ===\n")
cat("This may take a moment...\n")

R <- 2000     # Number of replications
sample_sizes <- c(20, 50, 100)

scenarios <- list(
  "Scenario 1" = c(k = 1, theta = 2),  # High Skewness
  "Scenario 2" = c(k = 5, theta = 1)   # Moderate Skewness
)

results_df <- data.frame()

for (scen_name in names(scenarios)) {

  true_k <- scenarios[[scen_name]]["k"]
  true_theta <- scenarios[[scen_name]]["theta"]

  for (n in sample_sizes) {

    k_mom_est <- numeric(R); t_mom_est <- numeric(R)
    k_mle_est <- numeric(R); t_mle_est <- numeric(R)

    for (i in 1:R) {
      # Generate Data
      x <- rgamma(n, shape = true_k, scale = true_theta)

      # MoM Estimation
      x_bar <- mean(x)
      s2 <- var(x)

      k_hat_mom <- (x_bar^2) / s2
      t_hat_mom <- s2 / x_bar

      k_mom_est[i] <- k_hat_mom
      t_mom_est[i] <- t_hat_mom

      # MLE Estimation
      target <- log(x_bar) - mean(log(x))

      try({
        mle_sol <- uniroot(function(k) log(k) - digamma(k) - target,
                           interval = c(1e-5, 100), extendInt = "yes")
        k_hat_mle <- mle_sol$root
        t_hat_mle <- x_bar / k_hat_mle

        k_mle_est[i] <- k_hat_mle
        t_mle_est[i] <- t_hat_mle
      }, silent = TRUE)
    }

    # Performance Metrics Calculation
    calc_metrics <- function(est, true_val) {
      bias <- mean(est, na.rm=TRUE) - true_val
      variance <- var(est, na.rm=TRUE)
      mse <- variance + bias^2
      return(c(bias, variance, mse))
    }

    mom_k <- calc_metrics(k_mom_est, true_k)
    mle_k <- calc_metrics(k_mle_est, true_k)
    mom_t <- calc_metrics(t_mom_est, true_theta)
    mle_t <- calc_metrics(t_mle_est, true_theta)

    temp_res <- data.frame(
      Scenario = scen_name,
      n = n,
      Method = rep(c("MoM", "MLE"), 2),
      Parameter = rep(c("k", "theta"), each = 2),
      Bias = c(mom_k[1], mle_k[1], mom_t[1], mle_t[1]),
      Variance = c(mom_k[2], mle_k[2], mom_t[2], mle_t[2]),
      MSE = c(mom_k[3], mle_k[3], mom_t[3], mle_t[3])
    )
    results_df <- rbind(results_df, temp_res)
  }
  cat("Completed:", scen_name, "\n")
}

# Save results to CSV for LaTeX table
write.csv(results_df, "figures/simulation_results.csv", row.names = FALSE)
cat("Simulation results saved: figures/simulation_results.csv\n")

# ============================================================
# Figure 3: MSE Plots for Scenario 1 (High Skewness)
# ============================================================

pdf("figures/mse_scenario1.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Plot 1: MSE of k (Scenario 1)
plot_data <- subset(results_df, Scenario == "Scenario 1" & Parameter == "k")
mom_data <- subset(plot_data, Method == "MoM")
mle_data <- subset(plot_data, Method == "MLE")

plot(mom_data$n, mom_data$MSE, type = "b", col = "blue", pch = 19, lty = 2,
     ylim = c(0, max(plot_data$MSE)*1.1),
     xlab = "Sample Size (n)", ylab = "MSE",
     main = expression(paste("MSE of ", hat(k), " (High Skewness, k=1)")))
lines(mle_data$n, mle_data$MSE, type = "b", col = "red", pch = 19, lty = 1)
legend("topright", legend = c("MoM", "MLE"), col = c("blue", "red"), lty = c(2, 1), pch = 19)

# Plot 2: MSE of theta (Scenario 1)
plot_data_t <- subset(results_df, Scenario == "Scenario 1" & Parameter == "theta")
mom_data_t <- subset(plot_data_t, Method == "MoM")
mle_data_t <- subset(plot_data_t, Method == "MLE")

plot(mom_data_t$n, mom_data_t$MSE, type = "b", col = "blue", pch = 19, lty = 2,
     ylim = c(0, max(plot_data_t$MSE)*1.1),
     xlab = "Sample Size (n)", ylab = "MSE",
     main = expression(paste("MSE of ", hat(theta), " (High Skewness, k=1)")))
lines(mle_data_t$n, mle_data_t$MSE, type = "b", col = "red", pch = 19, lty = 1)
legend("topright", legend = c("MoM", "MLE"), col = c("blue", "red"), lty = c(2, 1), pch = 19)

dev.off()
cat("Figure 3 saved: figures/mse_scenario1.pdf\n")

# ============================================================
# Figure 4: MSE Plots for Scenario 2 (Moderate Skewness)
# ============================================================

pdf("figures/mse_scenario2.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Plot 1: MSE of k (Scenario 2)
plot_data <- subset(results_df, Scenario == "Scenario 2" & Parameter == "k")
mom_data <- subset(plot_data, Method == "MoM")
mle_data <- subset(plot_data, Method == "MLE")

plot(mom_data$n, mom_data$MSE, type = "b", col = "blue", pch = 19, lty = 2,
     ylim = c(0, max(plot_data$MSE)*1.1),
     xlab = "Sample Size (n)", ylab = "MSE",
     main = expression(paste("MSE of ", hat(k), " (Moderate Skewness, k=5)")))
lines(mle_data$n, mle_data$MSE, type = "b", col = "red", pch = 19, lty = 1)
legend("topright", legend = c("MoM", "MLE"), col = c("blue", "red"), lty = c(2, 1), pch = 19)

# Plot 2: MSE of theta (Scenario 2)
plot_data_t <- subset(results_df, Scenario == "Scenario 2" & Parameter == "theta")
mom_data_t <- subset(plot_data_t, Method == "MoM")
mle_data_t <- subset(plot_data_t, Method == "MLE")

plot(mom_data_t$n, mom_data_t$MSE, type = "b", col = "blue", pch = 19, lty = 2,
     ylim = c(0, max(plot_data_t$MSE)*1.1),
     xlab = "Sample Size (n)", ylab = "MSE",
     main = expression(paste("MSE of ", hat(theta), " (Moderate Skewness, k=5)")))
lines(mle_data_t$n, mle_data_t$MSE, type = "b", col = "red", pch = 19, lty = 1)
legend("topright", legend = c("MoM", "MLE"), col = c("blue", "red"), lty = c(2, 1), pch = 19)

dev.off()
cat("Figure 4 saved: figures/mse_scenario2.pdf\n")

# ============================================================
# Print Summary Table
# ============================================================

cat("\n=== Simulation Results Summary ===\n")
print(results_df, digits = 4)

cat("\n=== All figures generated successfully! ===\n")
